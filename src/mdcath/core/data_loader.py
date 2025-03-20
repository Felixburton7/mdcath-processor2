#!/usr/bin/env python3
"""
Core functionality for loading and processing H5 data from mdCATH dataset.
"""

import os
import h5py
import logging
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union, Any
from concurrent.futures import ProcessPoolExecutor, as_completed

class H5DataLoader:
    """
    Class for efficiently loading and extracting data from mdCATH H5 files.
    Uses chunking/streaming to handle large files.
    """

    def __init__(self, h5_path: str, config: Dict[str, Any]):
        """
        Initialize the H5 data loader.

        Args:
            h5_path: Path to H5 file
            config: Configuration dictionary
        """
        self.h5_path = h5_path
        self.config = config
        self.domain_id = os.path.basename(h5_path).replace("mdcath_dataset_", "").replace(".h5", "")
        self._validate_h5()

    def _validate_h5(self) -> bool:
        """
        Validate that the H5 file has the expected structure.
        
        Returns:
            Boolean indicating if the file is valid
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                # Check if domain exists
                if self.domain_id not in f:
                    logging.error(f"Domain {self.domain_id} not found in {self.h5_path}")
                    return False
                
                # Check for required metadata fields
                required_metadata = ["resid", "resname"]
                for field in required_metadata:
                    if field not in f[self.domain_id]:
                        logging.error(f"Required metadata field '{field}' not found for domain {self.domain_id}")
                        return False
                
                # Check for required temperature groups
                temps = [str(t) for t in self.config.get("temperatures", [320, 348, 379, 413, 450])]
                num_replicas = self.config.get("num_replicas", 5)
                
                temp_found = False
                for temp in temps:
                    if temp in f[self.domain_id]:
                        temp_found = True
                        # Check for replica groups
                        for r in range(num_replicas):
                            replica = str(r)
                            if replica in f[self.domain_id][temp]:
                                # Check for specific datasets
                                required_datasets = ['rmsf', 'dssp', 'coords']
                                for dataset in required_datasets:
                                    if dataset not in f[self.domain_id][temp][replica]:
                                        logging.warning(f"Dataset {dataset} not found for temperature {temp}, " 
                                                    f"replica {replica} in domain {self.domain_id}")
                            else:
                                logging.warning(f"Replica {replica} not found for temperature {temp} in domain {self.domain_id}")
                    else:
                        logging.warning(f"Temperature {temp} not found for domain {self.domain_id}")
                
                if not temp_found:
                    logging.error(f"No valid temperature groups found for domain {self.domain_id}")
                    return False
                    
                return True
        except Exception as e:
            logging.error(f"Failed to validate H5 file {self.h5_path}: {e}")
            return False

    def extract_rmsf(self, temperature: str, replica: str) -> Optional[pd.DataFrame]:
        """
        Extract RMSF data for a specific temperature and replica.
        
        Args:
            temperature: Temperature (e.g., "320")
            replica: Replica (e.g., "0")
        
        Returns:
            DataFrame with RMSF data or None if extraction fails
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                # Check if temperature and replica exist
                if temperature not in f[self.domain_id] or replica not in f[self.domain_id][temperature]:
                    logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
                    return None
                
                # Check if RMSF dataset exists
                if 'rmsf' not in f[self.domain_id][temperature][replica]:
                    logging.warning(f"RMSF data not found for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None
                
                # Get dataset and check size
                rmsf_dataset = f[self.domain_id][temperature][replica]['rmsf']
                rmsf_data = rmsf_dataset[:]
                
                # Extract residue information
                resid_dataset = f[self.domain_id]['resid']
                resname_dataset = f[self.domain_id]['resname']
                
                resids = resid_dataset[:]
                resnames_bytes = resname_dataset[:]
                
                # Decode resnames from bytes to strings
                resnames = [name.decode('utf-8') if isinstance(name, bytes) else str(name) for name in resnames_bytes]
                
                # Handle dimension mismatch - RMSF is typically per residue (not per atom)
                if len(resids) != len(rmsf_data):
                    logging.info(f"Dimension mismatch: resids {len(resids)}, rmsf_data {len(rmsf_data)}")
                    
                    # Build a proper mapping between atoms and residues
                    residue_dict = {}
                    for i, resid in enumerate(resids):
                        if resid not in residue_dict:
                            residue_dict[resid] = resnames[i]
                    
                    # Get the unique residue IDs in order
                    unique_resids = sorted(residue_dict.keys())
                    unique_resnames = [residue_dict[resid] for resid in unique_resids]
                    
                    if len(unique_resids) == len(rmsf_data):
                        logging.info(f"Using unique residue IDs for RMSF data alignment")
                        # Use unique residues for RMSF data
                        resids = unique_resids
                        resnames = unique_resnames
                    elif len(unique_resids) > len(rmsf_data):
                        # If we have more unique residues than RMSF data points,
                        # match by position in sequence
                        logging.warning(f"More unique residues ({len(unique_resids)}) than RMSF data points ({len(rmsf_data)})")
                        resids = unique_resids[:len(rmsf_data)]
                        resnames = unique_resnames[:len(rmsf_data)]
                    elif len(unique_resids) < len(rmsf_data):
                        # If we have fewer unique residues than RMSF data points,
                        # truncate the RMSF data
                        logging.warning(f"Fewer unique residues ({len(unique_resids)}) than RMSF data points ({len(rmsf_data)})")
                        rmsf_data = rmsf_data[:len(unique_resids)]
                    
                # Create DataFrame
                df = pd.DataFrame({
                    'domain_id': self.domain_id,
                    'resid': resids,
                    'resname': resnames,
                    f'rmsf_{temperature}': rmsf_data
                })
                
                return df
        except Exception as e:
            logging.error(f"Failed to extract RMSF data: {e}")
            return None

    def extract_pdb(self) -> Optional[str]:
        """
        Extract PDB data from the H5 file.

        Returns:
            PDB string or None if extraction fails
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                pdb_data = f[self.domain_id]['pdb'][()]
                if isinstance(pdb_data, bytes):
                    return pdb_data.decode('utf-8')
                return str(pdb_data)
        except Exception as e:
            logging.error(f"Failed to extract PDB data: {e}")
            return None

    def extract_dssp(self, temperature: str, replica: str, frame: int = -1) -> Optional[pd.DataFrame]:
        """
        Extract DSSP data for a specific temperature, replica, and frame.

        Args:
            temperature: Temperature (e.g., "320")
            replica: Replica (e.g., "0")
            frame: Frame index (default: -1 for last frame)

        Returns:
            DataFrame with DSSP data or None if extraction fails
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                # Check if temperature and replica exist
                if temperature not in f[self.domain_id] or replica not in f[self.domain_id][temperature]:
                    logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
                    return None

                # Check if DSSP dataset exists
                if 'dssp' not in f[self.domain_id][temperature][replica]:
                    logging.warning(f"DSSP data not found for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None
                    
                # Extract DSSP data
                dssp_dataset = f[self.domain_id][temperature][replica]['dssp']
                
                # Handle the case where frame is out of bounds
                num_frames = dssp_dataset.shape[0] if len(dssp_dataset.shape) > 0 else 0
                if num_frames == 0:
                    logging.warning(f"Empty DSSP dataset for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None
                    
                if frame < 0:
                    # Convert negative indices to positive
                    frame = num_frames + frame
                    
                if frame < 0 or frame >= num_frames:
                    logging.warning(f"Frame index {frame} out of bounds (0-{num_frames-1}) for domain {self.domain_id}")
                    frame = min(max(0, frame), num_frames - 1)  # Clamp to valid range
                    
                dssp_data = dssp_dataset[frame]

                # Extract residue information
                resids = f[self.domain_id]['resid'][:]
                resnames = [name.decode('utf-8') if isinstance(name, bytes) else str(name) for name in f[self.domain_id]['resname'][:]]

                # Decode DSSP codes
                dssp_codes = [code.decode('utf-8') if isinstance(code, bytes) else str(code) for code in dssp_data]

                # Handle dimension mismatch - DSSP is per residue, not per atom
                if len(resids) != len(dssp_codes):
                    logging.info(f"Dimension mismatch in DSSP: resids {len(resids)}, dssp_codes {len(dssp_codes)}")
                    
                    # Build a proper mapping between atoms and residues
                    residue_dict = {}
                    for i, resid in enumerate(resids):
                        if resid not in residue_dict:
                            residue_dict[resid] = resnames[i]
                    
                    # Get the unique residue IDs in order
                    unique_resids = sorted(residue_dict.keys())
                    unique_resnames = [residue_dict[resid] for resid in unique_resids]
                    
                    if len(unique_resids) == len(dssp_codes):
                        logging.info(f"Using unique residue IDs for DSSP data alignment")
                        # Use unique residues for DSSP data
                        resids = unique_resids
                        resnames = unique_resnames
                    elif len(unique_resids) > len(dssp_codes):
                        # If we have more unique residues than DSSP data points, 
                        # match by position in sequence
                        logging.warning(f"More unique residues ({len(unique_resids)}) than DSSP data points ({len(dssp_codes)})")
                        resids = unique_resids[:len(dssp_codes)]
                        resnames = unique_resnames[:len(dssp_codes)]
                    elif len(unique_resids) < len(dssp_codes):
                        # If we have fewer unique residues than DSSP data points,
                        # truncate the DSSP data
                        logging.warning(f"Fewer unique residues ({len(unique_resids)}) than DSSP data points ({len(dssp_codes)})")
                        dssp_codes = dssp_codes[:len(unique_resids)]

                # Create DataFrame
                df = pd.DataFrame({
                    'domain_id': self.domain_id,
                    'resid': resids,
                    'resname': resnames,
                    'dssp': dssp_codes
                })

                return df
        except Exception as e:
            logging.error(f"Failed to extract DSSP data: {e}")
            return None

    def extract_coordinates(self, temperature: str, replica: str, frame: int = -1) -> Optional[Tuple[np.ndarray, List[int], List[str]]]:
        """
        Extract coordinate data for a specific temperature, replica, and frame.

        Args:
            temperature: Temperature (e.g., "320")
            replica: Replica (e.g., "0")
            frame: Frame index (default: -1 for last frame)

        Returns:
            Tuple of (coordinates, residue IDs, residue names) or None if extraction fails
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                # Check if temperature and replica exist
                if temperature not in f[self.domain_id] or replica not in f[self.domain_id][temperature]:
                    logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
                    return None

                # Extract coordinate data for specified frame
                if 'coords' not in f[self.domain_id][temperature][replica]:
                    logging.warning(f"Coordinate data not found for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None
                    
                coords_dataset = f[self.domain_id][temperature][replica]['coords']
                
                # Handle the case where frame is out of bounds
                num_frames = coords_dataset.shape[0] if len(coords_dataset.shape) > 0 else 0
                if num_frames == 0:
                    logging.warning(f"Empty coordinates dataset for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None
                
                if frame < 0:
                    # Convert negative indices to positive
                    frame = num_frames + frame
                    
                if frame < 0 or frame >= num_frames:
                    logging.warning(f"Frame index {frame} out of bounds (0-{num_frames-1}) for domain {self.domain_id}")
                    frame = min(max(0, frame), num_frames - 1)  # Clamp to valid range
                
                # Extract the full coordinate array for the specified frame
                # coords should have shape (num_atoms, 3)
                coords = coords_dataset[frame]
                
                # Check correct shape - should be 2D with second dimension of 3
                if coords.ndim != 2 or coords.shape[1] != 3:
                    logging.error(f"Unexpected coordinate shape: {coords.shape} for domain {self.domain_id}")
                    return None

                # Extract residue information
                resids = f[self.domain_id]['resid'][:].tolist()
                resnames = [name.decode('utf-8') if isinstance(name, bytes) else str(name) for name in f[self.domain_id]['resname'][:]]

                # Verify we have correct number of atoms
                if len(resids) != coords.shape[0]:
                    logging.warning(f"Mismatch between residue IDs ({len(resids)}) and coordinates ({coords.shape[0]})")
                    # Truncate to the smaller size if needed
                    min_size = min(len(resids), coords.shape[0])
                    resids = resids[:min_size]
                    resnames = resnames[:min_size]
                    coords = coords[:min_size]

                return coords, resids, resnames
        except Exception as e:
            logging.error(f"Failed to extract coordinate data: {e}")
            import traceback
            logging.error(traceback.format_exc())
            return None

def process_domains(domain_ids: List[str], data_dir: str, config: Dict[str, Any],
                    num_cores: int = 1) -> Dict[str, Any]:
    """
    Process multiple domains in parallel.

    Args:
        domain_ids: List of domain IDs to process
        data_dir: Directory containing H5 files
        config: Configuration dictionary
        num_cores: Number of CPU cores to use

    Returns:
        Dictionary with processing results
    """
    # Determine number of cores to use
    max_cores = os.cpu_count() - 2 if os.cpu_count() > 2 else 1
    n_cores = min(num_cores if num_cores > 0 else max_cores, max_cores)

    results = {}
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        future_to_domain = {}
        for domain_id in domain_ids:
            h5_path = os.path.join(data_dir, f"mdcath_dataset_{domain_id}.h5")
            if not os.path.exists(h5_path):
                logging.warning(f"H5 file not found for domain {domain_id}")
                continue

            future = executor.submit(_process_single_domain, h5_path, config)
            future_to_domain[future] = domain_id

        for future in as_completed(future_to_domain):
            domain_id = future_to_domain[future]
            try:
                result = future.result()
                results[domain_id] = result
            except Exception as e:
                logging.error(f"Error processing domain {domain_id}: {e}")
                results[domain_id] = {"success": False, "error": str(e)}

    return results

def _process_single_domain(h5_path: str, config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process a single domain (helper function for parallel processing).

    Args:
        h5_path: Path to H5 file
        config: Configuration dictionary

    Returns:
        Dictionary with processing results
    """
    loader = H5DataLoader(h5_path, config)
    domain_id = loader.domain_id

    results = {"domain_id": domain_id, "success": False}

    # Extract RMSF data for all temperatures and replicas
    temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
    num_replicas = config.get("num_replicas", 5)

    rmsf_data = {}
    for temp in temps:
        rmsf_data[temp] = {}
        for r in range(num_replicas):
            replica = str(r)
            df = loader.extract_rmsf(temp, replica)
            if df is not None:
                rmsf_data[temp][replica] = df

    results["rmsf_data"] = rmsf_data

    # Extract PDB data
    pdb_str = loader.extract_pdb()
    if pdb_str:
        results["pdb_data"] = pdb_str

    # Extract DSSP data
    dssp_data = {}
    for temp in temps:
        dssp_data[temp] = {}
        for r in range(num_replicas):
            replica = str(r)
            df = loader.extract_dssp(temp, replica)
            if df is not None:
                dssp_data[temp][replica] = df

    results["dssp_data"] = dssp_data
    results["success"] = True

    return results
