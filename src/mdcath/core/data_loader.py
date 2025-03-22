

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
        RMSF is per-residue, so we build a unique residue-level list
        from the full 'resid'/'resname' arrays (which may be per-atom).
        
        Args:
            temperature: Temperature (e.g., "320")
            replica: Replica (e.g., "0")
        
        Returns:
            DataFrame with columns: [domain_id, resid, resname, rmsf_{temperature}]
            or None if extraction fails
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                # Check if temperature and replica exist
                if (temperature not in f[self.domain_id]) or (replica not in f[self.domain_id][temperature]):
                    logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
                    return None
                
                if 'rmsf' not in f[self.domain_id][temperature][replica]:
                    logging.warning(f"RMSF data not found for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None
                
                # RMSF data is typically length = number_of_residues
                rmsf_data = f[self.domain_id][temperature][replica]['rmsf'][:]

                # Extract the full, per-atom arrays
                resids_all = f[self.domain_id]['resid'][:]
                resnames_all = f[self.domain_id]['resname'][:]

                # Convert bytes -> string if needed
                resnames_all = [
                    rn.decode("utf-8") if isinstance(rn, bytes) else str(rn)
                    for rn in resnames_all
                ]

                # Build unique residue-level list
                # Map resid -> resname (the first occurrence of that resid)
                # This ensures one row per residue
                residue_dict = {}
                for i, resid_val in enumerate(resids_all):
                    if resid_val not in residue_dict:
                        residue_dict[resid_val] = resnames_all[i]

                unique_resids = sorted(residue_dict.keys())
                unique_resnames = [residue_dict[rid] for rid in unique_resids]

                # Check dimension mismatch
                if len(unique_resids) != len(rmsf_data):
                    logging.info(
                        f"Dimension mismatch: unique_resids {len(unique_resids)}, "
                        f"rmsf_data {len(rmsf_data)}"
                    )
                    # Attempt to align by length
                    if len(unique_resids) > len(rmsf_data):
                        logging.warning(
                            f"More unique residues ({len(unique_resids)}) than RMSF points ({len(rmsf_data)}) -- truncating residues"
                        )
                        unique_resids = unique_resids[:len(rmsf_data)]
                        unique_resnames = unique_resnames[:len(rmsf_data)]
                    else:
                        logging.warning(
                            f"Fewer unique residues ({len(unique_resids)}) than RMSF points ({len(rmsf_data)}) -- truncating RMSF"
                        )
                        rmsf_data = rmsf_data[:len(unique_resids)]
                    logging.info("Using unique residue-level alignment for RMSF data")

                # Create DataFrame with final 1:1 alignment
                df = pd.DataFrame({
                    'domain_id': self.domain_id,
                    'resid': unique_resids,
                    'resname': unique_resnames,
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
        DSSP is per-residue, so we build a unique residue-level list
        from the full 'resid'/'resname' arrays. Then align to DSSP codes.
        
        Args:
            temperature: Temperature (e.g., "320")
            replica: Replica (e.g., "0")
            frame: Frame index (default: -1 for last frame)

        Returns:
            DataFrame [domain_id, resid, resname, dssp] or None if extraction fails
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                if (temperature not in f[self.domain_id]) or (replica not in f[self.domain_id][temperature]):
                    logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
                    return None

                if 'dssp' not in f[self.domain_id][temperature][replica]:
                    logging.warning(f"DSSP data not found for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None
                    
                dssp_dataset = f[self.domain_id][temperature][replica]['dssp']

                # Number of frames
                num_frames = dssp_dataset.shape[0] if len(dssp_dataset.shape) > 0 else 0
                if num_frames == 0:
                    logging.warning(f"Empty DSSP dataset for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None

                # Convert negative frame index
                if frame < 0:
                    frame = num_frames + frame
                if frame < 0 or frame >= num_frames:
                    logging.warning(f"Frame index {frame} out of bounds (0-{num_frames-1}) for {self.domain_id}")
                    frame = max(0, min(frame, num_frames-1))  # clamp

                dssp_data = dssp_dataset[frame]

                # Full, per-atom arrays
                resids_all = f[self.domain_id]['resid'][:]
                resnames_all = f[self.domain_id]['resname'][:]
                resnames_all = [
                    rn.decode("utf-8") if isinstance(rn, bytes) else str(rn)
                    for rn in resnames_all
                ]

                # Build unique residue-level list
                residue_dict = {}
                for i, resid_val in enumerate(resids_all):
                    if resid_val not in residue_dict:
                        residue_dict[resid_val] = resnames_all[i]

                unique_resids = sorted(residue_dict.keys())
                unique_resnames = [residue_dict[rid] for rid in unique_resids]

                # DSSP codes might already be length = # of residues
                dssp_codes = [
                    c.decode("utf-8") if isinstance(c, bytes) else str(c)
                    for c in dssp_data
                ]

                if len(unique_resids) != len(dssp_codes):
                    logging.info(
                        f"Dimension mismatch in DSSP: unique_resids {len(unique_resids)}, dssp_codes {len(dssp_codes)}"
                    )
                    if len(unique_resids) > len(dssp_codes):
                        logging.warning(
                            f"More unique residues ({len(unique_resids)}) than DSSP codes ({len(dssp_codes)}) -- truncating residues"
                        )
                        unique_resids = unique_resids[:len(dssp_codes)]
                        unique_resnames = unique_resnames[:len(dssp_codes)]
                    else:
                        logging.warning(
                            f"Fewer unique residues ({len(unique_resids)}) than DSSP codes ({len(dssp_codes)}) -- truncating DSSP codes"
                        )
                        dssp_codes = dssp_codes[:len(unique_resids)]
                    logging.info("Using unique residue-level alignment for DSSP data")

                # Create final DataFrame
                df = pd.DataFrame({
                    'domain_id': self.domain_id,
                    'resid': unique_resids,
                    'resname': unique_resnames,
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
            Tuple of (coords, resids, resnames) where coords shape is (n_atoms, 3),
            or None if extraction fails.
        """
        try:
            with h5py.File(self.h5_path, 'r') as f:
                if (temperature not in f[self.domain_id]) or (replica not in f[self.domain_id][temperature]):
                    logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
                    return None

                if 'coords' not in f[self.domain_id][temperature][replica]:
                    logging.warning(f"Coordinate data not found for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None

                coords_dataset = f[self.domain_id][temperature][replica]['coords']
                num_frames = coords_dataset.shape[0] if coords_dataset.ndim > 0 else 0
                if num_frames == 0:
                    logging.warning(f"Empty coords dataset for domain {self.domain_id}, temperature {temperature}, replica {replica}")
                    return None

                # Convert negative frame index
                if frame < 0:
                    frame = num_frames + frame
                if frame < 0 or frame >= num_frames:
                    logging.warning(f"Frame index {frame} out of bounds (0-{num_frames-1}) for domain {self.domain_id}")
                    frame = max(0, min(frame, num_frames - 1))

                coords = coords_dataset[frame]  # shape (n_atoms, 3)
                if coords.ndim != 2 or coords.shape[1] != 3:
                    logging.error(f"Unexpected coordinate shape: {coords.shape} for domain {self.domain_id}")
                    return None

                resids_all = f[self.domain_id]['resid'][:].tolist()
                resnames_all = f[self.domain_id]['resname'][:]
                resnames_all = [
                    rn.decode("utf-8") if isinstance(rn, bytes) else str(rn)
                    for rn in resnames_all
                ]

                # Check shape alignment
                if len(resids_all) != coords.shape[0]:
                    logging.warning(
                        f"Mismatch between residue IDs ({len(resids_all)}) and coords ({coords.shape[0]})"
                    )
                    min_size = min(len(resids_all), coords.shape[0])
                    resids_all = resids_all[:min_size]
                    resnames_all = resnames_all[:min_size]
                    coords = coords[:min_size]

                return coords, resids_all, resnames_all

        except Exception as e:
            logging.error(f"Failed to extract coordinate data: {e}")
            import traceback
            logging.error(traceback.format_exc())
            return None


def process_domains(domain_ids: List[str], data_dir: str, config: Dict[str, Any],
                    num_cores: int = 1) -> Dict[str, Any]:
    """
    Process multiple domains in parallel.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

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
    """
    loader = H5DataLoader(h5_path, config)
    domain_id = loader.domain_id

    results = {"domain_id": domain_id, "success": False}

    # Extract RMSF data
    temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
    num_replicas = config.get("num_replicas", 5)

    rmsf_data = {}
    for temp in temps:
        rmsf_data[temp] = {}
        for r in range(num_replicas):
            replica = str(r)
            df_rmsf = loader.extract_rmsf(temp, replica)
            if df_rmsf is not None:
                rmsf_data[temp][replica] = df_rmsf
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
            df_dssp = loader.extract_dssp(temp, replica)
            if df_dssp is not None:
                dssp_data[temp][replica] = df_dssp
    results["dssp_data"] = dssp_data

    results["success"] = True
    return results
