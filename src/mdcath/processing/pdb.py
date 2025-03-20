#!/usr/bin/env python3
"""
Processing module for PDB data extraction and cleaning.
"""

import os
import logging
import numpy as np
from typing import Dict, Any, Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

# Import pdbUtils for better PDB handling
try:
    from pdbUtils import pdbUtils
    PDBUTILS_AVAILABLE = True
except ImportError:
    logging.warning("pdbUtils library not found. Installing fallback method for PDB cleaning.")
    PDBUTILS_AVAILABLE = False

def save_pdb_file(pdb_string: str, output_path: str, config: Dict[str, Any]) -> bool:
    """
    Save a PDB string to a file with cleaning applied.

    Args:
        pdb_string: PDB data as a string
        output_path: Path to save the cleaned PDB
        config: Configuration dictionary

    Returns:
        Boolean indicating if saving was successful
    """
    try:
        # Write original PDB to temporary file
        temp_path = output_path + ".temp"
        with open(temp_path, 'w') as f:
            f.write(pdb_string)

        # Clean the PDB file
        success = fix_pdb(temp_path, output_path, config)

        # Remove temporary file
        if os.path.exists(temp_path):
            os.remove(temp_path)

        return success
    except Exception as e:
        logging.error(f"Failed to save PDB file: {e}")
        return False

def fix_pdb(input_pdb_path: str, output_pdb_path: str, config: Dict[str, Any]) -> bool:
    """
    Clean and fix a PDB file for downstream processing using pdbUtils.

    Args:
        input_pdb_path: Path to input PDB file
        output_pdb_path: Path to save the cleaned PDB file
        config: Configuration dictionary

    Returns:
        Boolean indicating if cleaning was successful
    """
    if not os.path.isfile(input_pdb_path):
        logging.error(f"PDB file not found: {input_pdb_path}")
        return False

    try:
        if PDBUTILS_AVAILABLE:
            return fix_pdb_with_pdbutils(input_pdb_path, output_pdb_path, config)
        else:
            return fix_pdb_fallback(input_pdb_path, output_pdb_path, config)
    except Exception as e:
        logging.error(f"Failed to clean PDB {input_pdb_path}: {e}")
        return False

def fix_pdb_with_pdbutils(input_pdb_path: str, output_pdb_path: str, config: Dict[str, Any]) -> bool:
    """
    Clean PDB file using the pdbUtils library for better compatibility with aposteriori.

    Args:
        input_pdb_path: Path to input PDB file
        output_pdb_path: Path to save the cleaned PDB file
        config: Configuration dictionary

    Returns:
        Boolean indicating if cleaning was successful
    """
    try:
        # Convert PDB to DataFrame for easier manipulation
        pdb_df = pdbUtils.pdb2df(input_pdb_path)
        initial_atoms = len(pdb_df)
        
        # Extract cleaning configuration
        clean_config = config.get("pdb_cleaning", {})
        
        # Replace HSD/HSE/HSP with HIS
        if clean_config.get("correct_unusual_residue_names", True):
            if "RES_NAME" in pdb_df.columns:
                pdb_df["RES_NAME"] = pdb_df["RES_NAME"].apply(
                    lambda x: 'HIS' if str(x).strip() in ['HSD', 'HSE', 'HSP'] else x
                )
            elif "resName" in pdb_df.columns:
                pdb_df["resName"] = pdb_df["resName"].apply(
                    lambda x: 'HIS' if str(x).strip() in ['HSD', 'HSE', 'HSP'] else x
                )
        
        # Replace chain 0 with A
        if clean_config.get("replace_chain_0_with_A", True):
            chain_col = "CHAIN_ID" if "CHAIN_ID" in pdb_df.columns else "chainID"
            if chain_col in pdb_df.columns:
                pdb_df[chain_col] = pdb_df[chain_col].apply(
                    lambda x: 'A' if str(x).strip() == '0' else x
                )
        
        # Fix atom numbering
        if clean_config.get("fix_atom_numbering", True):
            atom_num_col = "ATOM_NUM" if "ATOM_NUM" in pdb_df.columns else "atomNum"
            if atom_num_col in pdb_df.columns:
                pdb_df[atom_num_col] = range(1, len(pdb_df) + 1)
        
        # Remove hydrogens if specified
        if clean_config.get("remove_hydrogens", False):
            element_col = "ELEMENT" if "ELEMENT" in pdb_df.columns else "element"
            if element_col in pdb_df.columns:
                pdb_df = pdb_df[pdb_df[element_col] != "H"]
        
        # Check for CRYST1 record and add if missing
        has_cryst1 = False
        if "RECORD_NAME" in pdb_df.columns:
            has_cryst1 = any(pdb_df["RECORD_NAME"].str.contains("CRYST1"))
        
        # Save cleaned PDB
        pdbUtils.df2pdb(pdb_df, output_pdb_path)
        
        # Add CRYST1 if missing and specified in config
        if not has_cryst1 and clean_config.get("add_cryst1_record", True):
            with open(output_pdb_path, 'r') as f:
                content = f.readlines()
            
            with open(output_pdb_path, 'w') as f:
                f.write("CRYST1   100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n")
                f.writelines(content)
        
        final_atoms = len(pdb_df)
        logging.info(f"Cleaned PDB {os.path.basename(input_pdb_path)}: {initial_atoms} atoms → {final_atoms} atoms")
        
        return True
    except Exception as e:
        logging.error(f"Failed to clean PDB with pdbUtils: {e}")
        return False

def fix_pdb_fallback(input_pdb_path: str, output_pdb_path: str, config: Dict[str, Any]) -> bool:
    """
    Fallback method to clean a PDB file when pdbUtils is not available.

    Args:
        input_pdb_path: Path to input PDB file
        output_pdb_path: Path to save the cleaned PDB file
        config: Configuration dictionary

    Returns:
        Boolean indicating if cleaning was successful
    """
    try:
        with open(input_pdb_path, 'r') as f:
            lines = f.readlines()

        # Check for and add CRYST1 if needed
        has_cryst1 = False
        for line in lines:
            if line.strip() and line.startswith("CRYST1"):
                has_cryst1 = True
                break

        cleaned_lines = []
        if not has_cryst1 and config.get("pdb_cleaning", {}).get("add_cryst1_record", True):
            cleaned_lines.append("CRYST1   100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n")

        # Process atom records
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    # Extract fields
                    record_type = line[0:6].strip()
                    atom_num = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    alt_loc = line[16:17].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21:22].strip()
                    res_num = line[22:26].strip()
                    ins_code = line[26:27].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    occupancy = line[54:60].strip() if line[54:60].strip() else "0.00"
                    temp_factor = line[60:66].strip() if line[60:66].strip() else "0.00"
                    element = line[76:78].strip() if len(line) >= 78 else ""

                    # Apply cleaning rules
                    clean_config = config.get("pdb_cleaning", {})
                    if clean_config.get("replace_chain_0_with_A", True) and chain_id == "0":
                        chain_id = "A"

                    if clean_config.get("correct_unusual_residue_names", True):
                        if res_name in ["HSD", "HSE", "HSP"]:
                            res_name = "HIS"

                    # Skip hydrogens if configured
                    if clean_config.get("remove_hydrogens", False) and element == "H":
                        continue

                    # Format the line according to PDB standard
                    new_line = f"{record_type:<6s}{atom_num:5d} {atom_name:<4s}{alt_loc:1s}{res_name:3s} {chain_id:1s}{res_num:4s}{ins_code:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{float(occupancy):6.2f}{float(temp_factor):6.2f}           {element:>2s}  \n"
                    cleaned_lines.append(new_line)
                except Exception as e:
                    logging.warning(f"Error processing line: {line.strip()} - {e}")
                    # Keep the original line in case of parsing errors
                    cleaned_lines.append(line)
            else:
                # Keep non-ATOM lines as they are
                cleaned_lines.append(line)

        # Write cleaned PDB
        with open(output_pdb_path, 'w') as f:
            f.writelines(cleaned_lines)

        return True
    except Exception as e:
        logging.error(f"Failed to clean PDB {input_pdb_path} with fallback method: {e}")
        return False

def extract_frames(coords: np.ndarray, resids: List[int], resnames: List[str],
                  domain_id: str, output_dir: str, temperature: str, replica: str,
                  config: Dict[str, Any]) -> bool:
    """
    Extract frames from coordinate data and save as PDB files.

    Args:
        coords: Coordinate data
        resids: Residue IDs
        resnames: Residue names
        domain_id: Domain identifier
        output_dir: Directory to save frame PDBs
        temperature: Temperature
        replica: Replica
        config: Configuration dictionary

    Returns:
        Boolean indicating if extraction was successful
    """
    frame_selection = config.get("processing", {}).get("frame_selection", {})
    method = frame_selection.get("method", "rmsd")
    num_frames = frame_selection.get("num_frames", 1)

    try:
        # Create output directory
        frame_dir = os.path.join(output_dir, "frames", f"replica_{replica}", temperature)
        os.makedirs(frame_dir, exist_ok=True)

        # For simplicity, extract the last frame if num_frames is 1
        if num_frames == 1:
            frame_idx = -1
            
            # Debug log to check coords shape
            logging.debug(f"Coordinate shape for domain {domain_id}: {coords.shape}")
            
            # Handle different coordinate shapes
            if coords.ndim == 1:
                logging.warning(f"Coordinates have shape {coords.shape} - expected 2D array. Reshaping if possible.")
                if len(coords) == 3:  # Single XYZ coordinate
                    # This is a single coordinate that needs to be reshaped
                    frame_coords = np.array([coords])  # Make it a 2D array with one row
                else:
                    logging.error(f"Cannot reshape coordinates with length {len(coords)}")
                    return False
            elif coords.ndim == 2:
                frame_coords = coords  # Already in the correct shape
            elif coords.ndim == 3:
                # If this is a multi-frame array, select the specific frame
                frame_coords = coords[frame_idx]
            else:
                logging.error(f"Unsupported coordinate dimensions: {coords.ndim}")
                return False
            
            # Verify coordinates have correct shape now
            if frame_coords.ndim != 2 or frame_coords.shape[1] != 3:
                logging.error(f"Coordinates have invalid shape {frame_coords.shape} after processing")
                return False
            
            # Create a map from residue IDs to their indices
            residue_map = {}
            for i, resid in enumerate(resids):
                if i < len(frame_coords):  # Ensure we have coordinates for this index
                    if resid not in residue_map:
                        residue_map[resid] = {
                            'indices': [],
                            'resname': resnames[i]
                        }
                    residue_map[resid]['indices'].append(i)
            
            # Create PDB string for the frame
            pdb_lines = []
            # Add CRYST1 record
            pdb_lines.append("CRYST1   100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n")
            atom_num = 1
            
            # Select representative atoms for each residue
            for resid, info in residue_map.items():
                # Use the first atom of each residue
                if not info['indices']:
                    continue
                    
                atom_idx = info['indices'][0]
                
                if atom_idx < len(frame_coords):
                    try:
                        x, y, z = frame_coords[atom_idx]
                        resname = info['resname']
                        atom_name = "CA"  # Use CA as representative
                        chain_id = "A"
                        pdb_lines.append(f"ATOM  {atom_num:5d} {atom_name:<4s} {resname:3s} {chain_id:1s}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n")
                        atom_num += 1
                    except (IndexError, TypeError, ValueError) as e:
                        logging.warning(f"Error processing coordinates for residue {resid}: {e}")
                        continue
            
            # Check if we have any atoms processed
            if atom_num <= 1:
                logging.warning(f"No valid atoms found for {domain_id}, skipping frame extraction")
                return False

            # Write frame PDB
            frame_path = os.path.join(frame_dir, f"{domain_id}_frame.pdb")
            with open(frame_path, 'w') as f:
                f.writelines(pdb_lines)

            logging.info(f"Extracted frame for domain {domain_id}, temperature {temperature}, replica {replica}")
            return True
        else:
            logging.warning(f"Multiple frame extraction not implemented yet")
            return False
    except Exception as e:
        logging.error(f"Failed to extract frames: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return False

def process_pdb_data(domain_results: Dict[str, Dict[str, Any]], config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process PDB data for all domains.
    
    Args:
        domain_results: Dictionary with processing results for all domains
        config: Configuration dictionary
    
    Returns:
        Dictionary with PDB processing results
    """
    from src.mdcath.core.data_loader import H5DataLoader
    from tqdm import tqdm
    
    output_dir = config.get("output", {}).get("base_dir", "./outputs")
    input_dir = config.get("input", {}).get("mdcath_folder", "/mnt/datasets/MD_CATH/data")
    
    # Create output directories
    pdb_dir = os.path.join(output_dir, "pdbs")
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(pdb_dir, exist_ok=True)
    os.makedirs(frames_dir, exist_ok=True)
    
    results = {}
    
    # Check if pdbUtils is available
    if PDBUTILS_AVAILABLE:
        logging.info("Using pdbUtils for PDB cleaning (recommended for aposteriori compatibility)")
    else:
        logging.warning("pdbUtils not available, using fallback PDB cleaning method")
    
    # Process domains
    logging.info("Processing PDB data for domains")
    for domain_id, result in tqdm(domain_results.items(), desc="Processing PDB files"):
        if not result.get("success", False):
            continue
        
        # Save cleaned PDB
        pdb_data = result.get("pdb_data")
        if pdb_data:
            pdb_path = os.path.join(pdb_dir, f"{domain_id}.pdb")
            success = save_pdb_file(pdb_data, pdb_path, config)
            if success:
                logging.info(f"Saved cleaned PDB for domain {domain_id}")
                results[domain_id] = {"pdb_saved": True, "pdb_path": pdb_path}
            else:
                logging.error(f"Failed to save cleaned PDB for domain {domain_id}")
                results[domain_id] = {"pdb_saved": False}
    
    # Extract frames from trajectories
    temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
    num_replicas = config.get("num_replicas", 5)
    
    logging.info("Extracting frames from trajectories")
    for domain_id, result in tqdm(results.items(), desc="Extracting frames"):
        if not result.get("pdb_saved", False):
            continue
        
        # Create data loader for this domain
        h5_path = os.path.join(input_dir, f"mdcath_dataset_{domain_id}.h5")
        if not os.path.exists(h5_path):
            logging.warning(f"H5 file not found for domain {domain_id}: {h5_path}")
            continue
        
        loader = H5DataLoader(h5_path, config)
        
        # Extract frames for each temperature and replica
        for temp in temps:
            for r in range(num_replicas):
                replica = str(r)
                coords_result = loader.extract_coordinates(temp, replica)
                if coords_result is not None:
                    coords, resids, resnames = coords_result
                    extract_success = extract_frames(coords, resids, resnames, domain_id, 
                                                   output_dir, temp, replica, config)
                    if extract_success:
                        if "frames" not in results[domain_id]:
                            results[domain_id]["frames"] = []
                        results[domain_id]["frames"].append((temp, replica))
    
    return results