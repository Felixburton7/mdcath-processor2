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
        # Write original PDB to a temporary file
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
    Clean and fix a PDB file for downstream processing using pdbUtils or fallback.

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
    Clean PDB file using the pdbUtils library for better compatibility with your pipeline.

    Args:
        input_pdb_path: Path to input PDB file
        output_pdb_path: Path to save the cleaned PDB file
        config: Configuration dictionary

    Returns:
        Boolean indicating if cleaning was successful
    """
    try:
        # First, check if we need to stop after TER
        stop_after_ter = config.get("pdb_cleaning", {}).get("stop_after_ter", True)
        
        # If we need to stop after TER, we need to process the file line by line first
        if stop_after_ter:
            with open(input_pdb_path, 'r') as f:
                lines = f.readlines()
                
            # Find the first TER marker
            ter_index = -1
            for i, line in enumerate(lines):
                if line.startswith("TER"):
                    ter_index = i
                    break
            
            # Create a temporary file with content up to TER marker (or all if no TER found)
            temp_input_path = input_pdb_path + ".preproc"
            with open(temp_input_path, 'w') as f:
                if ter_index > 0:
                    f.writelines(lines[:ter_index+1])  # Include the TER line
                    if not any(line.startswith("END") for line in lines[:ter_index+1]):
                        f.write("END\n")  # Add END if not present
                else:
                    f.writelines(lines)  # Use all lines if no TER found
                    
            # Now use this as our input
            actual_input = temp_input_path
        else:
            actual_input = input_pdb_path
            
        # Convert PDB to DataFrame using pdbUtils
        pdb_df = pdbUtils.pdb2df(actual_input)
        initial_atoms = len(pdb_df)

        clean_config = config.get("pdb_cleaning", {})

        # Fix ONLY the first CAY to N - more targeted approach
        if clean_config.get("fix_first_cay", True):
            atom_name_col = "ATOM_NAME" if "ATOM_NAME" in pdb_df.columns else "atom_name"
            if atom_name_col in pdb_df.columns:
                # Find the first CAY atom
                cay_rows = pdb_df[pdb_df[atom_name_col] == "CAY"]
                if not cay_rows.empty:
                    first_cay_idx = cay_rows.index[0]
                    # Change only the first CAY to N
                    pdb_df.at[first_cay_idx, atom_name_col] = "N"
                    logging.info("Changed first CAY atom to N")

        # Replace HSD/HSE/HSP with HIS
        if clean_config.get("correct_unusual_residue_names", True):
            if "RES_NAME" in pdb_df.columns:
                pdb_df["RES_NAME"] = pdb_df["RES_NAME"].apply(
                    lambda x: "HIS" if str(x).strip() in ["HSD", "HSE", "HSP"] else x
                )
            elif "resName" in pdb_df.columns:
                pdb_df["resName"] = pdb_df["resName"].apply(
                    lambda x: "HIS" if str(x).strip() in ["HSD", "HSE", "HSP"] else x
                )

        # Replace chain 0 with A
        if clean_config.get("replace_chain_0_with_A", True):
            chain_col = "CHAIN_ID" if "CHAIN_ID" in pdb_df.columns else "chainID"
            if chain_col in pdb_df.columns:
                pdb_df[chain_col] = pdb_df[chain_col].apply(
                    lambda x: "A" if str(x).strip() == "0" else x
                )

        # Fix atom numbering
        if clean_config.get("fix_atom_numbering", True):
            atom_num_col = "ATOM_NUM" if "ATOM_NUM" in pdb_df.columns else "atomNum"
            if atom_num_col in pdb_df.columns:
                pdb_df[atom_num_col] = range(1, len(pdb_df) + 1)

        # Remove hydrogens if specified
        if clean_config.get("remove_hydrogens", False):
            elem_col = "ELEMENT" if "ELEMENT" in pdb_df.columns else "element"
            if elem_col in pdb_df.columns:
                pdb_df = pdb_df[pdb_df[elem_col] != "H"]

        # Remove water/ions if config says so (TIP, WAT, HOH, SOD, CLA, chain W)
        if clean_config.get("remove_solvent_ions", False):
            skip_resnames = {"TIP", "WAT", "HOH", "SOD", "CLA"}
            res_col = "RES_NAME" if "RES_NAME" in pdb_df.columns else "resName"
            chain_col = "CHAIN_ID" if "CHAIN_ID" in pdb_df.columns else "chainID"
            if res_col in pdb_df.columns:
                pdb_df = pdb_df[~pdb_df[res_col].isin(skip_resnames)]
            if chain_col in pdb_df.columns:
                pdb_df = pdb_df[pdb_df[chain_col] != "W"]

        # Save cleaned PDB
        pdbUtils.df2pdb(pdb_df, output_pdb_path)

        # Add properly formatted CRYST1 record for DSSP compatibility
        # This is critical as the error logs show DSSP is failing due to malformed CRYST1
        with open(output_pdb_path, "r") as f:
            content = f.readlines()
        
        # Properly formatted CRYST1 record with proper spacing and floating point values
        cryst1_line = "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\n"
        
        # Check if CRYST1 already exists and modify/add as needed
        has_cryst1 = False
        for i, line in enumerate(content):
            if line.startswith("CRYST1"):
                has_cryst1 = True
                content[i] = cryst1_line
                break
        
        if not has_cryst1:
            content.insert(0, cryst1_line)
        
        with open(output_pdb_path, "w") as f:
            f.writelines(content)

        # Cleanup temporary file if created
        if stop_after_ter and os.path.exists(temp_input_path):
            os.remove(temp_input_path)

        final_atoms = len(pdb_df)
        logging.info(
            f"Cleaned PDB {os.path.basename(input_pdb_path)}: {initial_atoms} atoms → {final_atoms} atoms"
        )
        return True

    except Exception as e:
        logging.error(f"Failed to clean PDB with pdbUtils: {e}")
        import traceback
        logging.error(traceback.format_exc())
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

        clean_config = config.get("pdb_cleaning", {})
        has_cryst1 = any(line.strip().startswith("CRYST1") for line in lines)

        cleaned_lines = []
        # Add CRYST1 if missing
        if not has_cryst1 and clean_config.get("add_cryst1_record", True):
            # Properly formatted CRYST1 record with correct spacing
            cleaned_lines.append("CRYST1    100.000   100.000   100.000  90.00  90.00  90.00 P 1           1\n")

        # We can skip TIP, WAT, HOH, SOD, CLA, chain 'W' if remove_solvent_ions is True
        skip_solvent = clean_config.get("remove_solvent_ions", False)
        skip_resnames = {"TIP", "WAT", "HOH", "SOD", "CLA"}
        stop_after_ter = clean_config.get("stop_after_ter", True)
        
        # Track whether we've encountered TER
        ter_encountered = False
        
        # Fix atom numbering
        atom_num = 1

        for line in lines:
            # Check if we should stop processing (after TER marker)
            if stop_after_ter and ter_encountered:
                # Only keep certain record types after TER (like END)
                if line.startswith("END"):
                    cleaned_lines.append(line)
                continue
                
            # Process TER line if encountered
            if line.startswith("TER"):
                cleaned_lines.append(f"TER   {atom_num:5d}\n")
                if stop_after_ter:
                    ter_encountered = True
                continue
                
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    record_type = line[0:6].strip()
                    # atom_num from the original is ignored, we use our counter
                    atom_name = line[12:16].strip()
                    alt_loc = line[16:17].strip()
                    res_name = line[17:20].strip()
                    chain_id = line[21:22].strip()
                    res_num = line[22:26].strip()
                    ins_code = line[26:27].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    occ = line[54:60].strip() or "0.00"
                    temp_factor = line[60:66].strip() or "0.00"
                    # Derive element from atom name if not present
                    element = line[76:78].strip() if len(line) >= 78 else atom_name[0:1]
                    
                    # Fix chain '0' → 'A'
                    if clean_config.get("replace_chain_0_with_A", True) and chain_id == "0":
                        chain_id = "A"

                    # Convert HSD/HSE/HSP → HIS
                    if clean_config.get("correct_unusual_residue_names", True):
                        if res_name in ["HSD", "HSE", "HSP"]:
                            res_name = "HIS"

                    # Remove hydrogens if desired
                    if clean_config.get("remove_hydrogens", False) and (
                        element == "H" or atom_name.startswith("H")
                    ):
                        continue

                    # Remove water/ions if configured
                    if skip_solvent:
                        # If residue name is in skip set OR chain is W, skip
                        if res_name in skip_resnames or chain_id == "W":
                            continue

                    # Format line
                    new_line = (
                        f"{record_type:<6s}{atom_num:5d} {atom_name:<4s}{alt_loc:1s}"
                        f"{res_name:3s} {chain_id:1s}{res_num:4s}{ins_code:1s}"
                        f"   {x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{float(occ):6.2f}{float(temp_factor):6.2f}"
                        f"           {element:>2s}  \n"
                    )
                    cleaned_lines.append(new_line)
                    atom_num += 1
                    
                except ValueError as ve:
                    logging.warning(f"Error parsing line: {line.strip()} -> {ve}")
                    # Skip problematic lines
                    continue
                except Exception as e:
                    logging.warning(f"Other error processing line: {line.strip()} -> {e}")
                    # Skip problematic lines
                    continue
            elif not line.startswith("ATOM") and not line.startswith("HETATM"):
                # Keep non-ATOM lines as is (e.g. REMARK, unless it's after TER)
                cleaned_lines.append(line)

        # Add END record if not present
        if not any(line.startswith("END") for line in cleaned_lines):
            cleaned_lines.append("END\n")
            
        # Write cleaned PDB
        with open(output_pdb_path, 'w') as f:
            f.writelines(cleaned_lines)

        logging.info(f"Cleaned PDB {os.path.basename(input_pdb_path)} → {os.path.basename(output_pdb_path)}, {atom_num-1} atoms")
        return True

    except Exception as e:
        logging.error(f"Failed to clean PDB {input_pdb_path} with fallback method: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return False
    
    
def extract_frames(coords: np.ndarray,
                   resids: List[int],
                   resnames: List[str],
                   domain_id: str,
                   output_dir: str,
                   temperature: str,
                   replica: str,
                   config: Dict[str, Any]) -> bool:
    """
    Extract frames from coordinate data and save as PDB files.
    Uses the cleaned PDB as a template and updates only the coordinates.

    Args:
        coords: Coordinate data
        resids: Residue IDs
        resnames: Residue names
        domain_id: Domain identifier
        output_dir: Directory to save frame PDBs
        temperature: Temperature
        replica: Replica index
        config: Configuration dictionary

    Returns:
        Boolean indicating if extraction was successful
    """
    frame_selection = config.get("processing", {}).get("frame_selection", {})
    method = frame_selection.get("method", "rmsd")
    num_frames = frame_selection.get("num_frames", 1)

    try:
        # Create output directory for frames
        frame_dir = os.path.join(output_dir, "frames", f"replica_{replica}", temperature)
        os.makedirs(frame_dir, exist_ok=True)
        
        # Get the source cleaned PDB path to use as a template
        pdb_dir = os.path.join(output_dir, "pdbs")
        pdb_path = os.path.join(pdb_dir, f"{domain_id}.pdb")
        
        # Check if the cleaned PDB exists
        if not os.path.exists(pdb_path):
            logging.error(f"Cleaned PDB file not found for domain {domain_id}: {pdb_path}")
            return False
            
        # Read the cleaned PDB template
        with open(pdb_path, 'r') as f:
            pdb_lines = f.readlines()
        
        # For simplicity, extract only the last frame if num_frames == 1
        if num_frames == 1:
            frame_idx = -1
            logging.debug(f"Coordinate shape for domain {domain_id}: {coords.shape}")

            # Handle shape differences
            if coords.ndim == 1:
                logging.warning(f"Coordinates shape {coords.shape} - expected 2D array.")
                if len(coords) == 3:  # Single XYZ
                    frame_coords = np.array([coords])
                else:
                    logging.error(f"Cannot reshape coords with length {len(coords)}")
                    return False
            elif coords.ndim == 2:
                # shape: (atoms, xyz)
                frame_coords = coords
            elif coords.ndim == 3:
                # shape: (frames, atoms, xyz)
                frame_coords = coords[frame_idx]
            else:
                logging.error(f"Unsupported coordinate dims: {coords.ndim}")
                return False

            # Verify final shape is (atoms, 3)
            if frame_coords.ndim != 2 or frame_coords.shape[1] != 3:
                logging.error(f"Coordinates have invalid shape {frame_coords.shape}")
                return False

            # Create a mapping from resid to coordinates index
            resid_to_coord = {}
            for i, resid in enumerate(resids):
                if i < len(frame_coords):
                    if resid not in resid_to_coord:
                        resid_to_coord[resid] = []
                    resid_to_coord[resid].append(i)
            
            # Track which atom indices in the PDB have been mapped to coordinates
            atom_coord_mapping = {}
            
            # First pass: Process lines and build atom_resid mapping
            atom_resid_mapping = {}
            for i, line in enumerate(pdb_lines):
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        atom_num = int(line[6:11].strip())
                        res_num = int(line[22:26].strip())
                        atom_resid_mapping[atom_num] = res_num
                    except ValueError:
                        continue
            
            # Second pass: Create new PDB with updated coordinates where possible
            new_pdb_lines = []
            for line in pdb_lines:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        # Extract parts of the ATOM line we want to keep
                        record_type = line[0:6].strip()
                        atom_num = int(line[6:11].strip())
                        atom_name = line[12:16].strip()
                        alt_loc = line[16:17].strip()
                        res_name = line[17:20].strip()
                        chain_id = line[21:22].strip()
                        res_num = int(line[22:26].strip())
                        ins_code = line[26:27].strip()
                        
                        # Extract values from original line for fields we're not updating
                        occ = line[54:60].strip() or "0.00"
                        temp_factor = line[60:66].strip() or "0.00"
                        element = line[76:78].strip() if len(line) >= 78 else ""
                        
                        # Look up coordinates based on resid
                        if res_num in resid_to_coord and len(resid_to_coord[res_num]) > 0:
                            coord_idx = resid_to_coord[res_num][0]
                            # Update coordinates
                            x, y, z = frame_coords[coord_idx]
                            # Format new ATOM line with updated coordinates
                            new_line = (
                                f"{record_type:<6s}{atom_num:5d} {atom_name:<4s}{alt_loc:1s}"
                                f"{res_name:3s} {chain_id:1s}{res_num:4d}{ins_code:1s}"
                                f"   {x:8.3f}{y:8.3f}{z:8.3f}"
                                f"{float(occ):6.2f}{float(temp_factor):6.2f}"
                                f"           {element:>2s}  \n"
                            )
                            atom_coord_mapping[atom_num] = coord_idx
                        else:
                            # Keep original coordinates if resid not in coordinate data
                            new_line = line
                        
                        new_pdb_lines.append(new_line)
                    except ValueError:
                        # If we can't parse the line, keep it as is
                        new_pdb_lines.append(line)
                else:
                    # Non-ATOM lines (HEADER, REMARK, etc.) stay the same
                    new_pdb_lines.append(line)
            
            # Check if we've used all coordinates
            coords_used = len(atom_coord_mapping)
            logging.info(f"Used {coords_used}/{len(frame_coords)} coordinates for {domain_id}")
            
            # Write out the single-frame PDB
            frame_path = os.path.join(frame_dir, f"{domain_id}_frame.pdb")
            with open(frame_path, 'w') as f:
                f.writelines(new_pdb_lines)

            logging.info(f"Extracted 1-frame PDB for domain {domain_id} at T={temperature}, rep={replica}")
            return True

        else:
            # Not implemented
            logging.warning(f"Multiple-frame extraction not implemented; requested {num_frames}.")
            return False

    except Exception as e:
        logging.error(f"Failed to extract frames for domain {domain_id}: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return False

def process_pdb_data(domain_results: Dict[str, Dict[str, Any]], config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process PDB data for all domains: save cleaned PDB, optionally extract frames, etc.

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

    # Check for pdbUtils
    if PDBUTILS_AVAILABLE:
        logging.info("Using pdbUtils for PDB cleaning (recommended).")
    else:
        logging.warning("pdbUtils not available, using fallback cleaning method.")

    # Save cleaned PDBs
    logging.info("Processing PDB data for domains...")
    for domain_id, result in tqdm(domain_results.items(), desc="Processing PDB files"):
        if not result.get("success", False):
            continue

        pdb_data = result.get("pdb_data", "")
        if pdb_data:
            pdb_path = os.path.join(pdb_dir, f"{domain_id}.pdb")
            success = save_pdb_file(pdb_data, pdb_path, config)
            if success:
                logging.info(f"Saved cleaned PDB for domain {domain_id}")
                results[domain_id] = {"pdb_saved": True, "pdb_path": pdb_path}
            else:
                logging.error(f"Failed to save cleaned PDB for domain {domain_id}")
                results[domain_id] = {"pdb_saved": False}

    # Extract frames from HDF5 trajectories
    temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
    num_replicas = config.get("num_replicas", 5)

    logging.info("Extracting frames from trajectories...")
    for domain_id, result in tqdm(results.items(), desc="Extracting frames"):
        if not result.get("pdb_saved", False):
            continue

        h5_path = os.path.join(input_dir, f"mdcath_dataset_{domain_id}.h5")
        if not os.path.exists(h5_path):
            logging.warning(f"H5 file not found for domain {domain_id}: {h5_path}")
            continue

        loader = H5DataLoader(h5_path, config)

        # For each temperature and replica
        for temp in temps:
            for r in range(num_replicas):
                replica = str(r)
                coords_result = loader.extract_coordinates(temp, replica)
                if coords_result is not None:
                    coords, resids, resnames = coords_result
                    extract_success = extract_frames(
                        coords, resids, resnames,
                        domain_id, output_dir, temp, replica, config
                    )
                    if extract_success:
                        if "frames" not in results[domain_id]:
                            results[domain_id]["frames"] = []
                        results[domain_id]["frames"].append((temp, replica))

    return results
