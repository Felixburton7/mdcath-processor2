#!/usr/bin/env python3
"""
Processing module for voxelizing protein structures using aposteriori.
"""

import os
import logging
import subprocess
import traceback
import sys
from typing import Dict, Any, Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

def voxelize_domains(pdb_results: Dict[str, Any], config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Voxelize multiple domains by processing the PDB directory directly.
    Handles the interactive confirmation prompt automatically.
    
    Args:
        pdb_results: Dictionary with PDB processing results
        config: Configuration dictionary
    
    Returns:
        Dictionary with voxelization results
    """
    from tqdm import tqdm
    import glob
    import shutil
    
    output_dir = config.get("output", {}).get("base_dir", "./outputs")
    voxel_dir = os.path.join(output_dir, "voxelized")
    pdb_dir = os.path.join(output_dir, "pdbs")
    
    os.makedirs(voxel_dir, exist_ok=True)
    
    # Check if aposteriori's make-frame-dataset is available
    aposteriori_path = shutil.which("make-frame-dataset")
    aposteriori_available = aposteriori_path is not None
    
    if not aposteriori_available:
        logging.error("make-frame-dataset command not found. Please install aposteriori with: pip install aposteriori")
        return {"success": False, "error": "aposteriori not available"}
    
    logging.info(f"Using aposteriori's make-frame-dataset from: {aposteriori_path}")
    
    # Count valid PDB files
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    
    if not pdb_files:
        logging.error(f"No PDB files found for voxelization in directory: {pdb_dir}")
        return {"success": False, "error": "No PDB files found"}
    
    logging.info(f"Found {len(pdb_files)} PDB files for voxelization in {pdb_dir}")
    
    # Get voxelization parameters from config
    voxel_config = config.get("processing", {}).get("voxelization", {})
    frame_edge_length = voxel_config.get("frame_edge_length", 12.0)
    voxels_per_side = voxel_config.get("voxels_per_side", 21)
    atom_encoder = voxel_config.get("atom_encoder", "CNOCBCA")
    encode_cb = voxel_config.get("encode_cb", True)
    compression_gzip = voxel_config.get("compression_gzip", True)
    voxelise_all_states = voxel_config.get("voxelise_all_states", False)
    
    # Create output name
    output_name = "mdcath_voxelized"
    
    # Build aposteriori command to process the PDB directory
    cmd = [
        aposteriori_path,
        "-o", voxel_dir,  # Output directory
        "-n", output_name,  # Output filename
        "-v",  # Enable verbose output
        "-e", ".pdb",  # File extension to look for
        "--frame-edge-length", str(frame_edge_length),
        "--voxels-per-side", str(voxels_per_side),
        "-ae", atom_encoder,
        "-cb", str(encode_cb).lower(),
        "-comp", str(compression_gzip).lower(),
        "-vas", str(voxelise_all_states).lower(),
        pdb_dir  # Directory containing PDB files
    ]
    
    # Log the exact command being run
    logging.info(f"Running aposteriori command: {' '.join(cmd)}")
    
    try:
        # Use Popen instead of run to handle interactive prompts
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # Send 'y' to confirm when prompted, with a newline character
        stdout, stderr = process.communicate(input='y\n')
        
        # Log the output
        logging.info(f"aposteriori stdout:\n{stdout}")
        if stderr:
            logging.info(f"aposteriori stderr:\n{stderr}")
        
        # Check return code
        if process.returncode != 0:
            logging.error(f"aposteriori command failed with return code {process.returncode}")
        
        # Check for expected output files (try both .hdf5 and .h5 extensions)
        output_file_hdf5 = os.path.join(voxel_dir, f"{output_name}.hdf5")
        output_file_h5 = os.path.join(voxel_dir, f"{output_name}.h5")
        
        if os.path.exists(output_file_hdf5):
            logging.info(f"Voxelization complete. Output: {output_file_hdf5}")
            file_size_mb = os.path.getsize(output_file_hdf5) / (1024 * 1024)
            logging.info(f"Output file size: {file_size_mb:.2f} MB")
            
            # Count number of residue frames
            try:
                import h5py
                with h5py.File(output_file_hdf5, 'r') as f:
                    num_frames = 0
                    for pdb_code in f.keys():
                        for chain_id in f[pdb_code].keys():
                            num_frames += len(f[pdb_code][chain_id].keys())
                logging.info(f"Dataset contains {num_frames} residue frames")
            except Exception as e:
                logging.warning(f"Could not count frames in output file: {e}")
            
            return {
                "success": True,
                "output_file": output_file_hdf5,
                "file_size_mb": file_size_mb,
                "num_frames": num_frames if 'num_frames' in locals() else None
            }
        elif os.path.exists(output_file_h5):
            logging.info(f"Voxelization complete. Output: {output_file_h5}")
            file_size_mb = os.path.getsize(output_file_h5) / (1024 * 1024)
            logging.info(f"Output file size: {file_size_mb:.2f} MB")
            
            # Count number of residue frames
            try:
                import h5py
                with h5py.File(output_file_h5, 'r') as f:
                    num_frames = 0
                    for pdb_code in f.keys():
                        for chain_id in f[pdb_code].keys():
                            num_frames += len(f[pdb_code][chain_id].keys())
                logging.info(f"Dataset contains {num_frames} residue frames")
            except Exception as e:
                logging.warning(f"Could not count frames in output file: {e}")
            
            return {
                "success": True,
                "output_file": output_file_h5,
                "file_size_mb": file_size_mb,
                "num_frames": num_frames if 'num_frames' in locals() else None
            }
        else:
            logging.error(f"Voxelization failed: No output file generated")
            return {
                "success": False,
                "error": "Output file not found",
                "command": ' '.join(cmd),
                "stdout": stdout,
                "stderr": stderr
            }
        
    except Exception as e:
        logging.error(f"Voxelization error: {e}")
        logging.error(traceback.format_exc())
        return {
            "success": False, 
            "error": str(e)
        }