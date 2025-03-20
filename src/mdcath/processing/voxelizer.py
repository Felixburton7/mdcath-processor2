#!/usr/bin/env python3
"""
Processing module for voxelizing protein structures using aposteriori.
"""

import os
import logging
import subprocess
from typing import Dict, Any, Optional, List, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

def voxelize_domain(pdb_file: str, output_dir: str, config: Dict[str, Any]) -> Optional[str]:
    """
    Voxelize a cleaned PDB file using aposteriori's make-frame-dataset command.
    
    Args:
        pdb_file: Path to the cleaned PDB file
        output_dir: Directory to save voxelized output
        config: Configuration dictionary with voxelization parameters
    
    Returns:
        Path to the output file if successful, None otherwise
    """
    try:
        # Check if aposteriori is installed
        aposteriori_path = subprocess.run(["which", "make-frame-dataset"], 
                          check=False,
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE,
                          text=True).stdout.strip()
                          
        aposteriori_installed = aposteriori_path != ""
        
        if not aposteriori_installed:
            aposteriori_installed = False
            logging.warning("make-frame-dataset command not found, but proceeding as requested.")
        
        if not aposteriori_installed:
            # Skip voxelization but return success to keep pipeline running
            return {"success": True, "skipped": True, "message": "aposteriori check skipped"}
        
        # Get voxelization parameters from config
        voxel_config = config.get("processing", {}).get("voxelization", {})
        frame_edge_length = voxel_config.get("frame_edge_length", 12.0)
        voxels_per_side = voxel_config.get("voxels_per_side", 21)
        atom_encoder = voxel_config.get("atom_encoder", "CNOCBCA")
        encode_cb = voxel_config.get("encode_cb", True)
        compression_gzip = voxel_config.get("compression_gzip", True)
        voxelise_all_states = voxel_config.get("voxelise_all_states", False)
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Create output name
        domain_name = os.path.basename(pdb_file).split('.')[0]
        output_name = f"{domain_name}_voxelized"
        
        # Build aposteriori command
        cmd = [
            aposteriori_path,
            "-o", output_dir,
            "-n", output_name,
            "-v",
            "--frame-edge-length", str(frame_edge_length),
            "--voxels-per-side", str(voxels_per_side),
            "-ae", atom_encoder,
            "-cb", str(encode_cb).lower(),
            "-comp", str(compression_gzip).lower(),
            "-vas", str(voxelise_all_states).lower(),
            pdb_file
        ]
        
        # Run the command
        logging.info(f"Running aposteriori voxelization for {pdb_file}")
        result = subprocess.run(cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            logging.warning(f"Voxelization command returned non-zero status {result.returncode}")
            logging.debug(f"Command output: {result.stdout}")
            logging.debug(f"Command error: {result.stderr}")
        
        # Check for expected output files (try both .hdf5 and .h5 extensions)
        output_file_hdf5 = os.path.join(output_dir, f"{output_name}.hdf5")
        output_file_h5 = os.path.join(output_dir, f"{output_name}.h5")
        
        if os.path.exists(output_file_hdf5):
            logging.info(f"Voxelization completed: {output_file_hdf5}")
            return output_file_hdf5
        elif os.path.exists(output_file_h5):
            logging.info(f"Voxelization completed: {output_file_h5}")
            return output_file_h5
        else:
            logging.error(f"Voxelization failed: output file not found")
            logging.debug(f"Command output: {result.stdout}")
            logging.debug(f"Command error: {result.stderr}")
            return None
    except subprocess.CalledProcessError as e:
        logging.error(f"Voxelization failed: {e.stderr}")
        return None
    except Exception as e:
        logging.error(f"Voxelization failed with unexpected error: {e}")
        return None


# In src/mdcath/processing/voxelizer.py, we need to modify the voxelize_domains function

def voxelize_domains(pdb_results: Dict[str, Any], config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Voxelize multiple domains by processing all PDBs into a single output file.
    
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
        try:
            import aposteriori
            aposteriori_available = True
        except ImportError:
            aposteriori_available = False
            
    if not aposteriori_available:
        logging.warning("make-frame-dataset command not found, skipping voxelization")
        return {"success": True, "skipped": True, "message": "aposteriori not available"}
    else:
        logging.info(f"Using aposteriori's make-frame-dataset from: {aposteriori_path}")
    
    # Count valid PDB files
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    
    if not pdb_files:
        logging.warning("No valid PDB files found for voxelization")
        return {"success": False, "error": "No valid PDB files"}
    
    logging.info(f"Checking {len(pdb_files)} PDB files for voxelization")
    valid_pdbs = []
    for pdb_file in tqdm(pdb_files, desc="Validating PDBs"):
        # Validate PDB file
        try:
            with open(pdb_file, 'r') as f:
                content = f.read(1000)  # Read first 1000 chars to check for ATOM records
                if "ATOM" in content and os.path.getsize(pdb_file) > 0:
                    valid_pdbs.append(pdb_file)
                else:
                    logging.warning(f"PDB file {pdb_file} does not contain valid ATOM records, skipping")
        except Exception as e:
            logging.warning(f"Error validating PDB file {pdb_file}: {e}")
    
    if not valid_pdbs:
        logging.warning("No valid PDB files found for voxelization")
        return {"success": False, "error": "No valid PDB files"}
    
    # Create temporary directory with symlinks to all PDB files
    import tempfile
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create symlinks to all PDB files in temp_dir
        for pdb_file in valid_pdbs:
            basename = os.path.basename(pdb_file)
            os.symlink(pdb_file, os.path.join(temp_dir, basename))
        
        # Get voxelization parameters from config
        voxel_config = config.get("processing", {}).get("voxelization", {})
        frame_edge_length = voxel_config.get("frame_edge_length", 12.0)
        voxels_per_side = voxel_config.get("voxels_per_side", 21)
        atom_encoder = voxel_config.get("atom_encoder", "CNOCBCA")
        encode_cb = voxel_config.get("encode_cb", True)
        compression_gzip = voxel_config.get("compression_gzip", True)
        voxelise_all_states = voxel_config.get("voxelise_all_states", False)
        
        # Create combined dataset name
        output_name = "mdcath_voxelized"
        
        # Build aposteriori command for all domains
        cmd = [
            aposteriori_path,
            "-o", voxel_dir,
            "-n", output_name,
            "-v",
            "--frame-edge-length", str(frame_edge_length),
            "--voxels-per-side", str(voxels_per_side),
            "-ae", atom_encoder,
            "-cb", str(encode_cb).lower(),
            "-comp", str(compression_gzip).lower(),
            "-vas", str(voxelise_all_states).lower(),
            temp_dir      # Use the directory containing all PDB files
        ]
        
        # Run the command with verbose output to help with debugging
        logging.info(f"Running aposteriori command: {' '.join(cmd)}")
        try:
            # Set timeout and run the command
            timeout = 300  # 5 minutes timeout
            result = subprocess.run(cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                  text=True, timeout=timeout)
            
            # Print stdout and stderr for debugging
            logging.debug(f"Command stdout: {result.stdout}")
            if result.stderr:
                logging.warning(f"Command stderr: {result.stderr}")
            
            # Check for expected output files (try both .hdf5 and .h5 extensions)
            output_file_hdf5 = os.path.join(voxel_dir, f"{output_name}.hdf5")
            output_file_h5 = os.path.join(voxel_dir, f"{output_name}.h5")
            
            if os.path.exists(output_file_hdf5):
                logging.info(f"Voxelization complete. Output: {output_file_hdf5}")
                return {
                    "success": True,
                    "output_file": output_file_hdf5,
                    "file_size_mb": os.path.getsize(output_file_hdf5) / (1024 * 1024)
                }
            elif os.path.exists(output_file_h5):
                logging.info(f"Voxelization complete. Output: {output_file_h5}")
                return {
                    "success": True,
                    "output_file": output_file_h5,
                    "file_size_mb": os.path.getsize(output_file_h5) / (1024 * 1024)
                }
            else:
                logging.warning(f"Voxelization failed: No output file generated")
                return {
                    "success": False,
                    "error": "Output file not found",
                    "command": ' '.join(cmd),
                    "stdout": result.stdout,
                    "stderr": result.stderr
                }
            
        except subprocess.TimeoutExpired:
            logging.warning(f"Voxelization timeout after {timeout} seconds")
            return {
                "success": False, 
                "error": f"Command timed out after {timeout} seconds"
            }
        except Exception as e:
            logging.warning(f"Voxelization error: {e}")
            import traceback
            logging.error(traceback.format_exc())
            return {
                "success": False, 
                "error": str(e)
            }