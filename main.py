#!/usr/bin/env python3
"""
Main entry point for mdCATH dataset processing.
"""

import os
import sys
import logging
import argparse
import yaml
import time
from typing import Dict, Any, List, Optional

from src.mdcath.core.data_loader import process_domains
from src.mdcath.processing.rmsf import process_rmsf_data
from src.mdcath.processing.pdb import process_pdb_data
from src.mdcath.processing.features import process_ml_features
from src.mdcath.processing.visualization import generate_visualizations
from src.mdcath.processing.voxelizer import voxelize_domains


def setup_logging(config: Dict[str, Any]) -> None:
    """
    Set up logging based on configuration.
    
    Args:
        config: Configuration dictionary
    """
    log_level = config.get("logging", {}).get("level", "INFO")
    console_level = config.get("logging", {}).get("console_level", log_level)
    file_level = config.get("logging", {}).get("file_level", log_level)
    
    # Set up root logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)  # Allow all levels and filter in handlers
    
    # Remove existing handlers if any
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
    
    # Create handlers
    file_handler = logging.FileHandler("mdcath_processing.log")
    file_handler.setLevel(getattr(logging, file_level))
    
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(getattr(logging, console_level))
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load configuration from file or use default.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        Configuration dictionary
    """
    default_config_path = os.path.join("src", "mdcath", "config", "default_config.yaml")
    
    if config_path is None or not os.path.exists(config_path):
        config_path = default_config_path
    
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        return config
    except Exception as e:
        logging.error(f"Failed to load configuration: {e}")
        # Create minimal default config
        return {
            "input": {"mdcath_folder": "/mnt/datasets/MD_CATH/data", "domain_ids": ["12asA00"]},
            "temperatures": [320, 348, 379, 413, 450],
            "num_replicas": 5,
            "output": {"base_dir": "./outputs"},
            "logging": {"level": "INFO", "console_level": "INFO", "file_level": "DEBUG"}
        }


def create_output_directories(config: Dict[str, Any]) -> None:
    """
    Create output directories based on configuration.
    
    Args:
        config: Configuration dictionary
    """
    base_dir = config.get("output", {}).get("base_dir", "./outputs")
    
    directories = [
        os.path.join(base_dir, "pdbs"),
        os.path.join(base_dir, "frames"),
        os.path.join(base_dir, "RMSF", "replicas"),
        os.path.join(base_dir, "RMSF", "replica_average", "average"),
        os.path.join(base_dir, "ML_features"),
        os.path.join(base_dir, "visualizations"),
        os.path.join(base_dir, "voxelized")
    ]
    
    for directory in directories:
        os.makedirs(directory, exist_ok=True)


def process_mdcath(config: Dict[str, Any]) -> int:
    """
    Process mdCATH dataset based on configuration.
    
    Args:
        config: Configuration dictionary
        
    Returns:
        Exit code (0 for success, non-zero for error)
    """
    try:
        # Start time
        start_time = time.time()
        
        # Create output directories
        create_output_directories(config)
        
        # Get domain IDs and data directory
        domain_ids = config.get("input", {}).get("domain_ids", [])
        mdcath_folder = config.get("input", {}).get("mdcath_folder", "/mnt/datasets/MD_CATH/data")
        
        # Process domains
        logging.info(f"Processing {len(domain_ids)} domains from {mdcath_folder}")
        
        num_cores = config.get("performance", {}).get("num_cores", 0)
        domain_results = {}
        
        try:
            domain_results = process_domains(domain_ids, mdcath_folder, config, num_cores)
            domains_processed = sum(1 for r in domain_results.values() if r.get("success", False))
            logging.info(f"Successfully processed {domains_processed}/{len(domain_ids)} domains")
        except Exception as e:
            logging.error(f"Domain processing failed: {e}")
            import traceback
            traceback.print_exc()
            domain_results = {}
        
        # Process RMSF data
        rmsf_results = {"save_success": False}
        try:
            logging.info("Processing RMSF data")
            rmsf_results = process_rmsf_data(domain_results, config)
        except Exception as e:
            logging.error(f"RMSF processing failed: {e}")
            import traceback
            traceback.print_exc()
        
        # Process PDB data
        pdb_results = {}
        try:
            logging.info("Processing PDB data")
            pdb_results = process_pdb_data(domain_results, config)
        except Exception as e:
            logging.error(f"PDB processing failed: {e}")
            import traceback
            traceback.print_exc()
        
        # Generate ML features
        ml_results = {"success": False}
        try:
            logging.info("Generating ML features")
            ml_results = process_ml_features(rmsf_results, pdb_results, domain_results, config)
        except Exception as e:
            logging.error(f"ML feature generation failed: {e}")
            import traceback
            traceback.print_exc()
        
        # Voxelize domains
        voxel_results = {"success": False}
        try:
            logging.info("Voxelizing domains")
            voxel_results = voxelize_domains(pdb_results, config)
            # if config.get("processing", {}).get("voxelization", {}).get("process_frames", False):
            #     frame_voxel_results = voxelize_frames(config)
            #     processing_results["frame_voxelization"] = frame_voxel_results
        except Exception as e:
            logging.error(f"Voxelization failed: {e}")
            import traceback
            traceback.print_exc()
        
        # Generate visualizations
        vis_results = {}
        try:
            logging.info("Generating visualizations")
            vis_results = generate_visualizations(rmsf_results, ml_results, domain_results, config)
        except Exception as e:
            logging.error(f"Visualization generation failed: {e}")
            import traceback
            traceback.print_exc()
        
        # Summarize results
        logging.info("======== Processing Summary ========")
        domains_processed = sum(1 for r in domain_results.values() if r.get("success", False))
        logging.info(f"Domains processed: {domains_processed}/{len(domain_ids)}")
        logging.info(f"RMSF data saved: {rmsf_results.get('save_success', False)}")
        logging.info(f"PDB files saved: {sum(1 for r in pdb_results.values() if r.get('pdb_saved', False))}")
        logging.info(f"Feature files generated: {ml_results.get('success', False)}")
        
        # Handle voxel_results properly whether it's a dict or bool
        domains_voxelized = 0
        if isinstance(voxel_results, dict):
            if "success" in voxel_results and voxel_results["success"]:
                if "domain_ids" in voxel_results:
                    domains_voxelized = len(voxel_results.get("domain_ids", []))
                elif "num_domains" in voxel_results:
                    domains_voxelized = voxel_results.get("num_domains", 0)
                else:
                    domains_voxelized = 1  # At least something was successful
        elif isinstance(voxel_results, bool) and voxel_results:
            domains_voxelized = 1
        
        logging.info(f"Domains voxelized: {domains_voxelized}")
        
        # Calculate elapsed time
        elapsed_time = time.time() - start_time
        logging.info(f"Total processing time: {elapsed_time:.2f} seconds")
        
        return 0  # Success
    except Exception as e:
        logging.error(f"Processing failed with error: {str(e)}")
        import traceback
        traceback.print_exc()
        return 1  # Error

def voxelize_frames(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Voxelize frame directories for each temperature if specified in config.
    
    Args:
        config: Configuration dictionary
    
    Returns:
        Dictionary with voxelization results
    """
    import subprocess
    import shutil
    import os
    import logging
    from tqdm import tqdm
    
    output_dir = config.get("output", {}).get("base_dir", "./outputs")
    voxel_dir = os.path.join(output_dir, "voxelized")
    frames_dir = os.path.join(output_dir, "frames")
    
    os.makedirs(voxel_dir, exist_ok=True)
    
    # Check if aposteriori's make-frame-dataset is available
    aposteriori_path = shutil.which("make-frame-dataset")
    if not aposteriori_path:
        logging.warning("make-frame-dataset command not found, skipping frame voxelization")
        return {"success": False, "error": "aposteriori not available"}
    
    # Get voxelization parameters from config
    voxel_config = config.get("processing", {}).get("voxelization", {})
    process_frames = voxel_config.get("process_frames", False)
    
    if not process_frames:
        logging.info("Frame voxelization not requested in config, skipping")
        return {"success": True, "skipped": True, "message": "Frame voxelization not requested"}
    
    # Get temperatures and replicas to process
    temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
    num_replicas = config.get("num_replicas", 5)
    replicas = [str(r) for r in range(num_replicas)]
    
    results = {"success": True, "temp_results": {}}
    
    # Process each temperature and replica
    for temp in tqdm(temps, desc="Processing temperatures"):
        temp_dir = os.path.join(frames_dir, "replica_0", temp)  # Use replica_0 as reference
        if not os.path.exists(temp_dir) or not os.listdir(temp_dir):
            logging.warning(f"No frames found for temperature {temp}, skipping")
            continue
        
        # Create a temporary directory to hold symlinks to all PDB frame files for this temperature
        import tempfile
        with tempfile.TemporaryDirectory() as temp_frame_dir:
            # Create symlinks to all PDB frame files across all replicas for this temperature
            frame_count = 0
            for replica in replicas:
                replica_dir = os.path.join(frames_dir, f"replica_{replica}", temp)
                if not os.path.exists(replica_dir):
                    continue
                
                for frame_file in os.listdir(replica_dir):
                    if frame_file.endswith(".pdb"):
                        # Create a unique name for each frame file
                        frame_path = os.path.join(replica_dir, frame_file)
                        domain_id = frame_file.split("_")[0]
                        unique_name = f"{domain_id}_r{replica}_{frame_file}"
                        os.symlink(frame_path, os.path.join(temp_frame_dir, unique_name))
                        frame_count += 1
            
            if frame_count == 0:
                logging.warning(f"No valid frame files found for temperature {temp}")
                results["temp_results"][temp] = {"success": False, "error": "No valid frames"}
                continue
            
            # Get parameters from config
            frame_edge_length = voxel_config.get("frame_edge_length", 12.0)
            voxels_per_side = voxel_config.get("voxels_per_side", 21)
            atom_encoder = voxel_config.get("atom_encoder", "CNOCBCA")
            encode_cb = voxel_config.get("encode_cb", True)
            compression_gzip = voxel_config.get("compression_gzip", True)
            voxelise_all_states = voxel_config.get("voxelise_all_states", False)
            
            # Create output name for this temperature
            output_name = f"frames_temp{temp}_voxelized"
            
            # Build aposteriori command for this temperature's frames
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
                temp_frame_dir
            ]
            
            # Run the command with verbose output
            logging.info(f"Running aposteriori command for temperature {temp}: {' '.join(cmd)}")
            try:
                timeout = 300  # 5 minutes timeout
                result = subprocess.run(cmd, check=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                      text=True, timeout=timeout)
                
                logging.debug(f"Command stdout: {result.stdout}")
                if result.stderr:
                    logging.warning(f"Command stderr: {result.stderr}")
                
                # Check for expected output files
                output_file_hdf5 = os.path.join(voxel_dir, f"{output_name}.hdf5")
                output_file_h5 = os.path.join(voxel_dir, f"{output_name}.h5")
                
                if os.path.exists(output_file_hdf5):
                    logging.info(f"Voxelization complete for temperature {temp}. Output: {output_file_hdf5}")
                    results["temp_results"][temp] = {
                        "success": True,
                        "output_file": output_file_hdf5,
                        "file_size_mb": os.path.getsize(output_file_hdf5) / (1024 * 1024)
                    }
                elif os.path.exists(output_file_h5):
                    logging.info(f"Voxelization complete for temperature {temp}. Output: {output_file_h5}")
                    results["temp_results"][temp] = {
                        "success": True,
                        "output_file": output_file_h5,
                        "file_size_mb": os.path.getsize(output_file_h5) / (1024 * 1024)
                    }
                else:
                    logging.warning(f"Voxelization failed for temperature {temp}: No output file generated")
                    results["temp_results"][temp] = {
                        "success": False,
                        "error": "Output file not found"
                    }
                    results["success"] = False
                
            except subprocess.TimeoutExpired:
                logging.warning(f"Voxelization timeout for temperature {temp} after {timeout} seconds")
                results["temp_results"][temp] = {
                    "success": False, 
                    "error": f"Command timed out after {timeout} seconds"
                }
                results["success"] = False
            except Exception as e:
                logging.warning(f"Voxelization error for temperature {temp}: {e}")
                results["temp_results"][temp] = {
                    "success": False, 
                    "error": str(e)
                }
                results["success"] = False
    
    # Count successful voxelizations
    successful = sum(1 for t, r in results["temp_results"].items() if r.get("success", False))
    results["temps_processed"] = len(temps)
    results["temps_successful"] = successful
    
    if successful > 0:
        logging.info(f"Successfully voxelized frames for {successful}/{len(temps)} temperatures")
    else:
        logging.warning("Failed to voxelize frames for any temperature")
        
    return results


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Parsed arguments
    """
    parser = argparse.ArgumentParser(description="Process mdCATH dataset for ML feature extraction")
    parser.add_argument("--config", "-c", type=str, help="Path to configuration file")
    parser.add_argument("--domains", "-d", type=str, nargs="+", help="List of domain IDs to process")
    parser.add_argument("--output", "-o", type=str, help="Output directory")
    parser.add_argument("--verbose", "-v", action="store_true", help="Enable verbose logging")
    
    return parser.parse_args()


def main() -> int:
    """
    Main entry point.
    
    Returns:
        Exit code
    """
    # Parse arguments
    args = parse_arguments()
    
    # Load configuration
    config = load_config(args.config)
    
    # Override with command line arguments
    if args.domains:
        config["input"]["domain_ids"] = args.domains
    if args.output:
        config["output"]["base_dir"] = args.output
    if args.verbose:
        config["logging"]["console_level"] = "DEBUG"
    
    # Configure logging
    setup_logging(config)
    
    # Process mdCATH dataset
    return process_mdcath(config)


if __name__ == "__main__":
    sys.exit(main())