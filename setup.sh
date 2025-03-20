# #!/bin/bash

# # Make the script exit on error
# set -e

# echo "Generating code files for mdCATH project..."

# # ------------------------------------
# # INIT FILES
# # ------------------------------------

# # Create src/mdcath/__init__.py
# cat > src/mdcath/__init__.py << 'EOFMARKER'
# """
# mdCATH - A package for processing mdCATH dataset for ML applications
# """

# __version__ = '0.1.0'
# EOFMARKER
# echo "Created src/mdcath/__init__.py"

# # Create src/mdcath/core/__init__.py
# cat > src/mdcath/core/__init__.py << 'EOFMARKER'
# """
# Core data loading and processing functions
# """
# EOFMARKER
# echo "Created src/mdcath/core/__init__.py"

# # Create src/mdcath/processing/__init__.py
# cat > src/mdcath/processing/__init__.py << 'EOFMARKER'
# """
# Data processing modules for mdCATH
# """
# EOFMARKER
# echo "Created src/mdcath/processing/__init__.py"

# # Create src/mdcath/config/__init__.py
# cat > src/mdcath/config/__init__.py << 'EOFMARKER'
# """
# Configuration handling for mdCATH
# """
# EOFMARKER
# echo "Created src/mdcath/config/__init__.py"

# # ------------------------------------
# # CORE MODULE FILES
# # ------------------------------------

# # Create src/mdcath/core/data_loader.py
# cat > src/mdcath/core/data_loader.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Core functionality for loading and processing H5 data from mdCATH dataset.
# """

# import os
# import h5py
# import logging
# import numpy as np
# import pandas as pd
# from typing import List, Dict, Tuple, Optional, Union, Any
# from concurrent.futures import ProcessPoolExecutor, as_completed

# class H5DataLoader:
#     """
#     Class for efficiently loading and extracting data from mdCATH H5 files.
#     Uses chunking/streaming to handle large files.
#     """

#     def __init__(self, h5_path: str, config: Dict[str, Any]):
#         """
#         Initialize the H5 data loader.

#         Args:
#             h5_path: Path to H5 file
#             config: Configuration dictionary
#         """
#         self.h5_path = h5_path
#         self.config = config
#         self.domain_id = os.path.basename(h5_path).replace("mdcath_dataset_", "").replace(".h5", "")
#         self._validate_h5()

#     def _validate_h5(self) -> bool:
#         """
#         Validate that the H5 file has the expected structure.

#         Returns:
#             Boolean indicating if the file is valid
#         """
#         try:
#             with h5py.File(self.h5_path, 'r') as f:
#                 # Check if domain exists
#                 if self.domain_id not in f:
#                     logging.error(f"Domain {self.domain_id} not found in {self.h5_path}")
#                     return False

#                 # Check for required temperature groups
#                 temps = [str(t) for t in self.config.get("temperatures", [320, 348, 379, 413, 450])]
#                 for temp in temps:
#                     if temp not in f[self.domain_id]:
#                         logging.warning(f"Temperature {temp} not found for domain {self.domain_id}")

#                 # Check for replica groups
#                 num_replicas = self.config.get("num_replicas", 5)
#                 for temp in temps:
#                     if temp in f[self.domain_id]:
#                         for r in range(num_replicas):
#                             if str(r) not in f[self.domain_id][temp]:
#                                 logging.warning(f"Replica {r} not found for temperature {temp} in domain {self.domain_id}")
#             return True
#         except Exception as e:
#             logging.error(f"Failed to validate H5 file {self.h5_path}: {e}")
#             return False

#     def extract_rmsf(self, temperature: str, replica: str) -> Optional[pd.DataFrame]:
#         """
#         Extract RMSF data for a specific temperature and replica.

#         Args:
#             temperature: Temperature (e.g., "320")
#             replica: Replica (e.g., "0")

#         Returns:
#             DataFrame with RMSF data or None if extraction fails
#         """
#         try:
#             with h5py.File(self.h5_path, 'r') as f:
#                 # Check if temperature and replica exist
#                 if temperature not in f[self.domain_id] or replica not in f[self.domain_id][temperature]:
#                     logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
#                     return None

#                 # Extract RMSF data
#                 rmsf_data = f[self.domain_id][temperature][replica]['rmsf'][:]

#                 # Extract residue information
#                 resids = f[self.domain_id]['resid'][:]
#                 resnames = [name.decode('utf-8') for name in f[self.domain_id]['resname'][:]]

#                 # Create DataFrame
#                 df = pd.DataFrame({
#                     'domain_id': self.domain_id,
#                     'resid': resids,
#                     'resname': resnames,
#                     f'rmsf_{temperature}': rmsf_data
#                 })

#                 return df
#         except Exception as e:
#             logging.error(f"Failed to extract RMSF data: {e}")
#             return None

#     def extract_pdb(self) -> Optional[str]:
#         """
#         Extract PDB data from the H5 file.

#         Returns:
#             PDB string or None if extraction fails
#         """
#         try:
#             with h5py.File(self.h5_path, 'r') as f:
#                 pdb_data = f[self.domain_id]['pdb'][()]
#                 if isinstance(pdb_data, bytes):
#                     return pdb_data.decode('utf-8')
#                 return str(pdb_data)
#         except Exception as e:
#             logging.error(f"Failed to extract PDB data: {e}")
#             return None

#     def extract_dssp(self, temperature: str, replica: str, frame: int = -1) -> Optional[pd.DataFrame]:
#         """
#         Extract DSSP data for a specific temperature, replica, and frame.

#         Args:
#             temperature: Temperature (e.g., "320")
#             replica: Replica (e.g., "0")
#             frame: Frame index (default: -1 for last frame)

#         Returns:
#             DataFrame with DSSP data or None if extraction fails
#         """
#         try:
#             with h5py.File(self.h5_path, 'r') as f:
#                 # Check if temperature and replica exist
#                 if temperature not in f[self.domain_id] or replica not in f[self.domain_id][temperature]:
#                     logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
#                     return None

#                 # Extract DSSP data
#                 dssp_data = f[self.domain_id][temperature][replica]['dssp'][frame]

#                 # Extract residue information
#                 resids = f[self.domain_id]['resid'][:]
#                 resnames = [name.decode('utf-8') for name in f[self.domain_id]['resname'][:]]

#                 # Decode DSSP codes
#                 dssp_codes = [code.decode('utf-8') for code in dssp_data]

#                 # Create DataFrame
#                 df = pd.DataFrame({
#                     'domain_id': self.domain_id,
#                     'resid': resids,
#                     'resname': resnames,
#                     'dssp': dssp_codes
#                 })

#                 return df
#         except Exception as e:
#             logging.error(f"Failed to extract DSSP data: {e}")
#             return None

#     def extract_coordinates(self, temperature: str, replica: str, frame: int = -1) -> Optional[Tuple[np.ndarray, List[int], List[str]]]:
#         """
#         Extract coordinate data for a specific temperature, replica, and frame.

#         Args:
#             temperature: Temperature (e.g., "320")
#             replica: Replica (e.g., "0")
#             frame: Frame index (default: -1 for last frame)

#         Returns:
#             Tuple of (coordinates, residue IDs, residue names) or None if extraction fails
#         """
#         try:
#             with h5py.File(self.h5_path, 'r') as f:
#                 # Check if temperature and replica exist
#                 if temperature not in f[self.domain_id] or replica not in f[self.domain_id][temperature]:
#                     logging.warning(f"Temperature {temperature} or replica {replica} not found for domain {self.domain_id}")
#                     return None

#                 # Extract coordinate data for specified frame
#                 coords = f[self.domain_id][temperature][replica]['coords'][frame]

#                 # Extract residue information
#                 resids = f[self.domain_id]['resid'][:].tolist()
#                 resnames = [name.decode('utf-8') for name in f[self.domain_id]['resname'][:]]

#                 return coords, resids, resnames
#         except Exception as e:
#             logging.error(f"Failed to extract coordinate data: {e}")
#             return None

# def process_domains(domain_ids: List[str], data_dir: str, config: Dict[str, Any],
#                     num_cores: int = 1) -> Dict[str, Any]:
#     """
#     Process multiple domains in parallel.

#     Args:
#         domain_ids: List of domain IDs to process
#         data_dir: Directory containing H5 files
#         config: Configuration dictionary
#         num_cores: Number of CPU cores to use

#     Returns:
#         Dictionary with processing results
#     """
#     # Determine number of cores to use
#     max_cores = os.cpu_count() - 2 if os.cpu_count() > 2 else 1
#     n_cores = min(num_cores if num_cores > 0 else max_cores, max_cores)

#     results = {}
#     with ProcessPoolExecutor(max_workers=n_cores) as executor:
#         future_to_domain = {}
#         for domain_id in domain_ids:
#             h5_path = os.path.join(data_dir, f"mdcath_dataset_{domain_id}.h5")
#             if not os.path.exists(h5_path):
#                 logging.warning(f"H5 file not found for domain {domain_id}")
#                 continue

#             future = executor.submit(_process_single_domain, h5_path, config)
#             future_to_domain[future] = domain_id

#         for future in as_completed(future_to_domain):
#             domain_id = future_to_domain[future]
#             try:
#                 result = future.result()
#                 results[domain_id] = result
#             except Exception as e:
#                 logging.error(f"Error processing domain {domain_id}: {e}")
#                 results[domain_id] = {"success": False, "error": str(e)}

#     return results

# def _process_single_domain(h5_path: str, config: Dict[str, Any]) -> Dict[str, Any]:
#     """
#     Process a single domain (helper function for parallel processing).

#     Args:
#         h5_path: Path to H5 file
#         config: Configuration dictionary

#     Returns:
#         Dictionary with processing results
#     """
#     loader = H5DataLoader(h5_path, config)
#     domain_id = loader.domain_id

#     results = {"domain_id": domain_id, "success": False}

#     # Extract RMSF data for all temperatures and replicas
#     temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
#     num_replicas = config.get("num_replicas", 5)

#     rmsf_data = {}
#     for temp in temps:
#         rmsf_data[temp] = {}
#         for r in range(num_replicas):
#             replica = str(r)
#             df = loader.extract_rmsf(temp, replica)
#             if df is not None:
#                 rmsf_data[temp][replica] = df

#     results["rmsf_data"] = rmsf_data

#     # Extract PDB data
#     pdb_str = loader.extract_pdb()
#     if pdb_str:
#         results["pdb_data"] = pdb_str

#     # Extract DSSP data
#     dssp_data = {}
#     for temp in temps:
#         dssp_data[temp] = {}
#         for r in range(num_replicas):
#             replica = str(r)
#             df = loader.extract_dssp(temp, replica)
#             if df is not None:
#                 dssp_data[temp][replica] = df

#     results["dssp_data"] = dssp_data
#     results["success"] = True

#     return results
# EOFMARKER
# echo "Created src/mdcath/core/data_loader.py"

# # ------------------------------------
# # PROCESSING MODULE FILES
# # ------------------------------------

# # Create src/mdcath/processing/rmsf.py
# cat > src/mdcath/processing/rmsf.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Processing module for RMSF data extraction and averaging.
# """

# import os
# import logging
# import numpy as np
# import pandas as pd
# from typing import List, Dict, Optional, Any, Union
# from concurrent.futures import ProcessPoolExecutor

# def calculate_replica_averages(rmsf_data: Dict[str, Dict[str, pd.DataFrame]],
#                               temperature: str) -> Optional[pd.DataFrame]:
#     """
#     Calculate average RMSF across all replicas for a specific temperature.

#     Args:
#         rmsf_data: Dictionary with RMSF data for all replicas
#         temperature: Temperature to calculate average for

#     Returns:
#         DataFrame with average RMSF values or None if calculation fails
#     """
#     try:
#         # Collect all dataframes for this temperature
#         dfs = []
#         for replica, df in rmsf_data.get(temperature, {}).items():
#             if df is not None:
#                 dfs.append(df)

#         if not dfs:
#             logging.warning(f"No RMSF data found for temperature {temperature}")
#             return None

#         # Combine the first dataframe for residue information
#         result_df = dfs[0][['domain_id', 'resid', 'resname']].copy()

#         # Calculate average RMSF
#         rmsf_values = []
#         for df in dfs:
#             rmsf_col = f"rmsf_{temperature}"
#             if rmsf_col in df.columns:
#                 rmsf_values.append(df[rmsf_col].values)

#         if not rmsf_values:
#             logging.warning(f"No RMSF values found for temperature {temperature}")
#             return None

#         # Calculate average
#         avg_rmsf = np.mean(rmsf_values, axis=0)
#         result_df[f"rmsf_{temperature}"] = avg_rmsf

#         return result_df
#     except Exception as e:
#         logging.error(f"Failed to calculate replica averages for temperature {temperature}: {e}")
#         return None

# def calculate_temperature_average(replica_averages: Dict[str, pd.DataFrame]) -> Optional[pd.DataFrame]:
#     """
#     Calculate average RMSF across all temperatures.

#     Args:
#         replica_averages: Dictionary with replica average RMSF data for all temperatures

#     Returns:
#         DataFrame with average RMSF values across all temperatures or None if calculation fails
#     """
#     try:
#         if not replica_averages:
#             logging.warning("No replica averages found")
#             return None

#         # Get the first dataframe for base structure
#         temps = list(replica_averages.keys())
#         first_temp = temps[0]
#         result_df = replica_averages[first_temp][['domain_id', 'resid', 'resname']].copy()

#         # Collect RMSF values for all temperatures
#         rmsf_cols = []
#         for temp, df in replica_averages.items():
#             rmsf_col = f"rmsf_{temp}"
#             if rmsf_col in df.columns:
#                 result_df[rmsf_col] = df[rmsf_col]
#                 rmsf_cols.append(rmsf_col)

#         if not rmsf_cols:
#             logging.warning("No RMSF columns found")
#             return None

#         # Calculate average across all temperatures
#         result_df['rmsf_average'] = result_df[rmsf_cols].mean(axis=1)

#         return result_df
#     except Exception as e:
#         logging.error(f"Failed to calculate temperature average: {e}")
#         return None

# def save_rmsf_data(rmsf_data: Dict[str, Dict[str, pd.DataFrame]],
#                   replica_averages: Dict[str, pd.DataFrame],
#                   temperature_average: pd.DataFrame,
#                   output_dir: str) -> bool:
#     """
#     Save RMSF data to CSV files.

#     Args:
#         rmsf_data: Dictionary with RMSF data for all temperatures and replicas
#         replica_averages: Dictionary with replica average RMSF data for all temperatures
#         temperature_average: DataFrame with average RMSF values across all temperatures
#         output_dir: Directory to save CSV files

#     Returns:
#         Boolean indicating if saving was successful
#     """
#     try:
#         # Create output directory structure
#         os.makedirs(os.path.join(output_dir, "RMSF", "replicas"), exist_ok=True)
#         os.makedirs(os.path.join(output_dir, "RMSF", "replica_average", "average"), exist_ok=True)

#         # Save replica data
#         for temp, replicas in rmsf_data.items():
#             for replica, df in replicas.items():
#                 replica_dir = os.path.join(output_dir, "RMSF", "replicas", f"replica_{replica}", temp)
#                 os.makedirs(replica_dir, exist_ok=True)

#                 output_file = os.path.join(replica_dir, f"rmsf_replica{replica}_temperature{temp}.csv")
#                 df.to_csv(output_file, index=False)
#                 logging.info(f"Saved RMSF data to {output_file}")

#         # Save replica averages
#         for temp, df in replica_averages.items():
#             temp_dir = os.path.join(output_dir, "RMSF", "replica_average", temp)
#             os.makedirs(temp_dir, exist_ok=True)

#             output_file = os.path.join(temp_dir, f"rmsf_replica_average_temperature{temp}.csv")
#             df.to_csv(output_file, index=False)
#             logging.info(f"Saved replica average RMSF data to {output_file}")

#         # Save temperature average
#         output_file = os.path.join(output_dir, "RMSF", "replica_average", "average",
#                                   "rmsf_all_temperatures_all_replicas.csv")
#         temperature_average.to_csv(output_file, index=False)
#         logging.info(f"Saved temperature average RMSF data to {output_file}")

#         return True
#     except Exception as e:
#         logging.error(f"Failed to save RMSF data: {e}")
#         return False

# def process_rmsf_data(domain_results: Dict[str, Dict[str, Any]], config: Dict[str, Any]) -> Dict[str, Any]:
#     """
#     Process RMSF data for all domains.

#     Args:
#         domain_results: Dictionary with processing results for all domains
#         config: Configuration dictionary

#     Returns:
#         Dictionary with RMSF processing results
#     """
#     temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
#     output_dir = config.get("output", {}).get("base_dir", "./outputs")

#     # Combine RMSF data from all domains
#     combined_rmsf_data = {temp: {} for temp in temps}
#     for domain_id, result in domain_results.items():
#         if not result.get("success", False):
#             continue

#         rmsf_data = result.get("rmsf_data", {})
#         for temp in temps:
#             if temp in rmsf_data:
#                 for replica, df in rmsf_data[temp].items():
#                     if replica not in combined_rmsf_data[temp]:
#                         combined_rmsf_data[temp][replica] = []
#                     combined_rmsf_data[temp][replica].append(df)

#     # Concatenate dataframes for each temperature and replica
#     for temp in combined_rmsf_data:
#         for replica in combined_rmsf_data[temp]:
#             if combined_rmsf_data[temp][replica]:
#                 combined_rmsf_data[temp][replica] = pd.concat(combined_rmsf_data[temp][replica], ignore_index=True)

#     # Calculate replica averages
#     replica_averages = {}
#     for temp in temps:
#         avg_df = calculate_replica_averages(combined_rmsf_data, temp)
#         if avg_df is not None:
#             replica_averages[temp] = avg_df

#     # Calculate temperature average
#     temperature_average = calculate_temperature_average(replica_averages)

#     # Save RMSF data
#     save_success = save_rmsf_data(combined_rmsf_data, replica_averages, temperature_average, output_dir)

#     return {
#         "combined_rmsf_data": combined_rmsf_data,
#         "replica_averages": replica_averages,
#         "temperature_average": temperature_average,
#         "save_success": save_success
#     }
# EOFMARKER
# echo "Created src/mdcath/processing/rmsf.py"

# # Create src/mdcath/processing/pdb.py
# cat > src/mdcath/processing/pdb.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Processing module for PDB data extraction and cleaning.
# """

# import os
# import logging
# import numpy as np
# from typing import Dict, Any, Optional, List, Tuple
# from concurrent.futures import ProcessPoolExecutor, as_completed

# def save_pdb_file(pdb_string: str, output_path: str, config: Dict[str, Any]) -> bool:
#     """
#     Save a PDB string to a file with cleaning applied.

#     Args:
#         pdb_string: PDB data as a string
#         output_path: Path to save the cleaned PDB
#         config: Configuration dictionary

#     Returns:
#         Boolean indicating if saving was successful
#     """
#     try:
#         # Write original PDB to temporary file
#         temp_path = output_path + ".temp"
#         with open(temp_path, 'w') as f:
#             f.write(pdb_string)

#         # Clean the PDB file
#         success = fix_pdb(temp_path, output_path, config)

#         # Remove temporary file
#         if os.path.exists(temp_path):
#             os.remove(temp_path)

#         return success
#     except Exception as e:
#         logging.error(f"Failed to save PDB file: {e}")
#         return False

# def fix_pdb(input_pdb_path: str, output_pdb_path: str, config: Dict[str, Any]) -> bool:
#     """
#     Clean and fix a PDB file for downstream processing.

#     Args:
#         input_pdb_path: Path to input PDB file
#         output_pdb_path: Path to save the cleaned PDB file
#         config: Configuration dictionary

#     Returns:
#         Boolean indicating if cleaning was successful
#     """
#     if not os.path.isfile(input_pdb_path):
#         logging.error(f"PDB file not found: {input_pdb_path}")
#         return False

#     try:
#         with open(input_pdb_path, 'r') as f:
#             lines = f.readlines()

#         # Check for and add CRYST1 if needed
#         has_cryst1 = False
#         for line in lines:
#             if line.strip() and line.startswith("CRYST1"):
#                 has_cryst1 = True
#                 break

#         cleaned_lines = []
#         if not has_cryst1 and config.get("pdb_cleaning", {}).get("add_cryst1_record", True):
#             cleaned_lines.append("CRYST1          1\n")

#         # Process atom records
#         for line in lines:
#             if line.startswith("ATOM") or line.startswith("HETATM"):
#                 try:
#                     # Extract fields
#                     record_type = line[0:6].strip()
#                     atom_num = int(line[6:11].strip())
#                     atom_name = line[12:16].strip()
#                     alt_loc = line[16:17].strip()
#                     res_name = line[17:20].strip()
#                     chain_id = line[21:22].strip()
#                     res_num = line[22:26].strip()
#                     ins_code = line[26:27].strip()
#                     x = float(line[30:38].strip())
#                     y = float(line[38:46].strip())
#                     z = float(line[46:54].strip())
#                     occupancy = line[54:60].strip() if line[54:60].strip() else "0.00"
#                     temp_factor = line[60:66].strip() if line[60:66].strip() else "0.00"
#                     element = line[76:78].strip() if len(line) >= 78 else ""

#                     # Apply cleaning rules
#                     clean_config = config.get("pdb_cleaning", {})
#                     if clean_config.get("replace_chain_0_with_A", True) and chain_id == "0":
#                         chain_id = "A"

#                     if clean_config.get("correct_unusual_residue_names", True):
#                         if res_name in ["HSD", "HSE", "HSP"]:
#                             res_name = "HIS"

#                     # Skip hydrogens if configured
#                     if clean_config.get("remove_hydrogens", False) and element == "H":
#                         continue

#                     # Format the line according to PDB standard
#                     new_line = f"ATOM  {atom_num:5d} {atom_name:<4s}{alt_loc:1s}{res_name:3s} {chain_id:1s}{res_num:4s}{ins_code:1s}   {x:8.3f}{y:8.3f}{z:8.3f}{float(occupancy):6.2f}{float(temp_factor):6.2f}           {element:>2s}  \n"
#                     cleaned_lines.append(new_line)
#                 except Exception as e:
#                     logging.warning(f"Error processing line: {line.strip()} - {e}")
#                     # Keep the original line in case of parsing errors
#                     cleaned_lines.append(line)
#             else:
#                 # Keep non-ATOM lines as they are
#                 cleaned_lines.append(line)

#         # Write cleaned PDB
#         with open(output_pdb_path, 'w') as f:
#             f.writelines(cleaned_lines)

#         return True
#     except Exception as e:
#         logging.error(f"Failed to clean PDB {input_pdb_path}: {e}")
#         return False

# def extract_frames(coords: np.ndarray, resids: List[int], resnames: List[str],
#                  domain_id: str, output_dir: str, temperature: str, replica: str,
#                  config: Dict[str, Any]) -> bool:
#     """
#     Extract frames from coordinate data and save as PDB files.

#     Args:
#         coords: Coordinate data
#         resids: Residue IDs
#         resnames: Residue names
#         domain_id: Domain identifier
#         output_dir: Directory to save frame PDBs
#         temperature: Temperature
#         replica: Replica
#         config: Configuration dictionary

#     Returns:
#         Boolean indicating if extraction was successful
#     """
#     frame_selection = config.get("processing", {}).get("frame_selection", {})
#     method = frame_selection.get("method", "rmsd")
#     num_frames = frame_selection.get("num_frames", 1)

#     try:
#         # Create output directory
#         frame_dir = os.path.join(output_dir, "frames", f"replica_{replica}", temperature)
#         os.makedirs(frame_dir, exist_ok=True)

#         # For simplicity, extract the last frame if num_frames is 1
#         if num_frames == 1:
#             frame_idx = -1
#             frame_coords = coords[frame_idx]

#             # Create PDB string for the frame
#             pdb_lines = ["CRYST1          1\n"]
#             atom_num = 1
#             for i, (resid, resname) in enumerate(zip(resids, resnames)):
#                 x, y, z = frame_coords[i]
#                 atom_name = "CA"  # Simplification: assume each residue has one CA atom
#                 chain_id = "A"
#                 pdb_lines.append(f"ATOM  {atom_num:5d} {atom_name:<4s} {resname:3s} {chain_id:1s}{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C  \n")
#                 atom_num += 1

#             # Write frame PDB
#             frame_path = os.path.join(frame_dir, f"{domain_id}_frame.pdb")
#             with open(frame_path, 'w') as f:
#                 f.writelines(pdb_lines)

#             logging.info(f"Extracted frame for domain {domain_id}, temperature {temperature}, replica {replica}")
#             return True
#         else:
#             logging.warning(f"Multiple frame extraction not implemented yet")
#             return False
#     except Exception as e:
#         logging.error(f"Failed to extract frames: {e}")
#         return False

# def process_pdb_data(domain_results: Dict[str, Dict[str, Any]], config: Dict[str, Any]) -> Dict[str, Any]:
#     """
#     Process PDB data for all domains.

#     Args:
#         domain_results: Dictionary with processing results for all domains
#         config: Configuration dictionary

#     Returns:
#         Dictionary with PDB processing results
#     """
#     output_dir = config.get("output", {}).get("base_dir", "./outputs")

#     # Create output directories
#     pdb_dir = os.path.join(output_dir, "pdbs")
#     os.makedirs(pdb_dir, exist_ok=True)

#     results = {}

#     for domain_id, result in domain_results.items():
#         if not result.get("success", False):
#             continue

#         # Save cleaned PDB
#         pdb_data = result.get("pdb_data")
#         if pdb_data:
#             pdb_path = os.path.join(pdb_dir, f"{domain_id}.pdb")
#             success = save_pdb_file(pdb_data, pdb_path, config)
#             if success:
#                 logging.info(f"Saved cleaned PDB for domain {domain_id}")
#                 results[domain_id] = {"pdb_saved": True, "pdb_path": pdb_path}
#             else:
#                 logging.error(f"Failed to save cleaned PDB for domain {domain_id}")
#                 results[domain_id] = {"pdb_saved": False}

#     return results
# EOFMARKER
# echo "Created src/mdcath/processing/pdb.py"

# # Create src/mdcath/processing/core_exterior.py
# cat > src/mdcath/processing/core_exterior.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Processing module for core/exterior classification.
# """

# import os
# import logging
# import subprocess
# import tempfile
# import pandas as pd
# import numpy as np
# from typing import Dict, Any, Optional, List, Tuple

# def compute_core_exterior(pdb_file: str, config: Dict[str, Any]) -> Optional[pd.DataFrame]:
#     """
#     Classify residues as 'core' or 'exterior' based on solvent accessibility.

#     Args:
#         pdb_file: Path to the cleaned PDB file
#         config: Configuration dictionary

#     Returns:
#         DataFrame with columns 'resid' and 'core_exterior' or None if classification fails
#     """
#     method = config.get("core_exterior", {}).get("method", "msms")

#     if method == "msms":
#         return compute_core_exterior_msms(pdb_file, config)
#     else:
#         return compute_core_exterior_biopython(pdb_file, config)

# def compute_core_exterior_msms(pdb_file: str, config: Dict[str, Any]) -> Optional[pd.DataFrame]:
#     """
#     Use MSMS to classify residues as 'core' or 'exterior'.

#     Args:
#         pdb_file: Path to the cleaned PDB file
#         config: Configuration dictionary

#     Returns:
#         DataFrame with columns 'resid' and 'core_exterior' or None if MSMS fails
#     """
#     msms_dir = config.get("core_exterior", {}).get("msms_executable_dir", "./msms")
#     ses_threshold = config.get("core_exterior", {}).get("ses_threshold", 1.0)
#     protein_name = os.path.basename(pdb_file).split('.')[0]

#     try:
#         # Create temporary directory for MSMS files
#         with tempfile.TemporaryDirectory() as tmp_dir:
#             # Paths to MSMS executables and output files
#             pdb2xyzr_exe = os.path.join(msms_dir, "pdb_to_xyzr")
#             msms_exe = os.path.join(msms_dir, "msms.x86_64Linux2.2.6.1")
#             xyzr_file = os.path.join(tmp_dir, f"{protein_name}.xyzr")
#             area_base = os.path.join(tmp_dir, f"{protein_name}")
#             area_file = os.path.join(tmp_dir, f"{protein_name}.area")

#             # Check MSMS executables
#             if not os.path.exists(pdb2xyzr_exe) or not os.path.exists(msms_exe):
#                 logging.warning(f"MSMS executables not found in {msms_dir}, falling back to Biopython")
#                 return compute_core_exterior_biopython(pdb_file, config)

#             # Run pdb_to_xyzr
#             cmd_xyzr = f"{pdb2xyzr_exe} {pdb_file} > {xyzr_file}"
#             result = subprocess.run(cmd_xyzr, shell=True, check=False,
#                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)

#             if result.returncode != 0 or not os.path.exists(xyzr_file) or os.path.getsize(xyzr_file) == 0:
#                 logging.warning(f"pdb_to_xyzr failed: {result.stderr.decode()}, falling back to Biopython")
#                 return compute_core_exterior_biopython(pdb_file, config)

#             # Run MSMS
#             cmd_msms = f"{msms_exe} -if {xyzr_file} -af {area_base}"
#             result = subprocess.run(cmd_msms, shell=True, check=False,
#                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)

#             if result.returncode != 0 or not os.path.exists(area_file):
#                 logging.warning(f"MSMS failed: {result.stderr.decode()}, falling back to Biopython")
#                 return compute_core_exterior_biopython(pdb_file, config)

#             # Parse atom-level PDB data
#             per_atom_df = parse_pdb_atoms(pdb_file)
#             if per_atom_df.empty:
#                 logging.warning(f"Failed to parse atoms from PDB, falling back to Biopython")
#                 return compute_core_exterior_biopython(pdb_file, config)

#             # Parse MSMS area file
#             area_df = parse_area_file(area_file)
#             if area_df.empty:
#                 logging.warning(f"Failed to parse area file, falling back to Biopython")
#                 return compute_core_exterior_biopython(pdb_file, config)

#             # Combine atom data with MSMS results
#             if len(area_df) != len(per_atom_df):
#                 logging.warning(f"Atom count mismatch: {len(area_df)} vs {len(per_atom_df)}, falling back to Biopython")
#                 return compute_core_exterior_biopython(pdb_file, config)

#             # Merge data
#             per_atom_df = pd.concat([per_atom_df.reset_index(drop=True),
#                                     area_df.reset_index(drop=True)], axis=1)

#             # Calculate mean SES per residue
#             mean_ses_per_res = per_atom_df.groupby("resid")["SES"].mean()

#             # Classify residues as core or exterior
#             exterior_residues = mean_ses_per_res[mean_ses_per_res > ses_threshold].index
#             resids = mean_ses_per_res.index.tolist()
#             core_exterior = ["exterior" if r in exterior_residues else "core" for r in resids]

#             # Create final dataframe
#             result_df = pd.DataFrame({
#                 "resid": resids,
#                 "core_exterior": core_exterior
#             })

#             return result_df
#     except Exception as e:
#         logging.warning(f"MSMS processing failed: {e}, falling back to Biopython")
#         return compute_core_exterior_biopython(pdb_file, config)

# def compute_core_exterior_biopython(pdb_file: str, config: Dict[str, Any]) -> pd.DataFrame:
#     """
#     Use Biopython's SASA calculation to classify residues as 'core' or 'exterior'.

#     Args:
#         pdb_file: Path to the cleaned PDB file
#         config: Configuration dictionary

#     Returns:
#         DataFrame with columns 'resid' and 'core_exterior'
#     """
#     try:
#         from Bio.PDB import PDBParser, Selection
#         from Bio.PDB.SASA import ShrakeRupley

#         # Set SASA threshold
#         sasa_threshold = config.get("core_exterior", {}).get("sasa_threshold", 20.0)

#         # Parse PDB
#         parser = PDBParser(QUIET=True)
#         structure = parser.get_structure("protein", pdb_file)
#         model = structure[0]

#         # Calculate SASA
#         sr = ShrakeRupley()
#         sr.compute(model, level="R")  # Compute at residue level

#         # Extract results
#         results = []
#         for chain in model:
#             for residue in chain:
#                 if residue.id[0] == " ":  # Standard residue
#                     resid = residue.id[1]
#                     sasa = residue.sasa
#                     core_exterior = "exterior" if sasa > sasa_threshold else "core"
#                     results.append({"resid": resid, "core_exterior": core_exterior})

#         return pd.DataFrame(results)
#     except Exception as e:
#         logging.error(f"Biopython SASA calculation failed: {e}")
#         return fallback_core_exterior(pdb_file)

# def fallback_core_exterior(pdb_file: str) -> pd.DataFrame:
#     """
#     Fallback method to classify residues when other methods fail.
#     Classifies outer 1/3 of residues as exterior, inner 2/3 as core.

#     Args:
#         pdb_file: Path to the cleaned PDB file

#     Returns:
#         DataFrame with columns 'resid' and 'core_exterior'
#     """
#     try:
#         # Parse PDB to get residue information
#         residue_df = parse_pdb_residues(pdb_file)
#         if residue_df.empty:
#             # Create empty DataFrame with required columns
#             return pd.DataFrame(columns=["resid", "core_exterior"])

#         # Sort by residue ID
#         residue_df = residue_df.sort_values("resid")

#         # Simple classification: outer 1/3 of residues as exterior, inner 2/3 as core
#         total_residues = len(residue_df)
#         boundary = int(total_residues * 2/3)

#         residue_df["core_exterior"] = ["core"] * total_residues
#         residue_df.loc[boundary:, "core_exterior"] = "exterior"

#         return residue_df[["resid", "core_exterior"]]
#     except Exception as e:
#         logging.error(f"Fallback classification failed: {e}")
#         return pd.DataFrame(columns=["resid", "core_exterior"])

# def parse_pdb_residues(pdb_file: str) -> pd.DataFrame:
#     """
#     Parse a PDB file to extract residue-level information.

#     Args:
#         pdb_file: Path to the PDB file

#     Returns:
#         DataFrame with residue information
#     """
#     try:
#         from Bio.PDB import PDBParser

#         parser = PDBParser(QUIET=True)
#         structure = parser.get_structure("protein", pdb_file)

#         records = []
#         for model in structure:
#             for chain in model:
#                 chain_id = chain.id
#                 for residue in chain:
#                     if residue.id[0] == " ":  # Standard residue
#                         records.append({
#                             "resid": residue.id[1],
#                             "resname": residue.get_resname(),
#                             "chain": chain_id
#                         })

#         return pd.DataFrame(records)
#     except Exception as e:
#         logging.error(f"Failed to parse PDB residues: {e}")
#         return pd.DataFrame()

# def parse_pdb_atoms(pdb_file: str) -> pd.DataFrame:
#     """
#     Parse a PDB file to extract atom-level information.

#     Args:
#         pdb_file: Path to the PDB file

#     Returns:
#         DataFrame with atom information
#     """
#     try:
#         from Bio.PDB import PDBParser

#         parser = PDBParser(QUIET=True)
#         structure = parser.get_structure("protein", pdb_file)

#         records = []
#         atom_idx = 0
#         for model in structure:
#             for chain in model:
#                 for residue in chain:
#                     if residue.id[0] == " ":  # Standard residue
#                         res_id = residue.id[1]
#                         res_name = residue.get_resname()
#                         for atom in residue:
#                             atom_idx += 1
#                             records.append({
#                                 "atom_idx": atom_idx,
#                                 "resid": res_id,
#                                 "resname": res_name,
#                                 "atom_name": atom.get_name()
#                             })

#         return pd.DataFrame(records)
#     except Exception as e:
#         logging.error(f"Failed to parse PDB atoms: {e}")
#         return pd.DataFrame()

# def parse_area_file(area_file: str) -> pd.DataFrame:
#     """
#     Parse an MSMS .area file to extract SES values per atom.

#     Args:
#         area_file: Path to the MSMS .area file

#     Returns:
#         DataFrame with SES values
#     """
#     try:
#         atom_idx = []
#         ses = []

#         with open(area_file, "r") as f:
#             for line in f:
#                 if "Atom" in line or not line.strip():
#                     continue

#                 cols = line.split()
#                 if len(cols) >= 2:
#                     atom_idx.append(int(cols[0]))
#                     ses.append(float(cols[1]))

#         return pd.DataFrame({"atom_idx": atom_idx, "SES": ses})
#     except Exception as e:
#         logging.error(f"Failed to parse area file: {e}")
#         return pd.DataFrame()
# EOFMARKER
# echo "Created src/mdcath/processing/core_exterior.py"

# # Create src/mdcath/processing/features.py
# cat > src/mdcath/processing/features.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Processing module for generating ML features.
# """

# import os
# import logging
# import numpy as np
# import pandas as pd
# from typing import Dict, Any, Optional, List, Tuple, Union

# def generate_ml_features(rmsf_data: Dict[str, pd.DataFrame],
#                        core_exterior_data: Dict[str, pd.DataFrame],
#                        dssp_data: Dict[str, Dict[str, pd.DataFrame]],
#                        config: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
#     """
#     Generate ML features for all domains.

#     Args:
#         rmsf_data: Dictionary with RMSF data for all temperatures
#         core_exterior_data: Dictionary with core/exterior classification
#         dssp_data: Dictionary with DSSP data for all temperatures and replicas
#         config: Configuration dictionary

#     Returns:
#         Dictionary with ML feature dataframes for all temperatures
#     """
#     try:
#         # Get list of all domains
#         domain_ids = set()
#         for temp, df in rmsf_data.items():
#             domain_ids.update(df["domain_id"].unique())

#         domain_ids = list(domain_ids)
#         logging.info(f"Generating ML features for {len(domain_ids)} domains")

#         # Create feature dataframes for each temperature
#         temps = [t for t in rmsf_data.keys() if t != "average"]
#         feature_dfs = {}

#         for temp in temps:
#             # Start with RMSF data
#             if temp not in rmsf_data:
#                 logging.warning(f"RMSF data not found for temperature {temp}")
#                 continue

#             df = rmsf_data[temp].copy()

#             # Add protein size
#             df["protein_size"] = df.groupby("domain_id")["resid"].transform("count")

#             # Add normalized residue position
#             df["normalized_resid"] = df.groupby("domain_id")["resid"].transform(
#                 lambda x: (x - x.min()) / (x.max() - x.min())
#             )

#             # Add core/exterior classification
#             for domain_id in df["domain_id"].unique():
#                 if domain_id in core_exterior_data:
#                     core_ext_df = core_exterior_data[domain_id]

#                     # Merge core/exterior data
#                     domain_mask = df["domain_id"] == domain_id
#                     df_domain = df[domain_mask].copy()

#                     # Reset index for proper merging
#                     df_domain = df_domain.reset_index(drop=True)
#                     df_domain = pd.merge(df_domain, core_ext_df, on="resid", how="left")

#                     # Update the main dataframe
#                     df.loc[domain_mask] = df_domain

#             # Fill missing core/exterior values with 'unknown'
#             if "core_exterior" in df.columns:
#                 df["core_exterior"] = df["core_exterior"].fillna("unknown")
#             else:
#                 df["core_exterior"] = "unknown"

#             # Add DSSP data
#             if temp in dssp_data and "0" in dssp_data[temp]:
#                 for domain_id in df["domain_id"].unique():
#                     domain_dssp = None

#                     # Find DSSP data for this domain
#                     for replica, dssp_df in dssp_data[temp].items():
#                         domain_dssp_subset = dssp_df[dssp_df["domain_id"] == domain_id]
#                         if not domain_dssp_subset.empty:
#                             domain_dssp = domain_dssp_subset[["resid", "dssp"]]
#                             break

#                     if domain_dssp is not None:
#                         # Merge DSSP data
#                         domain_mask = df["domain_id"] == domain_id
#                         df_domain = df[domain_mask].copy()
#                         df_domain = df_domain.reset_index(drop=True)
#                         df_domain = pd.merge(df_domain, domain_dssp, on="resid", how="left")

#                         # Update the main dataframe
#                         df.loc[domain_mask] = df_domain

#             # Fill missing DSSP values with 'C' (coil)
#             if "dssp" in df.columns:
#                 df["dssp"] = df["dssp"].fillna("C")
#             else:
#                 df["dssp"] = "C"

#             # Encode categorical variables
#             # Resname encoding
#             unique_resnames = df["resname"].unique()
#             resname_mapping = {name: i for i, name in enumerate(sorted(unique_resnames))}
#             df["resname_encoded"] = df["resname"].map(resname_mapping)

#             # Core/exterior encoding
#             core_ext_mapping = {"core": 1, "exterior": 2, "unknown": 0}
#             df["core_exterior_encoded"] = df["core_exterior"].map(core_ext_mapping)

#             # Secondary structure encoding
#             # Simplified 3-state encoding: Helix (H,G,I), Sheet (E,B), Coil (others)
#             def encode_ss(ss):
#                 if ss in ["H", "G", "I"]:
#                     return 0  # Helix
#                 elif ss in ["E", "B"]:
#                     return 1  # Sheet
#                 else:
#                     return 2  # Coil

#             df["secondary_structure_encoded"] = df["dssp"].apply(encode_ss)

#             # Add dummy relative accessibility (placeholder for more advanced calculation)
#             df["relative_accessibility"] = 1.0

#             # Reorder columns to put domain_id first
#             cols = df.columns.tolist()
#             cols.remove("domain_id")
#             cols = ["domain_id"] + cols
#             df = df[cols]

#             # Store the feature dataframe
#             feature_dfs[temp] = df

#         # Calculate average features
#         if temps:
#             # Start with a copy of the first temperature's features
#             avg_df = feature_dfs[temps[0]].copy()

#             # Collect RMSF columns
#             rmsf_cols = [f"rmsf_{temp}" for temp in temps]
#             rmsf_vals = []

#             for temp in temps:
#                 if temp in feature_dfs and f"rmsf_{temp}" in feature_dfs[temp].columns:
#                     temp_df = feature_dfs[temp]
#                     # Extract the domain, resid, and RMSF columns
#                     temp_subset = temp_df[["domain_id", "resid", f"rmsf_{temp}"]].copy()
#                     rmsf_vals.append(temp_subset)

#             # Merge all RMSF values
#             if rmsf_vals:
#                 merged_df = rmsf_vals[0]
#                 for i in range(1, len(rmsf_vals)):
#                     merged_df = pd.merge(merged_df, rmsf_vals[i], on=["domain_id", "resid"])

#                 # Calculate average RMSF
#                 merged_df["rmsf_average"] = merged_df[[f"rmsf_{temp}" for temp in temps]].mean(axis=1)

#                 # Replace the RMSF columns in the average dataframe
#                 avg_df = pd.merge(avg_df.drop(rmsf_cols, axis=1, errors="ignore"),
#                                 merged_df[["domain_id", "resid", "rmsf_average"]],
#                                 on=["domain_id", "resid"])

#             feature_dfs["average"] = avg_df

#         return feature_dfs
#     except Exception as e:
#         logging.error(f"Failed to generate ML features: {e}")
#         return {}

# def save_ml_features(feature_dfs: Dict[str, pd.DataFrame], output_dir: str) -> bool:
#     """
#     Save ML features to CSV files.

#     Args:
#         feature_dfs: Dictionary with ML feature dataframes
#         output_dir: Directory to save CSV files

#     Returns:
#         Boolean indicating if saving was successful
#     """
#     try:
#         # Create output directory
#         os.makedirs(output_dir, exist_ok=True)

#         # Save each feature dataframe
#         for temp, df in feature_dfs.items():
#             if temp == "average":
#                 output_file = os.path.join(output_dir, "final_dataset_temperature_average.csv")
#             else:
#                 output_file = os.path.join(output_dir, f"final_dataset_temperature_{temp}.csv")

#             df.to_csv(output_file, index=False)
#             logging.info(f"Saved ML features to {output_file}")

#         return True
#     except Exception as e:
#         logging.error(f"Failed to save ML features: {e}")
#         return False

# def process_ml_features(rmsf_results: Dict[str, Any],
#                        pdb_results: Dict[str, Any],
#                        domain_results: Dict[str, Dict[str, Any]],
#                        config: Dict[str, Any]) -> Dict[str, Any]:
#     """
#     Process ML features for all domains.

#     Args:
#         rmsf_results: Dictionary with RMSF processing results
#         pdb_results: Dictionary with PDB processing results
#         domain_results: Dictionary with processing results for all domains
#         config: Configuration dictionary

#     Returns:
#         Dictionary with ML feature processing results
#     """
#     output_dir = config.get("output", {}).get("base_dir", "./outputs")

#     # Extract RMSF data
#     replica_averages = rmsf_results.get("replica_averages", {})
#     temperature_average = rmsf_results.get("temperature_average")

#     if not replica_averages:
#         logging.error("No RMSF data available for ML feature generation")
#         return {"success": False, "error": "No RMSF data available"}

#     # Create dictionary with all RMSF data
#     rmsf_data = replica_averages.copy()
#     if temperature_average is not None:
#         rmsf_data["average"] = temperature_average

#     # Extract core/exterior data
#     core_exterior_data = {}

#     # Collect DSSP data
#     dssp_data = {}
#     temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]

#     for domain_id, result in domain_results.items():
#         if not result.get("success", False):
#             continue

#         # Extract DSSP data
#         domain_dssp = result.get("dssp_data", {})
#         for temp in temps:
#             if temp in domain_dssp:
#                 if temp not in dssp_data:
#                     dssp_data[temp] = {}

#                 for replica, df in domain_dssp[temp].items():
#                     if replica not in dssp_data[temp]:
#                         dssp_data[temp][replica] = []

#                     dssp_data[temp][replica].append(df)

#     # Concatenate DSSP dataframes
#     for temp in dssp_data:
#         for replica in dssp_data[temp]:
#             if dssp_data[temp][replica]:
#                 dssp_data[temp][replica] = pd.concat(dssp_data[temp][replica], ignore_index=True)

#     # Generate and save ML features
#     feature_dfs = generate_ml_features(rmsf_data, core_exterior_data, dssp_data, config)

#     if not feature_dfs:
#         logging.error("Failed to generate ML features")
#         return {"success": False, "error": "Feature generation failed"}

#     # Save ML features
#     ml_dir = os.path.join(output_dir, "ML_features")
#     save_success = save_ml_features(feature_dfs, ml_dir)

#     return {
#         "success": save_success,
#         "feature_dfs": feature_dfs,
#         "output_dir": ml_dir
#     }
# EOFMARKER
# echo "Created src/mdcath/processing/features.py"

# # Create src/mdcath/processing/voxelizer.py
# cat > src/mdcath/processing/voxelizer.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Processing module for voxelizing protein structures using aposteriori.
# """

# import os
# import logging
# import subprocess
# from typing import Dict, Any, Optional, List, Tuple
# from concurrent.futures import ProcessPoolExecutor, as_completed

# def voxelize_domain(pdb_file: str, output_dir: str, config: Dict[str, Any]) -> Optional[str]:
#     """
#     Voxelize a cleaned PDB file using aposteriori's make-frame-dataset command.

#     Args:
#         pdb_file: Path to the cleaned PDB file
#         output_dir: Directory to save voxelized output
#         config: Configuration dictionary with voxelization parameters

#     Returns:
#         Path to the output file if successful, None otherwise
#     """
#     try:
#         # Get voxelization parameters from config
#         voxel_config = config.get("processing", {}).get("voxelization", {})
#         frame_edge_length = voxel_config.get("frame_edge_length", 12.0)
#         voxels_per_side = voxel_config.get("voxels_per_side", 21)
#         atom_encoder = voxel_config.get("atom_encoder", "CNOCBCA")
#         encode_cb = voxel_config.get("encode_cb", True)
#         compression_gzip = voxel_config.get("compression_gzip", True)
#         voxelise_all_states = voxel_config.get("voxelise_all_states", False)

#         # Create output directory
#         os.makedirs(output_dir, exist_ok=True)

#         # Create output name
#         output_name = os.path.basename(pdb_file).split('.')[0] + "_voxelized"

#         # Build aposteriori command
#         cmd = [
#             "make-frame-dataset",
#             "-o", output_dir,
#             "-n", output_name,
#             "-v",
#             "--frame-edge-length", str(frame_edge_length),
#             "--voxels-per-side", str(voxels_per_side),
#             "-ae", atom_encoder,
#             "-cb", str(encode_cb).lower(),
#             "-comp", str(compression_gzip).lower(),
#             "-vas", str(voxelise_all_states).lower(),
#             pdb_file
#         ]

#         # Run the command
#         logging.info(f"Running aposteriori voxelization for {pdb_file}")
#         result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

#         # Check for expected output file
#         output_file = os.path.join(output_dir, f"{output_name}.hdf5")
#         if os.path.exists(output_file):
#             logging.info(f"Voxelization completed: {output_file}")
#             return output_file
#         else:
#             logging.error(f"Voxelization failed: output file not found")
#             logging.debug(f"Command output: {result.stdout}")
#             logging.debug(f"Command error: {result.stderr}")
#             return None
#     except subprocess.CalledProcessError as e:
#         logging.error(f"Voxelization failed: {e.stderr}")
#         return None
#     except Exception as e:
#         logging.error(f"Voxelization failed with unexpected error: {e}")
#         return None

# def voxelize_domains(pdb_results: Dict[str, Any], config: Dict[str, Any]) -> Dict[str, Any]:
#     """
#     Voxelize multiple domains in parallel.

#     Args:
#         pdb_results: Dictionary with PDB processing results
#         config: Configuration dictionary

#     Returns:
#         Dictionary with voxelization results
#     """
#     output_dir = config.get("output", {}).get("base_dir", "./outputs")
#     voxel_dir = os.path.join(output_dir, "voxelized")
#     os.makedirs(voxel_dir, exist_ok=True)

#     # Determine number of cores to use
#     performance_config = config.get("performance", {})
#     num_cores = performance_config.get("num_cores", 0)
#     max_cores = os.cpu_count() - 2 if os.cpu_count() > 2 else 1
#     n_cores = min(num_cores if num_cores > 0 else max_cores, max_cores)

#     results = {}
#     with ProcessPoolExecutor(max_workers=n_cores) as executor:
#         future_to_domain = {}

#         for domain_id, result in pdb_results.items():
#             if not result.get("pdb_saved", False):
#                 continue

#             pdb_path = result.get("pdb_path")
#             if not pdb_path or not os.path.exists(pdb_path):
#                 logging.warning(f"PDB file not found for domain {domain_id}")
#                 continue

#             future = executor.submit(voxelize_domain, pdb_path, voxel_dir, config)
#             future_to_domain[future] = domain_id

#         for future in as_completed(future_to_domain):
#             domain_id = future_to_domain[future]
#             try:
#                 output_file = future.result()
#                 results[domain_id] = {
#                     "success": output_file is not None,
#                     "output_file": output_file
#                 }
#             except Exception as e:
#                 logging.error(f"Error voxelizing domain {domain_id}: {e}")
#                 results[domain_id] = {"success": False, "error": str(e)}

#     return results
# EOFMARKER
# echo "Created src/mdcath/processing/voxelizer.py"

# # Create src/mdcath/processing/visualization.py
# cat > src/mdcath/processing/visualization.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Module for generating visualizations of processed mdCATH data.
# """

# import os
# import logging
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# from typing import Dict, Any, Optional, List, Tuple

# def create_temperature_summary_heatmap(rmsf_data: Dict[str, pd.DataFrame], 
#                                      output_dir: str) -> Optional[str]:
#     """
#     Create a heatmap showing RMSF values across temperatures for all domains.
    
#     Args:
#         rmsf_data: Dictionary with RMSF data for all temperatures
#         output_dir: Directory to save visualization
        
#     Returns:
#         Path to the saved figure or None if creation fails
#     """
#     try:
#         # Ensure output directory exists
#         vis_dir = os.path.join(output_dir, "visualizations")
#         os.makedirs(vis_dir, exist_ok=True)
        
#         # Extract temperature values
#         temps = [temp for temp in rmsf_data.keys() if temp != "average"]
        
#         if not temps:
#             logging.warning("No temperature data available for heatmap")
#             return None
            
#         # Prepare data for heatmap
#         domain_ids = set()
#         for temp in temps:
#             if temp in rmsf_data:
#                 domain_ids.update(rmsf_data[temp]["domain_id"].unique())
        
#         domain_ids = sorted(list(domain_ids))
        
#         # Create a dataframe for the heatmap
#         heatmap_data = []
        
#         for domain_id in domain_ids:
#             domain_data = {"domain_id": domain_id}
            
#             for temp in temps:
#                 if temp in rmsf_data:
#                     domain_temp_data = rmsf_data[temp][rmsf_data[temp]["domain_id"] == domain_id]
#                     if not domain_temp_data.empty:
#                         domain_data[f"rmsf_{temp}"] = domain_temp_data[f"rmsf_{temp}"].mean()
            
#             heatmap_data.append(domain_data)
        
#         if not heatmap_data:
#             logging.warning("No data available for heatmap")
#             return None
            
#         # Create dataframe and pivot for heatmap
#         heatmap_df = pd.DataFrame(heatmap_data)
#         heatmap_pivot = heatmap_df.set_index("domain_id")
        
#         # Create heatmap
#         plt.figure(figsize=(12, len(domain_ids) * 0.4 + 2))
#         sns.heatmap(heatmap_pivot, annot=True, cmap="viridis", fmt=".3f")
#         plt.title("Average RMSF by Domain and Temperature")
#         plt.xlabel("Temperature (K)")
#         plt.ylabel("Domain ID")
#         plt.tight_layout()
        
#         # Save figure
#         output_path = os.path.join(vis_dir, "temperature_summary.png")
#         plt.savefig(output_path, dpi=300)
#         plt.close()
        
#         logging.info(f"Temperature summary heatmap saved to {output_path}")
#         return output_path
#     except Exception as e:
#         logging.error(f"Failed to create temperature summary heatmap: {e}")
#         return None

# def create_temperature_average_summary(temperature_average: pd.DataFrame, 
#                                      output_dir: str) -> Optional[str]:
#     """
#     Create a visualization showing average RMSF across temperatures.
    
#     Args:
#         temperature_average: DataFrame with average RMSF values across all temperatures
#         output_dir: Directory to save visualization
        
#     Returns:
#         Path to the saved figure or None if creation fails
#     """
#     try:
#         # Ensure output directory exists
#         vis_dir = os.path.join(output_dir, "visualizations")
#         os.makedirs(vis_dir, exist_ok=True)
        
#         if temperature_average is None or temperature_average.empty:
#             logging.warning("No temperature average data available for summary")
#             return None
            
#         # Group by domain_id and calculate statistics
#         domain_stats = temperature_average.groupby("domain_id")["rmsf_average"].agg(
#             ["mean", "std", "min", "max"]).reset_index()
        
#         # Sort by mean RMSF
#         domain_stats = domain_stats.sort_values("mean", ascending=False)
        
#         # Create bar plot
#         plt.figure(figsize=(12, 8))
#         plt.bar(domain_stats["domain_id"], domain_stats["mean"], yerr=domain_stats["std"])
#         plt.xticks(rotation=90)
#         plt.title("Average RMSF by Domain (Across All Temperatures)")
#         plt.xlabel("Domain ID")
#         plt.ylabel("Average RMSF (nm)")
#         plt.tight_layout()
        
#         # Save figure
#         output_path = os.path.join(vis_dir, "temperature_average_summary.png")
#         plt.savefig(output_path, dpi=300)
#         plt.close()
        
#         logging.info(f"Temperature average summary saved to {output_path}")
#         return output_path
#     except Exception as e:
#         logging.error(f"Failed to create temperature average summary: {e}")
#         return None

# def create_rmsf_distribution_plots(rmsf_data: Dict[str, pd.DataFrame], 
#                                   output_dir: str) -> Optional[str]:
#     """
#     Create distribution plots (violin plot and histogram) showing RMSF distribution by temperature.
    
#     Args:
#         rmsf_data: Dictionary with RMSF data for all temperatures
#         output_dir: Directory to save visualization
        
#     Returns:
#         Path to the saved figure or None if creation fails
#     """
#     try:
#         # Ensure output directory exists
#         vis_dir = os.path.join(output_dir, "visualizations")
#         os.makedirs(vis_dir, exist_ok=True)
        
#         # Extract temperature values
#         temps = [temp for temp in rmsf_data.keys() if temp != "average"]
        
#         if not temps:
#             logging.warning("No temperature data available for distribution plots")
#             return None
            
#         # Prepare data for plotting
#         dist_data = []
        
#         for temp in temps:
#             if temp in rmsf_data:
#                 temp_df = rmsf_data[temp]
#                 rmsf_col = f"rmsf_{temp}"
                
#                 if rmsf_col in temp_df.columns:
#                     for _, row in temp_df.iterrows():
#                         dist_data.append({
#                             "Temperature": temp,
#                             "RMSF": row[rmsf_col]
#                         })
        
#         if not dist_data:
#             logging.warning("No data available for distribution plots")
#             return None
            
#         # Create dataframe for plotting
#         dist_df = pd.DataFrame(dist_data)
        
#         # Create violin plot
#         plt.figure(figsize=(10, 6))
#         sns.violinplot(x="Temperature", y="RMSF", data=dist_df)
#         plt.title("RMSF Distribution by Temperature")
#         plt.xlabel("Temperature (K)")
#         plt.ylabel("RMSF (nm)")
#         plt.tight_layout()
        
#         # Save violin plot
#         violin_path = os.path.join(vis_dir, "rmsf_violin_plot.png")
#         plt.savefig(violin_path, dpi=300)
#         plt.close()
        
#         # Create histogram
#         plt.figure(figsize=(10, 6))
#         for temp in temps:
#             temp_data = dist_df[dist_df["Temperature"] == temp]["RMSF"]
#             if not temp_data.empty:
#                 sns.histplot(temp_data, kde=True, label=f"{temp}K")
        
#         plt.title("RMSF Histogram by Temperature")
#         plt.xlabel("RMSF (nm)")
#         plt.ylabel("Frequency")
#         plt.legend()
#         plt.tight_layout()
        
#         # Save histogram
#         hist_path = os.path.join(vis_dir, "rmsf_histogram.png")
#         plt.savefig(hist_path, dpi=300)
#         plt.close()
        
#         logging.info(f"RMSF distribution plots saved to {violin_path} and {hist_path}")
#         return violin_path
#     except Exception as e:
#         logging.error(f"Failed to create RMSF distribution plots: {e}")
#         return None

# def create_amino_acid_rmsf_plot(rmsf_data: Dict[str, pd.DataFrame], 
#                               output_dir: str) -> Optional[str]:
#     """
#     Create a violin plot showing RMSF distribution by amino acid type.
    
#     Args:
#         rmsf_data: Dictionary with RMSF data for all temperatures
#         output_dir: Directory to save visualization
        
#     Returns:
#         Path to the saved figure or None if creation fails
#     """
#     try:
#         # Ensure output directory exists
#         vis_dir = os.path.join(output_dir, "visualizations")
#         os.makedirs(vis_dir, exist_ok=True)
        
#         # Use temperature average if available
#         if "average" in rmsf_data and not rmsf_data["average"].empty:
#             aa_data = []
            
#             avg_df = rmsf_data["average"]
#             for _, row in avg_df.iterrows():
#                 aa_data.append({
#                     "Residue": row["resname"],
#                     "RMSF": row["rmsf_average"]
#                 })
                
#             # Create dataframe for plotting
#             aa_df = pd.DataFrame(aa_data)
            
#             # Create violin plot
#             plt.figure(figsize=(14, 8))
#             sns.violinplot(x="Residue", y="RMSF", data=aa_df, order=sorted(aa_df["Residue"].unique()))
#             plt.title("RMSF Distribution by Amino Acid Type")
#             plt.xlabel("Amino Acid")
#             plt.ylabel("RMSF (nm)")
#             plt.xticks(rotation=45)
#             plt.tight_layout()
            
#             # Save figure
#             output_path = os.path.join(vis_dir, "amino_acid_rmsf_violin_plot.png")
#             plt.savefig(output_path, dpi=300)
#             plt.close()
            
#             logging.info(f"Amino acid RMSF violin plot saved to {output_path}")
#             return output_path
#         else:
#             logging.warning("No average temperature data available for amino acid plot")
#             return None
#     except Exception as e:
#         logging.error(f"Failed to create amino acid RMSF plot: {e}")
#         return None

# def create_replica_variance_plot(rmsf_data: Dict[str, Dict[str, pd.DataFrame]],
#                                output_dir: str) -> Optional[str]:
#     """
#     Create a plot showing variance of RMSF values across different replicas.
    
#     Args:
#         rmsf_data: Dictionary with RMSF data for all temperatures and replicas
#         output_dir: Directory to save visualization
        
#     Returns:
#         Path to the saved figure or None if creation fails
#     """
#     try:
#         # Ensure output directory exists
#         vis_dir = os.path.join(output_dir, "visualizations")
#         os.makedirs(vis_dir, exist_ok=True)
        
#         # Extract temperatures
#         temps = list(rmsf_data.keys())
        
#         if not temps:
#             logging.warning("No temperature data available for replica variance plot")
#             return None
            
#         # Calculate variance for each temperature
#         variance_data = []
        
#         for temp in temps:
#             replicas = rmsf_data.get(temp, {})
            
#             if replicas:
#                 # Get all domain_ids and resids
#                 domain_resids = set()
                
#                 for replica, df in replicas.items():
#                     if df is not None and not df.empty:
#                         for _, row in df.iterrows():
#                             domain_resids.add((row["domain_id"], row["resid"]))
                
#                 # Calculate variance for each domain_id and resid
#                 for domain_id, resid in domain_resids:
#                     rmsf_values = []
                    
#                     for replica, df in replicas.items():
#                         if df is not None and not df.empty:
#                             mask = (df["domain_id"] == domain_id) & (df["resid"] == resid)
#                             if mask.any():
#                                 rmsf_values.append(df.loc[mask, f"rmsf_{temp}"].values[0])
                    
#                     if len(rmsf_values) > 1:
#                         variance_data.append({
#                             "Temperature": temp,
#                             "Domain": domain_id,
#                             "Resid": resid,
#                             "Variance": np.var(rmsf_values)
#                         })
        
#         if not variance_data:
#             logging.warning("No data available for replica variance plot")
#             return None
            
#         # Create dataframe for plotting
#         variance_df = pd.DataFrame(variance_data)
        
#         # Create box plot
#         plt.figure(figsize=(10, 6))
#         sns.boxplot(x="Temperature", y="Variance", data=variance_df)
#         plt.title("RMSF Variance Across Replicas")
#         plt.xlabel("Temperature (K)")
#         plt.ylabel("Variance of RMSF (nm²)")
#         plt.tight_layout()
        
#         # Save figure
#         output_path = os.path.join(vis_dir, "replica_variance_plot.png")
#         plt.savefig(output_path, dpi=300)
#         plt.close()
        
#         logging.info(f"Replica variance plot saved to {output_path}")
#         return output_path
#     except Exception as e:
#         logging.error(f"Failed to create replica variance plot: {e}")
#         return None

# def create_dssp_rmsf_correlation_plot(feature_dfs: Dict[str, pd.DataFrame],
#                                     output_dir: str) -> Optional[str]:
#     """
#     Create a visualization showing the relationship between secondary structure and RMSF values.
    
#     Args:
#         feature_dfs: Dictionary with ML feature dataframes
#         output_dir: Directory to save visualization
        
#     Returns:
#         Path to the saved figure or None if creation fails
#     """
#     try:
#         # Ensure output directory exists
#         vis_dir = os.path.join(output_dir, "visualizations")
#         os.makedirs(vis_dir, exist_ok=True)
        
#         # Use average temperature data if available
#         if "average" in feature_dfs and not feature_dfs["average"].empty:
#             avg_df = feature_dfs["average"]
            
#             if "dssp" in avg_df.columns and "rmsf_average" in avg_df.columns:
#                 # Group by DSSP code and calculate statistics
#                 dssp_stats = avg_df.groupby("dssp")["rmsf_average"].agg(
#                     ["mean", "std", "count"]).reset_index()
                
#                 # Sort by count (to prioritize common secondary structures)
#                 dssp_stats = dssp_stats.sort_values("count", ascending=False)
                
#                 # Create bar plot
#                 plt.figure(figsize=(12, 8))
#                 plt.bar(dssp_stats["dssp"], dssp_stats["mean"], yerr=dssp_stats["std"])
                
#                 # Add count as text on each bar
#                 for i, row in dssp_stats.iterrows():
#                     plt.text(i, row["mean"] + row["std"] + 0.01, 
#                             f"n={int(row['count'])}", 
#                             ha='center', va='bottom', rotation=0)
                
#                 plt.title("Average RMSF by Secondary Structure (DSSP)")
#                 plt.xlabel("DSSP Code")
#                 plt.ylabel("Average RMSF (nm)")
#                 plt.tight_layout()
                
#                 # Save figure
#                 output_path = os.path.join(vis_dir, "dssp_rmsf_correlation_plot.png")
#                 plt.savefig(output_path, dpi=300)
#                 plt.close()
                
#                 logging.info(f"DSSP vs RMSF correlation plot saved to {output_path}")
#                 return output_path
#             else:
#                 logging.warning("DSSP or RMSF data not found in feature dataframe")
#                 return None
#         else:
#             logging.warning("No average temperature data available for DSSP correlation plot")
#             return None
#     except Exception as e:
#         logging.error(f"Failed to create DSSP vs RMSF correlation plot: {e}")
#         return None

# def create_feature_correlation_plot(feature_dfs: Dict[str, pd.DataFrame],
#                                   output_dir: str) -> Optional[str]:
#     """
#     Create a visualization highlighting relationships between structural features and RMSF.
    
#     Args:
#         feature_dfs: Dictionary with ML feature dataframes
#         output_dir: Directory to save visualization
        
#     Returns:
#         Path to the saved figure or None if creation fails
#     """
#     try:
#         # Ensure output directory exists
#         vis_dir = os.path.join(output_dir, "visualizations")
#         os.makedirs(vis_dir, exist_ok=True)
        
#         # Use average temperature data if available
#         if "average" in feature_dfs and not feature_dfs["average"].empty:
#             avg_df = feature_dfs["average"]
            
#             # Select numerical columns for correlation
#             numerical_cols = []
#             for col in avg_df.columns:
#                 if col.startswith("rmsf_") or col == "normalized_resid" or col.endswith("_encoded"):
#                     numerical_cols.append(col)
            
#             if not numerical_cols:
#                 logging.warning("No numerical feature columns found for correlation plot")
#                 return None
                
#             # Calculate correlation
#             corr_df = avg_df[numerical_cols].corr()
            
#             # Create heatmap
#             plt.figure(figsize=(10, 8))
#             sns.heatmap(corr_df, annot=True, cmap="coolwarm", fmt=".2f", 
#                        vmin=-1, vmax=1, center=0)
#             plt.title("Correlation Between Features and RMSF")
#             plt.tight_layout()
            
#             # Save figure
#             output_path = os.path.join(vis_dir, "feature_correlation_plot.png")
#             plt.savefig(output_path, dpi=300)
#             plt.close()
            
#             logging.info(f"Feature correlation plot saved to {output_path}")
#             return output_path
#         else:
#             logging.warning("No average temperature data available for feature correlation plot")
#             return None
#     except Exception as e:
#         logging.error(f"Failed to create feature correlation plot: {e}")
#         return None

# def generate_visualizations(rmsf_results: Dict[str, Any], 
#                           ml_results: Dict[str, Any],
#                           domain_results: Dict[str, Dict[str, Any]],
#                           config: Dict[str, Any]) -> Dict[str, Any]:
#     """
#     Generate all required visualizations.
    
#     Args:
#         rmsf_results: Dictionary with RMSF processing results
#         ml_results: Dictionary with ML feature processing results
#         domain_results: Dictionary with processing results for all domains
#         config: Configuration dictionary
        
#     Returns:
#         Dictionary with visualization results
#     """
#     output_dir = config.get("output", {}).get("base_dir", "./outputs")
    
#     # Extract required data
#     replica_averages = rmsf_results.get("replica_averages", {})
#     temperature_average = rmsf_results.get("temperature_average")
#     combined_rmsf_data = rmsf_results.get("combined_rmsf_data", {})
#     feature_dfs = ml_results.get("feature_dfs", {})
    
#     # Generate visualizations
#     results = {}
    
#     # Temperature summary heatmap
#     results["temperature_summary"] = create_temperature_summary_heatmap(
#         replica_averages, output_dir)
        
#     # Temperature average summary
#     results["temperature_average_summary"] = create_temperature_average_summary(
#         temperature_average, output_dir)
        
#     # RMSF distribution plots
#     results["rmsf_distribution"] = create_rmsf_distribution_plots(
#         replica_averages, output_dir)
        
#     # Amino acid RMSF plot
#     results["amino_acid_rmsf"] = create_amino_acid_rmsf_plot(
#         {"average": temperature_average}, output_dir)
        
#     # Replica variance plot
#     results["replica_variance"] = create_replica_variance_plot(
#         combined_rmsf_data, output_dir)
        
#     # DSSP vs RMSF correlation plot
#     results["dssp_rmsf_correlation"] = create_dssp_rmsf_correlation_plot(
#         feature_dfs, output_dir)
        
#     # Feature correlation plot
#     results["feature_correlation"] = create_feature_correlation_plot(
#         feature_dfs, output_dir)
    
#     return results
# EOFMARKER
# echo "Created src/mdcath/processing/visualization.py"

# # ------------------------------------
# # CONFIG FILES
# # ------------------------------------

# # Create src/mdcath/config/default_config.yaml
# cat > src/mdcath/config/default_config.yaml << 'EOFMARKER'
# input:
#   mdcath_folder: "/mnt/datasets/MD_CATH/data"  # Path to the mdCATH folder
#   domain_ids: []   # Empty means process default domain (12asA00)

# temperatures: [320, 348, 379, 413, 450]
# num_replicas: 5  # Number of replicas to process per temperature

# output:
#   base_dir: "./outputs"

# processing:
#   frame_selection:
#     method: "rmsd"  # Options: regular, rmsd, gyration, random
#     num_frames: 1   # Number of frames to extract per domain/temperature
#     cluster_method: "kmeans"  # For RMSD-based selection

#   pdb_cleaning:
#     replace_chain_0_with_A: true
#     fix_atom_numbering: true
#     correct_unusual_residue_names: true
#     add_cryst1_record: true  # Add CRYST1 record for MSMS compatibility
#     remove_hydrogens: false  # Whether to remove hydrogen atoms

#   ml_feature_extraction:
#     min_residues_per_domain: 0
#     max_residues_per_domain: 50000
#     normalize_features: true
#     include_secondary_structure: true
#     include_core_exterior: true
#     include_dssp: true  # Extract and include per-residue DSSP data

#   core_exterior:
#     method: "msms"  # Options: msms, biopython, fallback
#     msms_executable_dir: "./msms"  # Path to MSMS executables
#     ses_threshold: 1.0  # Threshold for classifying residues (Å²)
#     sasa_threshold: 20.0  # Threshold for Biopython SASA (Å²)

#   voxelization:
#     frame_edge_length: 12.0  # Physical size of the voxel grid (Å)
#     voxels_per_side: 21  # Number of voxels along each dimension
#     atom_encoder: "CNOCBCA"  # Atom types to include (options: CNO, CNOCB, CNOCBCA)
#     encode_cb: true  # Whether to include CB atoms
#     compression_gzip: true  # Whether to compress the output files
#     voxelise_all_states: false  # Whether to voxelize all states in NMR structures

# performance:
#   num_cores: 0  # 0 means auto-detect (use max available cores - 2)
#   batch_size: 100
#   memory_limit_gb: 0  # 0 means no limit
#   use_gpu: false  # Whether to use GPU acceleration if available

# logging:
#   verbose: true
#   level: "INFO"
#   console_level: "INFO"
#   file_level: "DEBUG"
#   show_progress_bars: true
# EOFMARKER
# echo "Created src/mdcath/config/default_config.yaml"

# # ------------------------------------
# # MAIN PROJECT FILES
# # ------------------------------------

# # Create main.py
# cat > main.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Main entry point for mdCATH dataset processing.
# """

# import os
# import sys
# import logging
# import argparse
# import yaml
# from typing import Dict, Any, List, Optional

# from src.mdcath.core.data_loader import process_domains, scan_available_domains
# from src.mdcath.processing.rmsf import process_rmsf_data
# from src.mdcath.processing.pdb import process_pdb_data
# from src.mdcath.processing.features import process_ml_features
# from src.mdcath.processing.voxelizer import voxelize_domains
# from src.mdcath.processing.visualization import generate_visualizations

# def setup_logging(config: Dict[str, Any]) -> None:
#     """
#     Set up logging configuration.

#     Args:
#         config: Configuration dictionary
#     """
#     log_config = config.get("logging", {})
#     level = getattr(logging, log_config.get("level", "INFO"))
#     console_level = getattr(logging, log_config.get("console_level", "INFO"))

#     # Configure root logger
#     logging.basicConfig(level=level,
#                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
#                        handlers=[
#                            logging.FileHandler("mdcath_processing.log"),
#                            logging.StreamHandler()
#                        ])

#     # Configure console handler separately
#     console_handler = logging.StreamHandler()
#     console_handler.setLevel(console_level)
#     console_formatter = logging.Formatter("%(levelname)s: %(message)s")
#     console_handler.setFormatter(console_formatter)

#     # Replace the default console handler
#     root_logger = logging.getLogger()
#     for handler in root_logger.handlers:
#         if isinstance(handler, logging.StreamHandler) and handler.stream == sys.stdout:
#             root_logger.removeHandler(handler)

#     root_logger.addHandler(console_handler)

# def load_config(config_path: str) -> Dict[str, Any]:
#     """
#     Load configuration from YAML file.

#     Args:
#         config_path: Path to configuration file

#     Returns:
#         Configuration dictionary
#     """
#     try:
#         with open(config_path, 'r') as f:
#             config = yaml.safe_load(f)
#         return config
#     except Exception as e:
#         logging.error(f"Failed to load configuration: {e}")
#         return {}

# def scan_available_domains(mdcath_folder: str) -> List[str]:
#     """
#     Scan the mdCATH folder for available domains.

#     Args:
#         mdcath_folder: Path to the mdCATH folder

#     Returns:
#         List of domain IDs
#     """
#     domain_ids = []
#     try:
#         if not os.path.exists(mdcath_folder):
#             logging.error(f"mdCATH folder not found: {mdcath_folder}")
#             return []

#         # Find all H5 files in the directory
#         for file in os.listdir(mdcath_folder):
#             if file.endswith(".h5") and file.startswith("mdcath_dataset_"):
#                 # Extract domain ID from file name
#                 domain_id = file.replace("mdcath_dataset_", "").replace(".h5", "")
#                 domain_ids.append(domain_id)

#         return sorted(domain_ids)
#     except Exception as e:
#         logging.error(f"Failed to scan for available domains: {e}")
#         return []

# def process_mdcath(config: Dict[str, Any]) -> int:
#     """
#     Process mdCATH dataset according to configuration.

#     Args:
#         config: Configuration dictionary

#     Returns:
#         Exit code (0 for success, non-zero for failure)
#     """
#     try:
#         # Get input parameters
#         input_config = config.get("input", {})
#         mdcath_folder = input_config.get("mdcath_folder", "/mnt/datasets/MD_CATH/data")
#         domain_ids = input_config.get("domain_ids", [])

#         # If domain_ids is empty, process all available domains
#         if not domain_ids:
#             logging.info("No specific domain IDs provided, scanning for all available domains")
#             domain_ids = scan_available_domains(mdcath_folder)
#             logging.info(f"Found {len(domain_ids)} domains")
            
#             # If still empty, something is wrong with the data directory
#             if not domain_ids:
#                 logging.error(f"No domains found in {mdcath_folder}")
#                 return 1

#         # Process domains
#         logging.info(f"Processing {len(domain_ids)} domains")
#         domain_results = process_domains(domain_ids, mdcath_folder, config,
#                                         num_cores=config.get("performance", {}).get("num_cores", 0))

#         # Process RMSF data
#         logging.info("Processing RMSF data")
#         rmsf_results = process_rmsf_data(domain_results, config)

#         # Process PDB data
#         logging.info("Processing PDB data")
#         pdb_results = process_pdb_data(domain_results, config)

#         # Process ML features
#         logging.info("Generating ML features")
#         ml_results = process_ml_features(rmsf_results, pdb_results, domain_results, config)

#         # Voxelize domains
#         logging.info("Voxelizing domains")
#         voxel_results = voxelize_domains(pdb_results, config)

#         # Generate visualizations
#         logging.info("Generating visualizations")
#         vis_results = generate_visualizations(rmsf_results, ml_results, domain_results, config)

#         logging.info("Processing completed successfully")
#         return 0
#     except Exception as e:
#         logging.error(f"Processing failed: {e}")
#         return 1

# def main() -> int:
#     """
#     Main entry point.

#     Returns:
#         Exit code
#     """
#     parser = argparse.ArgumentParser(description="Process mdCATH dataset for ML applications")
#     parser.add_argument("-c", "--config", default="./src/mdcath/config/default_config.yaml",
#                        help="Path to configuration file")
#     args = parser.parse_args()

#     # Load configuration
#     config = load_config(args.config)
#     if not config:
#         return 1

#     # Set up logging
#     setup_logging(config)

#     # Process mdCATH dataset
#     return process_mdcath(config)

# if __name__ == "__main__":
#     sys.exit(main())
# EOFMARKER
# echo "Created main.py"

# # Create setup.py
# cat > setup.py << 'EOFMARKER'
# #!/usr/bin/env python3
# """
# Setup script for mdCATH
# """

# from setuptools import setup, find_packages

# setup(
#     name="mdcath",
#     version="0.1.0",
#     description="Process mdCATH dataset for ML applications",
#     author="Biochemistry Team",
#     author_email="info@example.com",
#     packages=find_packages("src"),
#     package_dir={"": "src"},
#     python_requires=">=3.7",
#     install_requires=[
#         "h5py",
#         "numpy",
#         "pandas",
#         "biopython",
#         "pyyaml",
#         "matplotlib",
#         "seaborn",
#         "tqdm",
#     ],
#     classifiers=[
#         "Development Status :: 3 - Alpha",
#         "Intended Audience :: Science/Research",
#         "License :: OSI Approved :: MIT License",
#         "Programming Language :: Python :: 3",
#         "Programming Language :: Python :: 3.7",
#         "Programming Language :: Python :: 3.8",
#         "Programming Language :: Python :: 3.9",
#     ],
# )
# EOFMARKER
# echo "Created setup.py"

#!/bin/bash

# Make the script exit on error
set -e

echo "Generating remaining code files for mdCATH project..."

# ------------------------------------
# FINISH requirements.txt
# ------------------------------------

# Since the requirements.txt was cut off in the middle, let's rewrite it completely
cat > requirements.txt << 'EOFMARKER'
h5py>=3.1.0
numpy>=1.19.0
pandas>=1.1.0
biopython>=1.78
pyyaml>=5.4.0
matplotlib>=3.3.0
seaborn>=0.11.0
tqdm>=4.50.0
# Install aposteriori separately via pip
EOFMARKER
echo "Created requirements.txt"

# ------------------------------------
# CREATE README.md
# ------------------------------------

cat > README.md << 'EOFMARKER'
# mdCATH Dataset Processing Project for RMSF Prediction

## Overview

This project provides a comprehensive data processing pipeline for the mdCATH protein dynamics dataset. It extracts, transforms, and organizes the data into formats optimized for training machine learning models that predict Root Mean Square Fluctuation (RMSF) from protein structure information.

The resulting datasets enable the development of ML architectures that can accurately predict protein dynamics from structural features.

## Features

- Extract RMSF data for each simulation and replica
- Process PDB data from H5 files with proper cleaning
- Calculate core/exterior classification for protein residues
- Generate comprehensive ML features including DSSP secondary structure
- Voxelize protein structures using aposteriori
- Create visualizations to analyze and validate the data

## Installation

### Prerequisites

- Python 3.7 or higher
- aposteriori (for voxelization)
- MSMS (optional, for improved core/exterior classification)

### Setup

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/mdcath-processing.git
   cd mdcath-processing
   ```

2. Create a virtual environment:
   ```
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install the required packages:
   ```
   pip install -r requirements.txt
   ```

4. Install aposteriori:
   ```
   pip install aposteriori
   ```

5. Install the package in development mode:
   ```
   pip install -e .
   ```

## Usage

### Basic Usage

Process all domains in the default location:

```
python main.py
```

### Configuration

You can customize the processing by creating your own configuration file:

```
python main.py -c path/to/your/config.yaml
```

See `src/mdcath/config/default_config.yaml` for available configuration options.

### Output Structure

The processed data will be organized in the following directory structure:

```
outputs/
├── RMSF/
│   ├── replicas/
│   │   ├── replica_0/
│   │   │   ├── 320/
│   │   │   ├── ...
│   │   ├── ...
│   └── replica_average/
│       ├── 320/
│       ├── ...
│       └── average/
├── pdbs/
├── frames/
├── voxelized/
├── ML_features/
└── visualizations/
```

## Processing Steps

1. **Data Loading**: Extract data from H5 files
2. **RMSF Processing**: Calculate and average RMSF values
3. **PDB Processing**: Clean and format PDB files
4. **Core/Exterior Classification**: Identify core and exterior residues
5. **ML Feature Generation**: Create features for machine learning
6. **Voxelization**: Create 3D grid representations for CNN models
7. **Visualization**: Generate plots and visualizations

## License

MIT License - See LICENSE file for details.

## Acknowledgments

- mdCATH dataset created by [Institution]
- MSMS software by [Creator]
- aposteriori voxelization tool by [Creator]
EOFMARKER
echo "Created README.md"

# ------------------------------------
# CREATE test_h5_loading.py
# ------------------------------------

cat > test_h5_loading.py << 'EOFMARKER'
#!/usr/bin/env python3
"""
Simple test script for verifying mdCATH processing
"""

import os
import sys
import logging
from src.mdcath.core.data_loader import H5DataLoader

def test_h5_loading():
    """
    Test loading data from an H5 file
    """
    # Set up logging
    logging.basicConfig(level=logging.INFO,
                       format="%(levelname)s: %(message)s")
    
    # Get mdCATH folder from command line or use default
    mdcath_folder = sys.argv[1] if len(sys.argv) > 1 else "/mnt/datasets/MD_CATH/data"
    
    # Test domain (hardcoded for simplicity)
    test_domain = "12asA00"
    
    # Test H5 path
    h5_path = os.path.join(mdcath_folder, f"mdcath_dataset_{test_domain}.h5")
    
    if not os.path.exists(h5_path):
        logging.error(f"Test file not found: {h5_path}")
        return False
    
    try:
        # Initialize H5DataLoader
        config = {"temperatures": [320, 348, 379, 413, 450], "num_replicas": 5}
        loader = H5DataLoader(h5_path, config)
        
        # Test validation
        logging.info("Testing H5 validation...")
        valid = loader._validate_h5()
        logging.info(f"H5 validation result: {valid}")
        
        # Test RMSF extraction
        logging.info("Testing RMSF extraction...")
        rmsf_df = loader.extract_rmsf("320", "0")
        if rmsf_df is not None:
            logging.info(f"Successfully extracted RMSF data: {len(rmsf_df)} rows")
            logging.info(f"First few rows: \n{rmsf_df.head()}")
        else:
            logging.error("Failed to extract RMSF data")
        
        # Test PDB extraction
        logging.info("Testing PDB extraction...")
        pdb_str = loader.extract_pdb()
        if pdb_str:
            pdb_lines = pdb_str.split('\n')
            logging.info(f"Successfully extracted PDB data: {len(pdb_lines)} lines")
            logging.info(f"First few lines: \n{pdb_lines[:5]}")
        else:
            logging.error("Failed to extract PDB data")
        
        # Test DSSP extraction
        logging.info("Testing DSSP extraction...")
        dssp_df = loader.extract_dssp("320", "0")
        if dssp_df is not None:
            logging.info(f"Successfully extracted DSSP data: {len(dssp_df)} rows")
            logging.info(f"First few rows: \n{dssp_df.head()}")
        else:
            logging.error("Failed to extract DSSP data")
        
        # Test coordinate extraction
        logging.info("Testing coordinate extraction...")
        coords_result = loader.extract_coordinates("320", "0")
        if coords_result:
            coords, resids, resnames = coords_result
            logging.info(f"Successfully extracted coordinates: {coords.shape} shape")
            logging.info(f"Number of residues: {len(resids)}")
        else:
            logging.error("Failed to extract coordinates")
        
        return True
    except Exception as e:
        logging.error(f"Test failed with error: {e}")
        return False

if __name__ == "__main__":
    success = test_h5_loading()
    sys.exit(0 if success else 1)
EOFMARKER
echo "Created test_h5_loading.py"

# ------------------------------------
# CREATE LICENSE FILE
# ------------------------------------

cat > LICENSE << 'EOFMARKER'
MIT License

Copyright (c) 2025 mdCATH Processing Project

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
EOFMARKER
echo "Created LICENSE file"

# Make the files executable
chmod +x main.py
chmod +x test_h5_loading.py

echo "Code generation completed successfully!"
echo "Run './setup_structure.sh' to create directories, then './generate_code.sh' for core code, and './generate_code2.sh' for remaining files."
echo "After that, you can install and run the project."