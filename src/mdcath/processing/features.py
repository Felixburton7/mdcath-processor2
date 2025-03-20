#!/usr/bin/env python3
"""
Processing module for generating ML features.
"""

import os
import logging
import numpy as np
import pandas as pd
from typing import Dict, Any, Optional, List, Tuple, Union
from tqdm import tqdm
from src.mdcath.processing.core_exterior import compute_core_exterior

def generate_ml_features(rmsf_data: Dict[str, pd.DataFrame],
                       core_exterior_data: Dict[str, pd.DataFrame],
                       dssp_data: Dict[str, Dict[str, pd.DataFrame]],
                       config: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    """
    Generate ML features for all domains.
    """
    try:
        # Get list of all domains
        domain_ids = set()
        for temp, df in rmsf_data.items():
            domain_ids.update(df["domain_id"].unique())

        domain_ids = list(domain_ids)
        logging.info(f"Generating ML features for {len(domain_ids)} domains")

        # Create feature dataframes for each temperature
        temps = [t for t in rmsf_data.keys() if t != "average"]
        feature_dfs = {}

        for temp in temps:
            # Start with RMSF data
            if temp not in rmsf_data:
                logging.warning(f"RMSF data not found for temperature {temp}")
                continue

            df = rmsf_data[temp].copy()

            # Ensure RMSF column is numeric
            rmsf_col = f"rmsf_{temp}"
            if rmsf_col in df.columns:
                # Convert to numeric, coerce errors to NaN
                df[rmsf_col] = pd.to_numeric(df[rmsf_col], errors='coerce')
                # Fill NaN values with 0
                df[rmsf_col] = df[rmsf_col].fillna(0.0)

            # Add protein size
            df["protein_size"] = df.groupby("domain_id")["resid"].transform("count")

            # Add normalized residue position - handle potential division by zero
            df["normalized_resid"] = df.groupby("domain_id")["resid"].transform(
                lambda x: (x - x.min()) / max(x.max() - x.min(), 1)
            )

            # Add core/exterior classification
            for domain_id in df["domain_id"].unique():
                if domain_id in core_exterior_data:
                    core_ext_df = core_exterior_data[domain_id]

                    # Merge core/exterior data
                    domain_mask = df["domain_id"] == domain_id
                    df_domain = df[domain_mask].copy()

                    # Reset index for proper merging
                    df_domain = df_domain.reset_index(drop=True)
                    df_domain = pd.merge(df_domain, core_ext_df, on="resid", how="left")

                    # Update the main dataframe
                    df.loc[domain_mask] = df_domain

            # Fill missing core/exterior values with 'unknown'
            if "core_exterior" in df.columns:
                df["core_exterior"] = df["core_exterior"].fillna("unknown")
            else:
                df["core_exterior"] = "unknown"

            # Add DSSP data
            if temp in dssp_data and "0" in dssp_data[temp]:
                for domain_id in df["domain_id"].unique():
                    domain_dssp = None

                    # Find DSSP data for this domain
                    for replica, dssp_df in dssp_data[temp].items():
                        domain_dssp_subset = dssp_df[dssp_df["domain_id"] == domain_id]
                        if not domain_dssp_subset.empty:
                            domain_dssp = domain_dssp_subset[["resid", "dssp"]]
                            break

                    if domain_dssp is not None:
                        # Ensure resid is numeric in both dataframes
                        domain_dssp["resid"] = pd.to_numeric(domain_dssp["resid"], errors='coerce')
                        df_domain = df[df["domain_id"] == domain_id].copy()
                        df_domain["resid"] = pd.to_numeric(df_domain["resid"], errors='coerce')
                        
                        # Merge DSSP data
                        df_domain = df_domain.reset_index(drop=True)
                        df_domain = pd.merge(df_domain, domain_dssp, on="resid", how="left")

                        # Update the main dataframe
                        df.loc[df["domain_id"] == domain_id] = df_domain

            # Fill missing DSSP values with 'C' (coil)
            if "dssp" in df.columns:
                df["dssp"] = df["dssp"].fillna("C")
            else:
                df["dssp"] = "C"

            # Encode categorical variables
            # Resname encoding
            unique_resnames = df["resname"].unique()
            resname_mapping = {name: i for i, name in enumerate(sorted(unique_resnames))}
            df["resname_encoded"] = df["resname"].map(resname_mapping)

            # Core/exterior encoding
            core_ext_mapping = {"core": 1, "exterior": 2, "unknown": 0}
            df["core_exterior_encoded"] = df["core_exterior"].map(core_ext_mapping)

            # Secondary structure encoding
            # Simplified 3-state encoding: Helix (H,G,I), Sheet (E,B), Coil (others)
            def encode_ss(ss):
                if ss in ["H", "G", "I"]:
                    return 0  # Helix
                elif ss in ["E", "B"]:
                    return 1  # Sheet
                else:
                    return 2  # Coil

            df["secondary_structure_encoded"] = df["dssp"].apply(encode_ss)

            # Add dummy relative accessibility (placeholder for more advanced calculation)
            df["relative_accessibility"] = 1.0

            # Reorder columns to put domain_id first
            cols = df.columns.tolist()
            cols.remove("domain_id")
            cols = ["domain_id"] + cols
            df = df[cols]

            # Store the feature dataframe
            feature_dfs[temp] = df

        # Calculate average features
        if temps:
            # Start with a copy of the first temperature's features
            avg_df = feature_dfs[temps[0]].copy()

            # Collect RMSF columns
            rmsf_cols = [f"rmsf_{temp}" for temp in temps]
            rmsf_vals = []

            for temp in temps:
                if temp in feature_dfs and f"rmsf_{temp}" in feature_dfs[temp].columns:
                    temp_df = feature_dfs[temp]
                    # Extract the domain, resid, and RMSF columns
                    temp_subset = temp_df[["domain_id", "resid", f"rmsf_{temp}"]].copy()
                    
                    # Ensure all RMSF values are numeric
                    temp_subset[f"rmsf_{temp}"] = pd.to_numeric(temp_subset[f"rmsf_{temp}"], errors='coerce').fillna(0.0)
                    
                    rmsf_vals.append(temp_subset)

            # Merge all RMSF values
            if rmsf_vals:
                merged_df = rmsf_vals[0]
                for i in range(1, len(rmsf_vals)):
                    merged_df = pd.merge(merged_df, rmsf_vals[i], on=["domain_id", "resid"])

                # Calculate average RMSF
                merged_df["rmsf_average"] = merged_df[[f"rmsf_{temp}" for temp in temps]].mean(axis=1)

                # Replace the RMSF columns in the average dataframe
                avg_df = pd.merge(avg_df.drop(rmsf_cols, axis=1, errors="ignore"),
                                merged_df[["domain_id", "resid", "rmsf_average"]],
                                on=["domain_id", "resid"])

            feature_dfs["average"] = avg_df

        return feature_dfs
    except Exception as e:
        logging.error(f"Failed to generate ML features: {e}")
        # Print the traceback for better debugging
        import traceback
        logging.error(traceback.format_exc())
        return {}

def save_ml_features(feature_dfs: Dict[str, pd.DataFrame], output_dir: str) -> bool:
    """
    Save ML features to CSV files.

    Args:
        feature_dfs: Dictionary with ML feature dataframes
        output_dir: Directory to save CSV files

    Returns:
        Boolean indicating if saving was successful
    """
    try:
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)

        # Save each feature dataframe
        for temp, df in feature_dfs.items():
            if temp == "average":
                output_file = os.path.join(output_dir, "final_dataset_temperature_average.csv")
            else:
                output_file = os.path.join(output_dir, f"final_dataset_temperature_{temp}.csv")

            df.to_csv(output_file, index=False)
            logging.info(f"Saved ML features to {output_file}")

        return True
    except Exception as e:
        logging.error(f"Failed to save ML features: {e}")
        return False
def process_ml_features(rmsf_results: Dict[str, Any],
                       pdb_results: Dict[str, Any],
                       domain_results: Dict[str, Dict[str, Any]],
                       config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process ML features for all domains.
    
    Args:
        rmsf_results: Dictionary with RMSF processing results
        pdb_results: Dictionary with PDB processing results
        domain_results: Dictionary with processing results for all domains
        config: Configuration dictionary
    
    Returns:
        Dictionary with ML feature processing results
    """
    
    
    output_dir = config.get("output", {}).get("base_dir", "./outputs")
    
    # Extract RMSF data
    replica_averages = rmsf_results.get("replica_averages", {})
    temperature_average = rmsf_results.get("temperature_average")
    
    if not replica_averages:
        logging.error("No RMSF data available for ML feature generation")
        return {"success": False, "error": "No RMSF data available"}
    
    # Create dictionary with all RMSF data
    rmsf_data = replica_averages.copy()
    if temperature_average is not None:
        rmsf_data["average"] = temperature_average
    
    # Extract core/exterior data
    core_exterior_data = {}
    
    # Compute core/exterior classification for each domain
    logging.info("Computing core/exterior classification for domains")
    for domain_id, result in tqdm(pdb_results.items(), desc="Core/exterior classification"):
        if not result.get("pdb_saved", False):
            continue
        
        pdb_path = result.get("pdb_path")
        if not pdb_path or not os.path.exists(pdb_path):
            logging.warning(f"PDB file not found for domain {domain_id}")
            continue
        
        core_ext_df = compute_core_exterior(pdb_path, config)
        if core_ext_df is not None:
            core_exterior_data[domain_id] = core_ext_df
    
    # Collect DSSP data
    dssp_data = {}
    temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
    
    logging.info("Collecting DSSP data")
    for domain_id, result in tqdm(domain_results.items(), desc="Processing DSSP data"):
        if not result.get("success", False):
            continue
        
        # Extract DSSP data
        domain_dssp = result.get("dssp_data", {})
        for temp in temps:
            if temp in domain_dssp:
                if temp not in dssp_data:
                    dssp_data[temp] = {}
                
                for replica, df in domain_dssp[temp].items():
                    if replica not in dssp_data[temp]:
                        dssp_data[temp][replica] = []
                    
                    dssp_data[temp][replica].append(df)
    
    # Concatenate DSSP dataframes
    logging.info("Concatenating DSSP data")
    for temp in temps:
        if temp in dssp_data:
            for replica in dssp_data[temp]:
                if dssp_data[temp][replica]:
                    dssp_data[temp][replica] = pd.concat(dssp_data[temp][replica], ignore_index=True)
    
    # Generate and save ML features
    logging.info("Generating ML features")
    feature_dfs = generate_ml_features(rmsf_data, core_exterior_data, dssp_data, config)
    
    if not feature_dfs:
        logging.error("Failed to generate ML features")
        return {"success": False, "error": "Feature generation failed"}
    
    # Save ML features
    ml_dir = os.path.join(output_dir, "ML_features")
    save_success = save_ml_features(feature_dfs, ml_dir)
    
    return {
        "success": save_success,
        "feature_dfs": feature_dfs,
        "output_dir": ml_dir
    }