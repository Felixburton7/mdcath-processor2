#!/usr/bin/env python3
"""
Processing module for RMSF data extraction and averaging.
"""

import os
import logging
import numpy as np
import pandas as pd
from typing import List, Dict, Optional, Any, Union
from concurrent.futures import ProcessPoolExecutor

def calculate_replica_averages(rmsf_data: Dict[str, Dict[str, pd.DataFrame]],
                              temperature: str) -> Optional[pd.DataFrame]:
    """
    Calculate average RMSF across all replicas for a specific temperature.

    Args:
        rmsf_data: Dictionary with RMSF data for all replicas
        temperature: Temperature to calculate average for

    Returns:
        DataFrame with average RMSF values or None if calculation fails
    """
    try:
        # Collect all dataframes for this temperature
        dfs = []
        for replica, df in rmsf_data.get(temperature, {}).items():
            if df is not None:
                dfs.append(df)

        if not dfs:
            logging.warning(f"No RMSF data found for temperature {temperature}")
            return None

        # Combine the first dataframe for residue information
        result_df = dfs[0][['domain_id', 'resid', 'resname']].copy()

        # Calculate average RMSF
        rmsf_values = []
        for df in dfs:
            rmsf_col = f"rmsf_{temperature}"
            if rmsf_col in df.columns:
                rmsf_values.append(df[rmsf_col].values)

        if not rmsf_values:
            logging.warning(f"No RMSF values found for temperature {temperature}")
            return None

        # Calculate average
        avg_rmsf = np.mean(rmsf_values, axis=0)
        result_df[f"rmsf_{temperature}"] = avg_rmsf

        return result_df
    except Exception as e:
        logging.error(f"Failed to calculate replica averages for temperature {temperature}: {e}")
        return None

def calculate_temperature_average(replica_averages: Dict[str, pd.DataFrame]) -> Optional[pd.DataFrame]:
    """
    Calculate average RMSF across all temperatures.

    Args:
        replica_averages: Dictionary with replica average RMSF data for all temperatures

    Returns:
        DataFrame with average RMSF values across all temperatures or None if calculation fails
    """
    try:
        if not replica_averages:
            logging.warning("No replica averages found")
            return None

        # Get the first dataframe for base structure
        temps = list(replica_averages.keys())
        first_temp = temps[0]
        result_df = replica_averages[first_temp][['domain_id', 'resid', 'resname']].copy()

        # Collect RMSF values for all temperatures
        rmsf_cols = []
        for temp, df in replica_averages.items():
            rmsf_col = f"rmsf_{temp}"
            if rmsf_col in df.columns:
                result_df[rmsf_col] = df[rmsf_col]
                rmsf_cols.append(rmsf_col)

        if not rmsf_cols:
            logging.warning("No RMSF columns found")
            return None

        # Calculate average across all temperatures
        result_df['rmsf_average'] = result_df[rmsf_cols].mean(axis=1)

        return result_df
    except Exception as e:
        logging.error(f"Failed to calculate temperature average: {e}")
        return None

def save_rmsf_data(rmsf_data: Dict[str, Dict[str, pd.DataFrame]],
                  replica_averages: Dict[str, pd.DataFrame],
                  temperature_average: pd.DataFrame,
                  output_dir: str) -> bool:
    """
    Save RMSF data to CSV files.

    Args:
        rmsf_data: Dictionary with RMSF data for all temperatures and replicas
        replica_averages: Dictionary with replica average RMSF data for all temperatures
        temperature_average: DataFrame with average RMSF values across all temperatures
        output_dir: Directory to save CSV files

    Returns:
        Boolean indicating if saving was successful
    """
    try:
        # Create output directory structure
        os.makedirs(os.path.join(output_dir, "RMSF", "replicas"), exist_ok=True)
        os.makedirs(os.path.join(output_dir, "RMSF", "replica_average", "average"), exist_ok=True)

        # Save replica data
        for temp, replicas in rmsf_data.items():
            for replica, df in replicas.items():
                replica_dir = os.path.join(output_dir, "RMSF", "replicas", f"replica_{replica}", temp)
                os.makedirs(replica_dir, exist_ok=True)

                output_file = os.path.join(replica_dir, f"rmsf_replica{replica}_temperature{temp}.csv")
                df.to_csv(output_file, index=False)
                logging.info(f"Saved RMSF data to {output_file}")

        # Save replica averages
        for temp, df in replica_averages.items():
            temp_dir = os.path.join(output_dir, "RMSF", "replica_average", temp)
            os.makedirs(temp_dir, exist_ok=True)

            output_file = os.path.join(temp_dir, f"rmsf_replica_average_temperature{temp}.csv")
            df.to_csv(output_file, index=False)
            logging.info(f"Saved replica average RMSF data to {output_file}")

        # Save temperature average
        output_file = os.path.join(output_dir, "RMSF", "replica_average", "average",
                                  "rmsf_all_temperatures_all_replicas.csv")
        temperature_average.to_csv(output_file, index=False)
        logging.info(f"Saved temperature average RMSF data to {output_file}")

        return True
    except Exception as e:
        logging.error(f"Failed to save RMSF data: {e}")
        return False

def process_rmsf_data(domain_results: Dict[str, Dict[str, Any]], config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Process RMSF data for all domains.

    Args:
        domain_results: Dictionary with processing results for all domains
        config: Configuration dictionary

    Returns:
        Dictionary with RMSF processing results
    """
    temps = [str(t) for t in config.get("temperatures", [320, 348, 379, 413, 450])]
    output_dir = config.get("output", {}).get("base_dir", "./outputs")

    # Combine RMSF data from all domains
    combined_rmsf_data = {temp: {} for temp in temps}
    for domain_id, result in domain_results.items():
        if not result.get("success", False):
            continue

        rmsf_data = result.get("rmsf_data", {})
        for temp in temps:
            if temp in rmsf_data:
                for replica, df in rmsf_data[temp].items():
                    if replica not in combined_rmsf_data[temp]:
                        combined_rmsf_data[temp][replica] = []
                    combined_rmsf_data[temp][replica].append(df)

    # Concatenate dataframes for each temperature and replica
    for temp in combined_rmsf_data:
        for replica in combined_rmsf_data[temp]:
            if combined_rmsf_data[temp][replica]:
                combined_rmsf_data[temp][replica] = pd.concat(combined_rmsf_data[temp][replica], ignore_index=True)

    # Calculate replica averages
    replica_averages = {}
    for temp in temps:
        avg_df = calculate_replica_averages(combined_rmsf_data, temp)
        if avg_df is not None:
            replica_averages[temp] = avg_df

    # Calculate temperature average
    temperature_average = calculate_temperature_average(replica_averages)

    # Save RMSF data
    save_success = save_rmsf_data(combined_rmsf_data, replica_averages, temperature_average, output_dir)

    return {
        "combined_rmsf_data": combined_rmsf_data,
        "replica_averages": replica_averages,
        "temperature_average": temperature_average,
        "save_success": save_success
    }
