#!/usr/bin/env python3
"""
Module for generating visualizations of processed mdCATH data.
"""

import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, Any, Optional, List, Tuple

def create_temperature_summary_heatmap(rmsf_data: Dict[str, pd.DataFrame], 
                                     output_dir: str) -> Optional[str]:
    """
    Create a heatmap showing RMSF values across temperatures for all domains.
    
    Args:
        rmsf_data: Dictionary with RMSF data for all temperatures
        output_dir: Directory to save visualization
        
    Returns:
        Path to the saved figure or None if creation fails
    """
    try:
        # Ensure output directory exists
        vis_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        
        # Extract temperature values
        temps = [temp for temp in rmsf_data.keys() if temp != "average"]
        
        if not temps:
            logging.warning("No temperature data available for heatmap")
            return None
            
        # Prepare data for heatmap
        domain_ids = set()
        for temp in temps:
            if temp in rmsf_data:
                domain_ids.update(rmsf_data[temp]["domain_id"].unique())
        
        domain_ids = sorted(list(domain_ids))
        
        # Create a dataframe for the heatmap
        heatmap_data = []
        
        for domain_id in domain_ids:
            domain_data = {"domain_id": domain_id}
            
            for temp in temps:
                if temp in rmsf_data:
                    domain_temp_data = rmsf_data[temp][rmsf_data[temp]["domain_id"] == domain_id]
                    if not domain_temp_data.empty:
                        domain_data[f"rmsf_{temp}"] = domain_temp_data[f"rmsf_{temp}"].mean()
            
            heatmap_data.append(domain_data)
        
        if not heatmap_data:
            logging.warning("No data available for heatmap")
            return None
            
        # Create dataframe and pivot for heatmap
        heatmap_df = pd.DataFrame(heatmap_data)
        heatmap_pivot = heatmap_df.set_index("domain_id")
        
        # Create heatmap
        plt.figure(figsize=(12, len(domain_ids) * 0.4 + 2))
        sns.heatmap(heatmap_pivot, annot=True, cmap="viridis", fmt=".3f")
        plt.title("Average RMSF by Domain and Temperature")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Domain ID")
        plt.tight_layout()
        
        # Save figure
        output_path = os.path.join(vis_dir, "temperature_summary.png")
        plt.savefig(output_path, dpi=300)
        plt.close()
        
        logging.info(f"Temperature summary heatmap saved to {output_path}")
        return output_path
    except Exception as e:
        logging.error(f"Failed to create temperature summary heatmap: {e}")
        return None

def create_temperature_average_summary(temperature_average: pd.DataFrame, 
                                     output_dir: str) -> Optional[str]:
    """
    Create a visualization showing average RMSF across temperatures.
    
    Args:
        temperature_average: DataFrame with average RMSF values across all temperatures
        output_dir: Directory to save visualization
        
    Returns:
        Path to the saved figure or None if creation fails
    """
    try:
        # Ensure output directory exists
        vis_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        
        if temperature_average is None or temperature_average.empty:
            logging.warning("No temperature average data available for summary")
            return None
            
        # Group by domain_id and calculate statistics
        domain_stats = temperature_average.groupby("domain_id")["rmsf_average"].agg(
            ["mean", "std", "min", "max"]).reset_index()
        
        # Sort by mean RMSF
        domain_stats = domain_stats.sort_values("mean", ascending=False)
        
        # Create bar plot
        plt.figure(figsize=(12, 8))
        plt.bar(domain_stats["domain_id"], domain_stats["mean"], yerr=domain_stats["std"])
        plt.xticks(rotation=90)
        plt.title("Average RMSF by Domain (Across All Temperatures)")
        plt.xlabel("Domain ID")
        plt.ylabel("Average RMSF (nm)")
        plt.tight_layout()
        
        # Save figure
        output_path = os.path.join(vis_dir, "temperature_average_summary.png")
        plt.savefig(output_path, dpi=300)
        plt.close()
        
        logging.info(f"Temperature average summary saved to {output_path}")
        return output_path
    except Exception as e:
        logging.error(f"Failed to create temperature average summary: {e}")
        return None

def create_rmsf_distribution_plots(rmsf_data: Dict[str, pd.DataFrame], 
                                  output_dir: str) -> Optional[str]:
    """
    Create distribution plots (violin plot and histogram) showing RMSF distribution by temperature.
    
    Args:
        rmsf_data: Dictionary with RMSF data for all temperatures
        output_dir: Directory to save visualization
        
    Returns:
        Path to the saved figure or None if creation fails
    """
    try:
        # Ensure output directory exists
        vis_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        
        # Extract temperature values
        temps = [temp for temp in rmsf_data.keys() if temp != "average"]
        
        if not temps:
            logging.warning("No temperature data available for distribution plots")
            return None
            
        # Prepare data for plotting
        dist_data = []
        
        for temp in temps:
            if temp in rmsf_data:
                temp_df = rmsf_data[temp]
                rmsf_col = f"rmsf_{temp}"
                
                if rmsf_col in temp_df.columns:
                    for _, row in temp_df.iterrows():
                        dist_data.append({
                            "Temperature": temp,
                            "RMSF": row[rmsf_col]
                        })
        
        if not dist_data:
            logging.warning("No data available for distribution plots")
            return None
            
        # Create dataframe for plotting
        dist_df = pd.DataFrame(dist_data)
        
        # Create violin plot
        plt.figure(figsize=(10, 6))
        sns.violinplot(x="Temperature", y="RMSF", data=dist_df)
        plt.title("RMSF Distribution by Temperature")
        plt.xlabel("Temperature (K)")
        plt.ylabel("RMSF (nm)")
        plt.tight_layout()
        
        # Save violin plot
        violin_path = os.path.join(vis_dir, "rmsf_violin_plot.png")
        plt.savefig(violin_path, dpi=300)
        plt.close()
        
        # Create histogram
        plt.figure(figsize=(10, 6))
        for temp in temps:
            temp_data = dist_df[dist_df["Temperature"] == temp]["RMSF"]
            if not temp_data.empty:
                sns.histplot(temp_data, kde=True, label=f"{temp}K")
        
        plt.title("RMSF Histogram by Temperature")
        plt.xlabel("RMSF (nm)")
        plt.ylabel("Frequency")
        plt.legend()
        plt.tight_layout()
        
        # Save histogram
        hist_path = os.path.join(vis_dir, "rmsf_histogram.png")
        plt.savefig(hist_path, dpi=300)
        plt.close()
        
        logging.info(f"RMSF distribution plots saved to {violin_path} and {hist_path}")
        return violin_path
    except Exception as e:
        logging.error(f"Failed to create RMSF distribution plots: {e}")
        return None

def create_amino_acid_rmsf_plot(rmsf_data: Dict[str, pd.DataFrame], 
                              output_dir: str) -> Optional[str]:
    """
    Create a violin plot showing RMSF distribution by amino acid type.
    
    Args:
        rmsf_data: Dictionary with RMSF data for all temperatures
        output_dir: Directory to save visualization
        
    Returns:
        Path to the saved figure or None if creation fails
    """
    try:
        # Ensure output directory exists
        vis_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        
        # Use temperature average if available
        if "average" in rmsf_data and not rmsf_data["average"].empty:
            aa_data = []
            
            avg_df = rmsf_data["average"]
            for _, row in avg_df.iterrows():
                aa_data.append({
                    "Residue": row["resname"],
                    "RMSF": row["rmsf_average"]
                })
                
            # Create dataframe for plotting
            aa_df = pd.DataFrame(aa_data)
            
            # Create violin plot
            plt.figure(figsize=(14, 8))
            sns.violinplot(x="Residue", y="RMSF", data=aa_df, order=sorted(aa_df["Residue"].unique()))
            plt.title("RMSF Distribution by Amino Acid Type")
            plt.xlabel("Amino Acid")
            plt.ylabel("RMSF (nm)")
            plt.xticks(rotation=45)
            plt.tight_layout()
            
            # Save figure
            output_path = os.path.join(vis_dir, "amino_acid_rmsf_violin_plot.png")
            plt.savefig(output_path, dpi=300)
            plt.close()
            
            logging.info(f"Amino acid RMSF violin plot saved to {output_path}")
            return output_path
        else:
            logging.warning("No average temperature data available for amino acid plot")
            return None
    except Exception as e:
        logging.error(f"Failed to create amino acid RMSF plot: {e}")
        return None

def create_replica_variance_plot(rmsf_data: Dict[str, Dict[str, pd.DataFrame]],
                               output_dir: str) -> Optional[str]:
    """
    Create a plot showing variance of RMSF values across different replicas.
    
    Args:
        rmsf_data: Dictionary with RMSF data for all temperatures and replicas
        output_dir: Directory to save visualization
        
    Returns:
        Path to the saved figure or None if creation fails
    """
    try:
        # Ensure output directory exists
        vis_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        
        # Extract temperatures
        temps = list(rmsf_data.keys())
        
        if not temps:
            logging.warning("No temperature data available for replica variance plot")
            return None
            
        # Calculate variance for each temperature
        variance_data = []
        
        for temp in temps:
            replicas = rmsf_data.get(temp, {})
            
            if replicas:
                # Get all domain_ids and resids
                domain_resids = set()
                
                for replica, df in replicas.items():
                    if df is not None and not df.empty:
                        for _, row in df.iterrows():
                            domain_resids.add((row["domain_id"], row["resid"]))
                
                # Calculate variance for each domain_id and resid
                for domain_id, resid in domain_resids:
                    rmsf_values = []
                    
                    for replica, df in replicas.items():
                        if df is not None and not df.empty:
                            mask = (df["domain_id"] == domain_id) & (df["resid"] == resid)
                            if mask.any():
                                rmsf_values.append(df.loc[mask, f"rmsf_{temp}"].values[0])
                    
                    if len(rmsf_values) > 1:
                        variance_data.append({
                            "Temperature": temp,
                            "Domain": domain_id,
                            "Resid": resid,
                            "Variance": np.var(rmsf_values)
                        })
        
        if not variance_data:
            logging.warning("No data available for replica variance plot")
            return None
            
        # Create dataframe for plotting
        variance_df = pd.DataFrame(variance_data)
        
        # Create box plot
        plt.figure(figsize=(10, 6))
        sns.boxplot(x="Temperature", y="Variance", data=variance_df)
        plt.title("RMSF Variance Across Replicas")
        plt.xlabel("Temperature (K)")
        plt.ylabel("Variance of RMSF (nm²)")
        plt.tight_layout()
        
        # Save figure
        output_path = os.path.join(vis_dir, "replica_variance_plot.png")
        plt.savefig(output_path, dpi=300)
        plt.close()
        
        logging.info(f"Replica variance plot saved to {output_path}")
        return output_path
    except Exception as e:
        logging.error(f"Failed to create replica variance plot: {e}")
        return None

def create_dssp_rmsf_correlation_plot(feature_dfs: Dict[str, pd.DataFrame],
                                    output_dir: str) -> Optional[str]:
    """
    Create a visualization showing the relationship between secondary structure and RMSF values.
    
    Args:
        feature_dfs: Dictionary with ML feature dataframes
        output_dir: Directory to save visualization
        
    Returns:
        Path to the saved figure or None if creation fails
    """
    try:
        # Ensure output directory exists
        vis_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        
        # Use average temperature data if available
        if "average" in feature_dfs and not feature_dfs["average"].empty:
            avg_df = feature_dfs["average"]
            
            if "dssp" in avg_df.columns and "rmsf_average" in avg_df.columns:
                # Group by DSSP code and calculate statistics
                dssp_stats = avg_df.groupby("dssp")["rmsf_average"].agg(
                    ["mean", "std", "count"]).reset_index()
                
                # Sort by count (to prioritize common secondary structures)
                dssp_stats = dssp_stats.sort_values("count", ascending=False)
                
                # Create bar plot
                plt.figure(figsize=(12, 8))
                plt.bar(dssp_stats["dssp"], dssp_stats["mean"], yerr=dssp_stats["std"])
                
                # Add count as text on each bar
                for i, row in dssp_stats.iterrows():
                    plt.text(i, row["mean"] + row["std"] + 0.01, 
                            f"n={int(row['count'])}", 
                            ha='center', va='bottom', rotation=0)
                
                plt.title("Average RMSF by Secondary Structure (DSSP)")
                plt.xlabel("DSSP Code")
                plt.ylabel("Average RMSF (nm)")
                plt.tight_layout()
                
                # Save figure
                output_path = os.path.join(vis_dir, "dssp_rmsf_correlation_plot.png")
                plt.savefig(output_path, dpi=300)
                plt.close()
                
                logging.info(f"DSSP vs RMSF correlation plot saved to {output_path}")
                return output_path
            else:
                logging.warning("DSSP or RMSF data not found in feature dataframe")
                return None
        else:
            logging.warning("No average temperature data available for DSSP correlation plot")
            return None
    except Exception as e:
        logging.error(f"Failed to create DSSP vs RMSF correlation plot: {e}")
        return None

def create_feature_correlation_plot(feature_dfs: Dict[str, pd.DataFrame],
                                  output_dir: str) -> Optional[str]:
    """
    Create a visualization highlighting relationships between structural features and RMSF.
    
    Args:
        feature_dfs: Dictionary with ML feature dataframes
        output_dir: Directory to save visualization
        
    Returns:
        Path to the saved figure or None if creation fails
    """
    try:
        # Ensure output directory exists
        vis_dir = os.path.join(output_dir, "visualizations")
        os.makedirs(vis_dir, exist_ok=True)
        
        # Use average temperature data if available
        if "average" in feature_dfs and not feature_dfs["average"].empty:
            avg_df = feature_dfs["average"]
            
            # Select numerical columns for correlation
            numerical_cols = []
            for col in avg_df.columns:
                if col.startswith("rmsf_") or col == "normalized_resid" or col.endswith("_encoded"):
                    numerical_cols.append(col)
            
            if not numerical_cols:
                logging.warning("No numerical feature columns found for correlation plot")
                return None
                
            # Calculate correlation
            corr_df = avg_df[numerical_cols].corr()
            
            # Create heatmap
            plt.figure(figsize=(10, 8))
            sns.heatmap(corr_df, annot=True, cmap="coolwarm", fmt=".2f", 
                       vmin=-1, vmax=1, center=0)
            plt.title("Correlation Between Features and RMSF")
            plt.tight_layout()
            
            # Save figure
            output_path = os.path.join(vis_dir, "feature_correlation_plot.png")
            plt.savefig(output_path, dpi=300)
            plt.close()
            
            logging.info(f"Feature correlation plot saved to {output_path}")
            return output_path
        else:
            logging.warning("No average temperature data available for feature correlation plot")
            return None
    except Exception as e:
        logging.error(f"Failed to create feature correlation plot: {e}")
        return None

def generate_visualizations(rmsf_results: Dict[str, Any], 
                          ml_results: Dict[str, Any],
                          domain_results: Dict[str, Dict[str, Any]],
                          config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate all required visualizations.
    
    Args:
        rmsf_results: Dictionary with RMSF processing results
        ml_results: Dictionary with ML feature processing results
        domain_results: Dictionary with processing results for all domains
        config: Configuration dictionary
        
    Returns:
        Dictionary with visualization results
    """
    output_dir = config.get("output", {}).get("base_dir", "./outputs")
    
    # Extract required data
    replica_averages = rmsf_results.get("replica_averages", {})
    temperature_average = rmsf_results.get("temperature_average")
    combined_rmsf_data = rmsf_results.get("combined_rmsf_data", {})
    feature_dfs = ml_results.get("feature_dfs", {})
    
    # Generate visualizations
    results = {}
    
    # Temperature summary heatmap
    results["temperature_summary"] = create_temperature_summary_heatmap(
        replica_averages, output_dir)
        
    # Temperature average summary
    results["temperature_average_summary"] = create_temperature_average_summary(
        temperature_average, output_dir)
        
    # RMSF distribution plots
    results["rmsf_distribution"] = create_rmsf_distribution_plots(
        replica_averages, output_dir)
        
    # Amino acid RMSF plot
    results["amino_acid_rmsf"] = create_amino_acid_rmsf_plot(
        {"average": temperature_average}, output_dir)
        
    # Replica variance plot
    results["replica_variance"] = create_replica_variance_plot(
        combined_rmsf_data, output_dir)
        
    # DSSP vs RMSF correlation plot
    results["dssp_rmsf_correlation"] = create_dssp_rmsf_correlation_plot(
        feature_dfs, output_dir)
        
    # Feature correlation plot
    results["feature_correlation"] = create_feature_correlation_plot(
        feature_dfs, output_dir)
    
    return results
