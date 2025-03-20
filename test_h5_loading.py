#!/usr/bin/env python3
"""
inspect_hdf5.py

Enhanced script to inspect and compare the contents of mdCATH .h5 files.
"""

import os
import sys
import h5py
import argparse
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Any, Set, Tuple
import matplotlib.pyplot as plt
from io import StringIO
from tabulate import tabulate  # For prettier table outputs


def print_section(title: str, width: int = 80):
    """Print a section header."""
    print("\n" + "=" * width)
    print(title.center(width))
    print("=" * width)


def print_subsection(title: str, width: int = 80):
    """Print a subsection header."""
    print("\n" + "-" * width)
    print(title)
    print("-" * width)


def print_group_attributes(group: h5py.Group, indent: int = 0):
    """Print attributes of an HDF5 group."""
    if len(group.attrs) == 0:
        print(" " * indent + "No attributes")
        return
    
    for attr_name, attr_value in group.attrs.items():
        print(" " * indent + f"{attr_name}: {attr_value}")


def explore_group(group: h5py.Group, max_depth: int = -1, indent: int = 0, current_depth: int = 0, 
                 path: str = "", details: bool = False):
    """
    Recursively explore and print the structure of an HDF5 group.
    
    Args:
        group: The HDF5 group to explore.
        max_depth: Maximum depth to explore. -1 means no limit.
        indent: Initial indentation level.
        current_depth: Current depth in the exploration.
        path: Current path in the hierarchy.
        details: Whether to show detailed information about datasets.
    """
    if max_depth != -1 and current_depth > max_depth:
        print(" " * indent + "...")
        return
    
    for key in sorted(group.keys()):
        item = group[key]
        current_path = f"{path}/{key}" if path else key
        
        if isinstance(item, h5py.Group):
            print(" " * indent + f"{key}/ (Group)")
            if len(item.attrs) > 0:
                print(" " * (indent + 2) + "Attributes:")
                print_group_attributes(item, indent + 4)
            explore_group(item, max_depth, indent + 2, current_depth + 1, current_path, details)
        elif isinstance(item, h5py.Dataset):
            print(" " * indent + f"{key} (Dataset): shape={item.shape}, dtype={item.dtype}")
            if len(item.attrs) > 0:
                print(" " * (indent + 2) + "Attributes:")
                print_group_attributes(item, indent + 4)
            
            if details:
                try:
                    print(" " * (indent + 2) + "Data Sample:")
                    print_dataset_sample(item, indent + 4)
                    
                    # Show statistics for numerical datasets
                    if np.issubdtype(item.dtype, np.number) and item.size > 0:
                        print(" " * (indent + 2) + "Statistics:")
                        print_dataset_statistics(item, indent + 4)
                except Exception as e:
                    print(" " * (indent + 2) + f"Error analyzing dataset: {e}")


def print_dataset_sample(dataset: h5py.Dataset, indent: int = 0, max_elements: int = 5):
    """Print a sample of a dataset."""
    if dataset.shape == ():
        print(" " * indent + f"Value: {dataset[()]}")
        return
    
    if dataset.ndim == 1:
        if len(dataset) <= max_elements:
            sample = dataset[:]
        else:
            sample = dataset[:max_elements]
        
        # Handle string decoding if needed
        if dataset.dtype.kind == 'S':
            try:
                sample = [s.decode('utf-8') if isinstance(s, bytes) else s for s in sample]
            except:
                pass  # Keep as bytes if decoding fails
        
        print(" " * indent + f"Values: {sample}" + 
              (f"... ({len(dataset)} total elements)" if len(dataset) > max_elements else ""))
    else:
        # For multi-dimensional arrays, show first slice
        sample = dataset[(0,) * (dataset.ndim - 1)][:max_elements]
        print(" " * indent + f"First slice: {sample}... (shape: {dataset.shape})")


def print_dataset_statistics(dataset: h5py.Dataset, indent: int = 0):
    """Print statistics for a numerical dataset."""
    try:
        # For large datasets, take a sample
        if dataset.size > 10000:
            # Take a random sample of 10,000 elements
            indices = tuple(np.random.choice(s, min(s, 100)) for s in dataset.shape)
            sample = dataset[indices]
            stats_source = "random sample"
        else:
            sample = dataset[:]
            stats_source = "entire dataset"
        
        # Calculate statistics
        print(" " * indent + f"Statistics ({stats_source}):")
        print(" " * indent + f"  Min: {np.min(sample)}")
        print(" " * indent + f"  Max: {np.max(sample)}")
        print(" " * indent + f"  Mean: {np.mean(sample)}")
        print(" " * indent + f"  Median: {np.median(sample)}")
        print(" " * indent + f"  Std Dev: {np.std(sample)}")
    except Exception as e:
        print(" " * indent + f"  Error calculating statistics: {e}")


def print_domain_info(h5_file: h5py.File, domain_id: str, detailed: bool = False):
    """Print detailed information about a specific domain."""
    if domain_id not in h5_file:
        print(f"Domain {domain_id} not found in the file.")
        return
    
    domain_group = h5_file[domain_id]
    print_section(f"Domain: {domain_id}")
    
    # Print domain attributes
    print_subsection("Domain Attributes")
    print_group_attributes(domain_group)
    
    # Print domain datasets
    print_subsection("Domain-level Datasets")
    for key in sorted(domain_group.keys()):
        if isinstance(domain_group[key], h5py.Dataset):
            dataset = domain_group[key]
            print(f"{key} (Dataset): shape={dataset.shape}, dtype={dataset.dtype}")
            
            if detailed:
                try:
                    print("  Sample Data:")
                    print_dataset_sample(dataset, indent=4)
                    
                    # Show statistics for numerical datasets
                    if np.issubdtype(dataset.dtype, np.number) and dataset.size > 0:
                        print("  Statistics:")
                        print_dataset_statistics(dataset, indent=4)
                except Exception as e:
                    print(f"  Error analyzing dataset: {e}")
    
    # Print temperature groups
    print_subsection("Temperature Groups")
    temp_groups = []
    for key in domain_group.keys():
        try:
            temp = int(key)
            temp_groups.append(key)
        except ValueError:
            continue
    
    if not temp_groups:
        print("No temperature groups found.")
        return
    
    temp_groups = sorted(temp_groups)
    print(f"Available temperatures: {', '.join(temp_groups)}")
    
    # Analyze and show structure for each temperature
    for temp in temp_groups[:2]:  # Limit to first two temperatures to avoid excessive output
        print_subsection(f"Temperature: {temp}")
        temp_group = domain_group[temp]
        
        # Print replica groups
        replica_groups = sorted(list(temp_group.keys()))
        if not replica_groups:
            print("No replica groups found.")
            continue
        
        print(f"Available replicas: {', '.join(replica_groups)}")
        
        # Print information for each replica
        for replica in replica_groups[:2]:  # Limit to first two replicas
            print(f"\nReplica: {replica}")
            replica_group = temp_group[replica]
            
            # Print replica attributes
            if len(replica_group.attrs) > 0:
                print("  Attributes:")
                print_group_attributes(replica_group, 4)
            
            # Print replica datasets
            print("\n  Datasets:")
            for key in sorted(replica_group.keys()):
                dataset = replica_group[key]
                print(f"    {key}: shape={dataset.shape}, dtype={dataset.dtype}")
                
                if detailed:
                    try:
                        print("      Sample Data:")
                        print_dataset_sample(dataset, indent=8)
                    except Exception as e:
                        print(f"      Error analyzing dataset: {e}")
    
    print("\n  ...")
    print("  (Output limited for brevity. Use --full-structure for complete exploration)")


def extract_rmsf(h5_file: h5py.File, domain_id: str, temperature: str, replica: str) -> Optional[pd.DataFrame]:
    """
    Extract RMSF data for a specific domain, temperature, and replica.
    
    Returns a DataFrame with RMSF data or None if extraction fails.
    """
    try:
        if domain_id not in h5_file:
            print(f"Domain {domain_id} not found.")
            return None
        
        domain_group = h5_file[domain_id]
        
        if temperature not in domain_group:
            print(f"Temperature {temperature} not found for domain {domain_id}.")
            return None
        
        temp_group = domain_group[temperature]
        
        if replica not in temp_group:
            print(f"Replica {replica} not found for domain {domain_id}, temperature {temperature}.")
            return None
        
        replica_group = temp_group[replica]
        
        if 'rmsf' not in replica_group:
            print(f"RMSF data not found for domain {domain_id}, temperature {temperature}, replica {replica}.")
            return None
        
        # Get dataset and check size
        rmsf_dataset = replica_group['rmsf']
        rmsf_data = rmsf_dataset[:]
        
        # Extract residue information
        resid_dataset = domain_group['resid']
        resname_dataset = domain_group['resname']
        
        resids = resid_dataset[:]
        resnames = [name.decode('utf-8') if isinstance(name, bytes) else str(name) for name in resname_dataset[:]]
        
        # Handle dimension mismatch - RMSF is typically per residue (not per atom)
        if len(resids) != len(rmsf_data):
            print(f"Dimension mismatch: resids {len(resids)}, rmsf_data {len(rmsf_data)}")
            
            # Get unique residue IDs with their first occurrence index
            unique_resids = []
            unique_resnames = []
            resid_indices = {}
            
            for i, resid in enumerate(resids):
                if resid not in resid_indices:
                    resid_indices[resid] = i
                    unique_resids.append(resid)
                    unique_resnames.append(resnames[i])
            
            if len(unique_resids) == len(rmsf_data):
                print(f"Using unique residue IDs for RMSF data alignment")
                # Use unique residues for RMSF data
                resids = unique_resids
                resnames = unique_resnames
            elif len(unique_resids) > len(rmsf_data) and len(rmsf_data) > 0:
                # If we have more unique residues than RMSF data points,
                # try to match the first N residues
                print(f"More unique residues ({len(unique_resids)}) than RMSF data points ({len(rmsf_data)})")
                resids = unique_resids[:len(rmsf_data)]
                resnames = unique_resnames[:len(rmsf_data)]
            elif len(unique_resids) < len(rmsf_data) and len(unique_resids) > 0:
                # If we have fewer unique residues than RMSF data points,
                # truncate the RMSF data
                print(f"Fewer unique residues ({len(unique_resids)}) than RMSF data points ({len(rmsf_data)})")
                rmsf_data = rmsf_data[:len(unique_resids)]
            else:
                # If things are really mismatched, fall back to using the smaller length
                min_len = min(len(unique_resids), len(rmsf_data))
                if min_len == 0:
                    print(f"No valid data for alignment between resids and RMSF")
                    return None
                
                print(f"Truncating data to length {min_len}")
                resids = unique_resids[:min_len]
                resnames = unique_resnames[:min_len]
                rmsf_data = rmsf_data[:min_len]
        
        # Create DataFrame
        df = pd.DataFrame({
            'domain_id': domain_id,
            'resid': resids,
            'resname': resnames,
            f'rmsf_{temperature}': rmsf_data
        })
        
        return df
    except Exception as e:
        print(f"Failed to extract RMSF data: {e}")
        return None


def extract_dssp(h5_file: h5py.File, domain_id: str, temperature: str, replica: str, frame: int = -1) -> Optional[pd.DataFrame]:
    """
    Extract DSSP data for a specific domain, temperature, replica, and frame.
    
    Returns a DataFrame with DSSP data or None if extraction fails.
    """
    try:
        if domain_id not in h5_file:
            print(f"Domain {domain_id} not found.")
            return None
        
        domain_group = h5_file[domain_id]
        
        if temperature not in domain_group:
            print(f"Temperature {temperature} not found for domain {domain_id}.")
            return None
        
        temp_group = domain_group[temperature]
        
        if replica not in temp_group:
            print(f"Replica {replica} not found for domain {domain_id}, temperature {temperature}.")
            return None
        
        replica_group = temp_group[replica]
        
        if 'dssp' not in replica_group:
            print(f"DSSP data not found for domain {domain_id}, temperature {temperature}, replica {replica}.")
            return None
        
        # Extract DSSP data
        dssp_dataset = replica_group['dssp']
        
        # Handle the case where frame is out of bounds
        num_frames = dssp_dataset.shape[0] if len(dssp_dataset.shape) > 0 else 0
        if num_frames == 0:
            print(f"Empty DSSP dataset for domain {domain_id}, temperature {temperature}, replica {replica}")
            return None
        
        if frame < 0:
            # Convert negative indices to positive
            frame = num_frames + frame
        
        if frame < 0 or frame >= num_frames:
            print(f"Frame index {frame} out of bounds (0-{num_frames-1}) for domain {domain_id}")
            frame = min(max(0, frame), num_frames - 1)  # Clamp to valid range
        
        dssp_data = dssp_dataset[frame]
        
        # Extract residue information
        resids = domain_group['resid'][:]
        resnames = [name.decode('utf-8') if isinstance(name, bytes) else str(name) for name in domain_group['resname'][:]]
        
        # Decode DSSP codes
        dssp_codes = [code.decode('utf-8') if isinstance(code, bytes) else str(code) for code in dssp_data]
        
        # Handle dimension mismatch - DSSP is per residue, not per atom
        if len(resids) != len(dssp_codes):
            print(f"Dimension mismatch in DSSP: resids {len(resids)}, dssp_codes {len(dssp_codes)}")
            
            # Get unique residue IDs with their first occurrence index
            unique_resids = []
            unique_resnames = []
            resid_indices = {}
            
            for i, resid in enumerate(resids):
                if resid not in resid_indices:
                    resid_indices[resid] = i
                    unique_resids.append(resid)
                    unique_resnames.append(resnames[i])
            
            if len(unique_resids) == len(dssp_codes):
                print(f"Using unique residue IDs for DSSP data alignment")
                # Use unique residues for DSSP data
                resids = unique_resids
                resnames = unique_resnames
            elif len(unique_resids) > len(dssp_codes) and len(dssp_codes) > 0:
                # If we have more unique residues than DSSP data points,
                # try to match the first N residues
                print(f"More unique residues ({len(unique_resids)}) than DSSP data points ({len(dssp_codes)})")
                resids = unique_resids[:len(dssp_codes)]
                resnames = unique_resnames[:len(dssp_codes)]
            elif len(unique_resids) < len(dssp_codes) and len(unique_resids) > 0:
                # If we have fewer unique residues than DSSP data points,
                # truncate the DSSP data
                print(f"Fewer unique residues ({len(unique_resids)}) than DSSP data points ({len(dssp_codes)})")
                dssp_codes = dssp_codes[:len(unique_resids)]
            else:
                # If things are really mismatched, fall back to using the smaller length
                min_len = min(len(unique_resids), len(dssp_codes))
                if min_len == 0:
                    print(f"No valid data for alignment between resids and DSSP")
                    return None
                
                print(f"Truncating DSSP data to length {min_len}")
                resids = unique_resids[:min_len]
                resnames = unique_resnames[:min_len]
                dssp_codes = dssp_codes[:min_len]
        
        # Create DataFrame
        df = pd.DataFrame({
            'domain_id': domain_id,
            'resid': resids,
            'resname': resnames,
            'dssp': dssp_codes
        })
        
        return df
    except Exception as e:
        print(f"Failed to extract DSSP data: {e}")
        return None


def extract_coordinates(h5_file: h5py.File, domain_id: str, temperature: str, replica: str, frame: int = -1) -> Optional[Tuple[np.ndarray, List[int], List[str]]]:
    """
    Extract coordinate data for a specific domain, temperature, replica, and frame.
    
    Returns a tuple of (coordinates, residue IDs, residue names) or None if extraction fails.
    """
    try:
        if domain_id not in h5_file:
            print(f"Domain {domain_id} not found.")
            return None
        
        domain_group = h5_file[domain_id]
        
        if temperature not in domain_group:
            print(f"Temperature {temperature} not found for domain {domain_id}.")
            return None
        
        temp_group = domain_group[temperature]
        
        if replica not in temp_group:
            print(f"Replica {replica} not found for domain {domain_id}, temperature {temperature}.")
            return None
        
        replica_group = temp_group[replica]
        
        if 'coords' not in replica_group:
            print(f"Coordinate data not found for domain {domain_id}, temperature {temperature}, replica {replica}.")
            return None
        
        # Extract coordinate data for specified frame
        coords_dataset = replica_group['coords']
        
        # Handle the case where frame is out of bounds
        num_frames = coords_dataset.shape[0] if len(coords_dataset.shape) > 0 else 0
        if num_frames == 0:
            print(f"Empty coordinates dataset for domain {domain_id}, temperature {temperature}, replica {replica}")
            return None
        
        if frame < 0:
            # Convert negative indices to positive
            frame = num_frames + frame
        
        if frame < 0 or frame >= num_frames:
            print(f"Frame index {frame} out of bounds (0-{num_frames-1}) for domain {domain_id}")
            frame = min(max(0, frame), num_frames - 1)  # Clamp to valid range
        
        coords = coords_dataset[frame]
        
        # Extract residue information
        resids = domain_group['resid'][:].tolist()
        resnames = [name.decode('utf-8') if isinstance(name, bytes) else str(name) for name in domain_group['resname'][:]]
        
        return coords, resids, resnames
    except Exception as e:
        print(f"Failed to extract coordinate data: {e}")
        return None


def extract_pdb(h5_file: h5py.File, domain_id: str) -> Optional[str]:
    """
    Extract PDB data for a specific domain.
    
    Returns PDB data as a string or None if extraction fails.
    """
    try:
        if domain_id not in h5_file:
            print(f"Domain {domain_id} not found.")
            return None
        
        domain_group = h5_file[domain_id]
        
        if 'pdb' not in domain_group:
            print(f"PDB data not found for domain {domain_id}.")
            return None
        
        pdb_data = domain_group['pdb'][()]
        if isinstance(pdb_data, bytes):
            return pdb_data.decode('utf-8')
        return str(pdb_data)
    except Exception as e:
        print(f"Failed to extract PDB data: {e}")
        return None


def compare_domains(h5_file: h5py.File, domain_id1: str, domain_id2: str):
    """Compare two domains side by side."""
    if domain_id1 not in h5_file:
        print(f"Domain {domain_id1} not found in the file.")
        return
    
    if domain_id2 not in h5_file:
        print(f"Domain {domain_id2} not found in the file.")
        return
    
    domain1 = h5_file[domain_id1]
    domain2 = h5_file[domain_id2]
    
    print_section(f"Comparing Domains: {domain_id1} vs {domain_id2}")
    
    # Compare metadata
    print_subsection("Domain Metadata Comparison")
    
    # Get domain-level datasets
    domain1_datasets = {key: domain1[key] for key in domain1.keys() if isinstance(domain1[key], h5py.Dataset)}
    domain2_datasets = {key: domain2[key] for key in domain2.keys() if isinstance(domain2[key], h5py.Dataset)}
    
    all_keys = sorted(set(list(domain1_datasets.keys()) + list(domain2_datasets.keys())))
    
    print("Domain-level Datasets:")
    data = []
    headers = ["Dataset", f"{domain_id1} Shape", f"{domain_id2} Shape"]
    
    for key in all_keys:
        shape1 = domain1_datasets[key].shape if key in domain1_datasets else "N/A"
        shape2 = domain2_datasets[key].shape if key in domain2_datasets else "N/A"
        data.append([key, shape1, shape2])
    
    print(tabulate(data, headers=headers, tablefmt="grid"))
    
    # Compare temperature groups
    print_subsection("Temperature Groups Comparison")
    
    temps1 = set()
    temps2 = set()
    
    for key in domain1.keys():
        try:
            temp = int(key)
            temps1.add(temp)
        except ValueError:
            continue
    
    for key in domain2.keys():
        try:
            temp = int(key)
            temps2.add(temp)
        except ValueError:
            continue
    
    print(f"Temperatures in {domain_id1}: {', '.join(map(str, sorted(temps1)))}")
    print(f"Temperatures in {domain_id2}: {', '.join(map(str, sorted(temps2)))}")
    print(f"Common temperatures: {', '.join(map(str, sorted(temps1 & temps2)))}")
    
    # Compare RMSF data for a common temperature and replica
    common_temps = sorted(list(temps1 & temps2))
    if common_temps:
        first_temp = str(common_temps[0])
        
        # Find common replicas
        replicas1 = set(domain1[first_temp].keys())
        replicas2 = set(domain2[first_temp].keys())
        common_replicas = sorted(list(replicas1 & replicas2))
        
        if common_replicas:
            first_replica = common_replicas[0]
            
            print_subsection(f"RMSF Data Comparison (Temperature: {first_temp}, Replica: {first_replica})")
            
            rmsf1 = extract_rmsf(h5_file, domain_id1, first_temp, first_replica)
            rmsf2 = extract_rmsf(h5_file, domain_id2, first_temp, first_replica)
            
            if rmsf1 is not None and rmsf2 is not None:
                print(f"\n{domain_id1} RMSF Sample (first 5 rows):")
                print(rmsf1.head().to_string())
                
                print(f"\n{domain_id2} RMSF Sample (first 5 rows):")
                print(rmsf2.head().to_string())
                
                # Compare statistics
                print("\nRMSF Statistics Comparison:")
                print(f"{domain_id1} - Count: {len(rmsf1)}, Min: {rmsf1[f'rmsf_{first_temp}'].min():.4f}, "
                      f"Max: {rmsf1[f'rmsf_{first_temp}'].max():.4f}, Mean: {rmsf1[f'rmsf_{first_temp}'].mean():.4f}")
                print(f"{domain_id2} - Count: {len(rmsf2)}, Min: {rmsf2[f'rmsf_{first_temp}'].min():.4f}, "
                      f"Max: {rmsf2[f'rmsf_{first_temp}'].max():.4f}, Mean: {rmsf2[f'rmsf_{first_temp}'].mean():.4f}")
            else:
                print("Could not extract RMSF data for comparison.")
            
            # Also compare DSSP data
            print_subsection(f"DSSP Data Comparison (Temperature: {first_temp}, Replica: {first_replica})")
            
            dssp1 = extract_dssp(h5_file, domain_id1, first_temp, first_replica)
            dssp2 = extract_dssp(h5_file, domain_id2, first_temp, first_replica)
            
            if dssp1 is not None and dssp2 is not None:
                print(f"\n{domain_id1} DSSP Sample (first 5 rows):")
                print(dssp1.head().to_string())
                
                print(f"\n{domain_id2} DSSP Sample (first 5 rows):")
                print(dssp2.head().to_string())
                
                # Compare secondary structure distribution
                print("\nSecondary Structure Distribution:")
                
                ss_dist1 = dssp1['dssp'].value_counts().sort_index()
                ss_dist2 = dssp2['dssp'].value_counts().sort_index()
                
                # Fill missing values with 0
                all_ss = sorted(set(list(ss_dist1.index) + list(ss_dist2.index)))
                for ss in all_ss:
                    if ss not in ss_dist1:
                        ss_dist1[ss] = 0
                    if ss not in ss_dist2:
                        ss_dist2[ss] = 0
                
                ss_dist1 = ss_dist1.sort_index()
                ss_dist2 = ss_dist2.sort_index()
                
                data = []
                for ss in all_ss:
                    data.append([ss, ss_dist1.get(ss, 0), ss_dist2.get(ss, 0)])
                
                print(tabulate(data, headers=["DSSP Code", domain_id1, domain_id2], tablefmt="grid"))
            else:
                print("Could not extract DSSP data for comparison.")
            
            # Compare coordinate data
            print_subsection(f"Coordinate Data Comparison (Temperature: {first_temp}, Replica: {first_replica})")
            
            coords1 = extract_coordinates(h5_file, domain_id1, first_temp, first_replica)
            coords2 = extract_coordinates(h5_file, domain_id2, first_temp, first_replica)
            
            if coords1 is not None and coords2 is not None:
                coords_data1, resids1, resnames1 = coords1
                coords_data2, resids2, resnames2 = coords2
                
                print(f"{domain_id1} Coordinates Shape: {coords_data1.shape}")
                print(f"{domain_id2} Coordinates Shape: {coords_data2.shape}")
                
                print(f"\n{domain_id1} First 3 Coordinates:")
                for i in range(min(3, len(resids1))):
                    print(f"  Residue {resids1[i]} ({resnames1[i]}): {coords_data1[i]}")
                
                print(f"\n{domain_id2} First 3 Coordinates:")
                for i in range(min(3, len(resids2))):
                    print(f"  Residue {resids2[i]} ({resnames2[i]}): {coords_data2[i]}")
            else:
                print("Could not extract coordinate data for comparison.")
        else:
            print(f"No common replicas found for temperature {first_temp}.")
    else:
        print("No common temperatures found for comparison.")


def summarize_file(h5_file: h5py.File):
    """Print a summary of the HDF5 file."""
    print_section("mdCATH .h5 File Summary")
    
    # Get all domain IDs
    domain_ids = list(h5_file.keys())
    num_domains = len(domain_ids)
    print(f"Number of domains: {num_domains}")
    
    if num_domains == 0:
        print("No domains found in the file.")
        return
    
    # Print a few sample domains
    max_samples = min(5, num_domains)
    print(f"\nSample domains: {', '.join(domain_ids[:max_samples])}" + 
          (f"... ({num_domains-max_samples} more)" if num_domains > max_samples else ""))
    
    # Get size information
    domain_sizes = {}
    for domain_id in domain_ids[:10]:  # Limit to 10 domains to avoid excessive processing
        domain_size = 0
        
        def get_size(name, obj):
            nonlocal domain_size
            if isinstance(obj, h5py.Dataset):
                domain_size += obj.size * obj.dtype.itemsize
        
        h5_file[domain_id].visititems(get_size)
        domain_sizes[domain_id] = domain_size
    
    print("\nDomain Sizes (bytes):")
    for domain_id, size in domain_sizes.items():
        print(f"  {domain_id}: {size:,} bytes ({size / (1024*1024):.2f} MB)")
    
    # Get available temperatures
    temperatures = set()
    for domain_id in domain_ids:
        domain_group = h5_file[domain_id]
        for key in domain_group.keys():
            try:
                temp = int(key)
                temperatures.add(temp)
            except ValueError:
                continue
    
    if temperatures:
        print(f"\nAvailable temperatures: {', '.join(map(str, sorted(temperatures)))}")
    else:
        print("\nNo temperature data found.")
    
    # Check for common domain datasets
    common_datasets = ["chain", "element", "resid", "resname", "pdb", "z"]
    
    print("\nCommon domain datasets:")
    dataset_presence = {ds: 0 for ds in common_datasets}
    
    for domain_id in domain_ids:
        domain_group = h5_file[domain_id]
        for ds_name in common_datasets:
            if ds_name in domain_group:
                dataset_presence[ds_name] += 1
    
    for ds_name, count in dataset_presence.items():
        percentage = (count / num_domains) * 100
        print(f"  {ds_name}: Found in {count}/{num_domains} domains ({percentage:.1f}%)")
    
    # Sample the first domain for detailed dataset info
    first_domain = h5_file[domain_ids[0]]
    print("\nDataset details for first domain:")
    for ds_name in common_datasets:
        if ds_name in first_domain:
            dataset = first_domain[ds_name]
            print(f"  {ds_name}: shape={dataset.shape}, dtype={dataset.dtype}")
            
            # Print sample data for smaller datasets
            if dataset.size < 100:
                try:
                    if dataset.dtype.kind == 'S':
                        sample = [s.decode('utf-8') if isinstance(s, bytes) else s for s in dataset[:]]
                    else:
                        sample = dataset[:]
                    print(f"    Sample: {sample}")
                except:
                    print("    Error showing sample")
            else:
                try:
                    if dataset.dtype.kind == 'S':
                        sample = [s.decode('utf-8') if isinstance(s, bytes) else s for s in dataset[:5]]
                    else:
                        sample = dataset[:5]
                    print(f"    Sample (first 5): {sample}...")
                except:
                    print("    Error showing sample")
        else:
            print(f"  {ds_name}: Not found")
    
    # Check for common attributes
    print("\nCommon domain attributes:")
    if len(first_domain.attrs) > 0:
        for attr_name, attr_value in first_domain.attrs.items():
            print(f"  {attr_name}: {attr_value}")
    else:
        print("  No attributes found.")
    
    # Sample structure for a temperature and replica
    first_domain = h5_file[domain_ids[0]]
    print("\nSample temperature and replica structure:")
    
    temp_groups = []
    for key in first_domain.keys():
        try:
            temp = int(key)
            temp_groups.append(key)
        except ValueError:
            continue
    
    if temp_groups:
        first_temp = sorted(temp_groups)[0]
        print(f"  Temperature: {first_temp}")
        
        temp_group = first_domain[first_temp]
        replica_groups = list(temp_group.keys())
        
        if replica_groups:
            first_replica = sorted(replica_groups)[0]
            print(f"  Replica: {first_replica}")
            
            replica_group = temp_group[first_replica]
            print("  Datasets:")
            
            for key in replica_group.keys():
                dataset = replica_group[key]
                print(f"    {key}: shape={dataset.shape}, dtype={dataset.dtype}")
        else:
            print("  No replica groups found.")
    else:
        print("  No temperature groups found.")


def main():
    parser = argparse.ArgumentParser(description="Enhanced inspection tool for mdCATH .h5 files")
    parser.add_argument("file_path", help="Path to the mdCATH .h5 file")
    parser.add_argument("--domain", help="Specific domain to inspect")
    parser.add_argument("--compare", nargs=2, metavar=("DOMAIN1", "DOMAIN2"), 
                        help="Compare two domains side by side")
    parser.add_argument("--full-structure", action="store_true", 
                        help="Print the full file structure")
    parser.add_argument("--max-depth", type=int, default=3, 
                        help="Maximum depth to explore in the file structure")
    parser.add_argument("--detailed", action="store_true", 
                        help="Show detailed information including data samples")
    parser.add_argument("--extract-rmsf", nargs=3, metavar=("DOMAIN", "TEMP", "REPLICA"), 
                        help="Extract RMSF data for a specific domain, temperature, and replica")
    parser.add_argument("--extract-dssp", nargs=3, metavar=("DOMAIN", "TEMP", "REPLICA"), 
                        help="Extract DSSP data for a specific domain, temperature, and replica")
    parser.add_argument("--extract-coords", nargs=3, metavar=("DOMAIN", "TEMP", "REPLICA"), 
                        help="Extract coordinate data for a specific domain, temperature, and replica")
    parser.add_argument("--extract-pdb", metavar="DOMAIN", 
                        help="Extract PDB data for a specific domain")

    args = parser.parse_args()

    if not os.path.exists(args.file_path):
        print(f"Error: File not found: {args.file_path}")
        sys.exit(1)

    try:
        with h5py.File(args.file_path, 'r') as h5_file:
            if args.domain:
                print_domain_info(h5_file, args.domain, args.detailed)
            elif args.compare:
                compare_domains(h5_file, args.compare[0], args.compare[1])
            elif args.full_structure:
                print_section("HDF5 File Structure")
                explore_group(h5_file, args.max_depth, details=args.detailed)
            elif args.extract_rmsf:
                domain, temp, replica = args.extract_rmsf
                df = extract_rmsf(h5_file, domain, temp, replica)
                if df is not None:
                    print(f"RMSF data for domain {domain}, temperature {temp}, replica {replica}:")
                    print(df.head(10).to_string())
                    print(f"\nTotal rows: {len(df)}")
            elif args.extract_dssp:
                domain, temp, replica = args.extract_dssp
                df = extract_dssp(h5_file, domain, temp, replica)
                if df is not None:
                    print(f"DSSP data for domain {domain}, temperature {temp}, replica {replica}:")
                    print(df.head(10).to_string())
                    print(f"\nTotal rows: {len(df)}")
            elif args.extract_coords:
                domain, temp, replica = args.extract_coords
                result = extract_coordinates(h5_file, domain, temp, replica)
                if result is not None:
                    coords, resids, resnames = result
                    print(f"Coordinate data for domain {domain}, temperature {temp}, replica {replica}:")
                    print(f"Shape: {coords.shape}")
                    print("\nFirst 10 residues:")
                    for i in range(min(10, len(resids))):
                        print(f"Residue {resids[i]} ({resnames[i]}): {coords[i]}")
            elif args.extract_pdb:
                pdb_data = extract_pdb(h5_file, args.extract_pdb)
                if pdb_data is not None:
                    print(f"PDB data for domain {args.extract_pdb}:")
                    lines = pdb_data.split('\n')
                    print(f"Total lines: {len(lines)}")
                    print("\nFirst 10 lines:")
                    for line in lines[:10]:
                        print(line)
            else:
                summarize_file(h5_file)
    except Exception as e:
        print(f"Error inspecting file: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()