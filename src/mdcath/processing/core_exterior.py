#!/usr/bin/env python3
"""
Processing module for core/exterior classification.
"""

import os
import logging
import subprocess
import tempfile
import pandas as pd
import numpy as np
import Bio
from Bio.PDB import PDBParser, DSSP, ShrakeRupley
from typing import Dict, Any, Optional, List, Tuple


def compute_core_exterior(pdb_file: str, config: Dict[str, Any]) -> Optional[pd.DataFrame]:
    """
    Classify residues as 'core' or 'exterior' based on solvent accessibility.

    Args:
        pdb_file: Path to the cleaned PDB file
        config: Configuration dictionary

    Returns:
        DataFrame with columns 'resid' and 'core_exterior' or None if classification fails
    """
    method = config.get("core_exterior", {}).get("method", "msms")

    if method == "msms":
        return compute_core_exterior_msms(pdb_file, config)
    else:
        return compute_core_exterior_biopython(pdb_file, config)

def compute_core_exterior_msms(pdb_file: str, config: Dict[str, Any]) -> Optional[pd.DataFrame]:
    """
    Use MSMS to classify residues as 'core' or 'exterior'.

    Args:
        pdb_file: Path to the cleaned PDB file
        config: Configuration dictionary

    Returns:
        DataFrame with columns 'resid' and 'core_exterior' or None if MSMS fails
    """
    msms_dir = config.get("core_exterior", {}).get("msms_executable_dir", "./msms")
    ses_threshold = config.get("core_exterior", {}).get("ses_threshold", 1.0)
    protein_name = os.path.basename(pdb_file).split('.')[0]

    try:
        # Create temporary directory for MSMS files
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Paths to MSMS executables and output files
            pdb2xyzr_exe = os.path.join(msms_dir, "pdb_to_xyzr")
            msms_exe = os.path.join(msms_dir, "msms.x86_64Linux2.2.6.1")
            xyzr_file = os.path.join(tmp_dir, f"{protein_name}.xyzr")
            area_base = os.path.join(tmp_dir, f"{protein_name}")
            area_file = os.path.join(tmp_dir, f"{protein_name}.area")

            # Check MSMS executables
            if not os.path.exists(pdb2xyzr_exe) or not os.path.exists(msms_exe):
                logging.warning(f"MSMS executables not found in {msms_dir}, falling back to Biopython")
                return compute_core_exterior_biopython(pdb_file, config)

            # Run pdb_to_xyzr
            cmd_xyzr = f"{pdb2xyzr_exe} {pdb_file} > {xyzr_file}"
            result = subprocess.run(cmd_xyzr, shell=True, check=False,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if result.returncode != 0 or not os.path.exists(xyzr_file) or os.path.getsize(xyzr_file) == 0:
                logging.warning(f"pdb_to_xyzr failed: {result.stderr.decode()}, falling back to Biopython")
                return compute_core_exterior_biopython(pdb_file, config)

            # Run MSMS
            cmd_msms = f"{msms_exe} -if {xyzr_file} -af {area_base}"
            result = subprocess.run(cmd_msms, shell=True, check=False,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if result.returncode != 0 or not os.path.exists(area_file):
                logging.warning(f"MSMS failed: {result.stderr.decode()}, falling back to Biopython")
                return compute_core_exterior_biopython(pdb_file, config)

            # Parse atom-level PDB data
            per_atom_df = parse_pdb_atoms(pdb_file)
            if per_atom_df.empty:
                logging.warning(f"Failed to parse atoms from PDB, falling back to Biopython")
                return compute_core_exterior_biopython(pdb_file, config)

            # Parse MSMS area file
            area_df = parse_area_file(area_file)
            if area_df.empty:
                logging.warning(f"Failed to parse area file, falling back to Biopython")
                return compute_core_exterior_biopython(pdb_file, config)

            # Combine atom data with MSMS results
            if len(area_df) != len(per_atom_df):
                logging.warning(f"Atom count mismatch: {len(area_df)} vs {len(per_atom_df)}, falling back to Biopython")
                return compute_core_exterior_biopython(pdb_file, config)

            # Merge data
            per_atom_df = pd.concat([per_atom_df.reset_index(drop=True),
                                    area_df.reset_index(drop=True)], axis=1)

            # Calculate mean SES per residue
            mean_ses_per_res = per_atom_df.groupby("resid")["SES"].mean()

            # Classify residues as core or exterior
            exterior_residues = mean_ses_per_res[mean_ses_per_res > ses_threshold].index
            resids = mean_ses_per_res.index.tolist()
            core_exterior = ["exterior" if r in exterior_residues else "core" for r in resids]

            # Create final dataframe
            result_df = pd.DataFrame({
                "resid": resids,
                "core_exterior": core_exterior
            })

            return result_df
    except Exception as e:
        logging.warning(f"MSMS processing failed: {e}, falling back to Biopython")
        return compute_core_exterior_biopython(pdb_file, config)
# In src/mdcath/processing/core_exterior.py, let's improve the Biopython fallback

def compute_core_exterior_biopython(pdb_file: str, config: Dict[str, Any]) -> pd.DataFrame:
    """
    Use Biopython's SASA calculation to classify residues as 'core' or 'exterior'.

    Args:
        pdb_file: Path to the cleaned PDB file
        config: Configuration dictionary

    Returns:
        DataFrame with columns 'resid' and 'core_exterior'
    """
    try:
        from Bio.PDB import PDBParser, Selection, DSSP
        from Bio.PDB.SASA import ShrakeRupley

        # Set SASA threshold
        sasa_threshold = config.get("core_exterior", {}).get("sasa_threshold", 20.0)

        # Parse PDB
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        model = structure[0]
        
        # Try DSSP first for better solvent accessibility calculation
        try:
            # Get the location of the DSSP executable
            # If it's in your path, this should be fine
            dssp_executable = shutil.which("dssp") or shutil.which("mkdssp")
            if dssp_executable:
                logging.info(f"Using DSSP executable: {dssp_executable}")
                # Run DSSP
                dssp = DSSP(model, pdb_file, dssp=dssp_executable)
                
                # Extract results - DSSP gives accessibility values directly
                results = []
                for chain in model:
                    for residue in chain:
                        if residue.id[0] == " ":  # Standard residue
                            resid = residue.id[1]
                            try:
                                # DSSP key is (chain_id, residue_id)
                                key = (chain.id, (' ', residue.id[1], ' '))
                                if key in dssp:
                                    # Get relative solvent accessibility (values range from 0-1)
                                    rel_acc = dssp[key][3]
                                    # Convert to absolute accessibility using max values for residue types
                                    # A value > 0.2 (20%) is generally considered accessible
                                    core_exterior = "exterior" if rel_acc > 0.2 else "core"
                                    results.append({"resid": resid, "core_exterior": core_exterior})
                                else:
                                    # If residue not found in DSSP, use default classification
                                    results.append({"resid": resid, "core_exterior": "unknown"})
                            except Exception as e:
                                logging.warning(f"Error processing DSSP for residue {resid}: {e}")
                                results.append({"resid": resid, "core_exterior": "unknown"})
                
                if results:
                    logging.info("Successfully used DSSP for core/exterior classification")
                    return pd.DataFrame(results)
                
            # If DSSP fails or no results, fall back to ShrakeRupley
            logging.info("DSSP failed or not available, falling back to ShrakeRupley SASA")
        
        except Exception as e:
            logging.warning(f"DSSP calculation failed: {e}, falling back to ShrakeRupley")
        
        # Fall back to ShrakeRupley SASA calculation
        sr = ShrakeRupley()
        sr.compute(model, level="R")  # Compute at residue level

        # Extract results
        results = []
        for chain in model:
            for residue in chain:
                if residue.id[0] == " ":  # Standard residue
                    resid = residue.id[1]
                    sasa = residue.sasa if hasattr(residue, 'sasa') else 0.0
                    core_exterior = "exterior" if sasa > sasa_threshold else "core"
                    results.append({"resid": resid, "core_exterior": core_exterior})

        return pd.DataFrame(results)
    except Exception as e:
        logging.error(f"Biopython SASA calculation failed: {e}")
        import traceback
        logging.error(traceback.format_exc())
        return fallback_core_exterior(pdb_file)

def fallback_core_exterior(pdb_file: str) -> pd.DataFrame:
    """
    Fallback method to classify residues when other methods fail.
    Classifies outer 1/3 of residues as exterior, inner 2/3 as core.

    Args:
        pdb_file: Path to the cleaned PDB file

    Returns:
        DataFrame with columns 'resid' and 'core_exterior'
    """
    try:
        # Parse PDB to get residue information
        residue_df = parse_pdb_residues(pdb_file)
        if residue_df.empty:
            # Create empty DataFrame with required columns
            return pd.DataFrame(columns=["resid", "core_exterior"])

        # Sort by residue ID
        residue_df = residue_df.sort_values("resid")

        # Simple classification: outer 1/3 of residues as exterior, inner 2/3 as core
        total_residues = len(residue_df)
        boundary = int(total_residues * 2/3)

        residue_df["core_exterior"] = ["core"] * total_residues
        residue_df.loc[boundary:, "core_exterior"] = "exterior"

        return residue_df[["resid", "core_exterior"]]
    except Exception as e:
        logging.error(f"Fallback classification failed: {e}")
        return pd.DataFrame(columns=["resid", "core_exterior"])

def parse_pdb_residues(pdb_file: str) -> pd.DataFrame:
    """
    Parse a PDB file to extract residue-level information.

    Args:
        pdb_file: Path to the PDB file

    Returns:
        DataFrame with residue information
    """
    try:
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)

        records = []
        for model in structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    if residue.id[0] == " ":  # Standard residue
                        records.append({
                            "resid": residue.id[1],
                            "resname": residue.get_resname(),
                            "chain": chain_id
                        })

        return pd.DataFrame(records)
    except Exception as e:
        logging.error(f"Failed to parse PDB residues: {e}")
        return pd.DataFrame()

def parse_pdb_atoms(pdb_file: str) -> pd.DataFrame:
    """
    Parse a PDB file to extract atom-level information.

    Args:
        pdb_file: Path to the PDB file

    Returns:
        DataFrame with atom information
    """
    try:
        from Bio.PDB import PDBParser

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)

        records = []
        atom_idx = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == " ":  # Standard residue
                        res_id = residue.id[1]
                        res_name = residue.get_resname()
                        for atom in residue:
                            atom_idx += 1
                            records.append({
                                "atom_idx": atom_idx,
                                "resid": res_id,
                                "resname": res_name,
                                "atom_name": atom.get_name()
                            })

        return pd.DataFrame(records)
    except Exception as e:
        logging.error(f"Failed to parse PDB atoms: {e}")
        return pd.DataFrame()

def parse_area_file(area_file: str) -> pd.DataFrame:
    """
    Parse an MSMS .area file to extract SES values per atom.

    Args:
        area_file: Path to the MSMS .area file

    Returns:
        DataFrame with SES values
    """
    try:
        atom_idx = []
        ses = []

        with open(area_file, "r") as f:
            for line in f:
                if "Atom" in line or not line.strip():
                    continue

                cols = line.split()
                if len(cols) >= 2:
                    atom_idx.append(int(cols[0]))
                    ses.append(float(cols[1]))

        return pd.DataFrame({"atom_idx": atom_idx, "SES": ses})
    except Exception as e:
        logging.error(f"Failed to parse area file: {e}")
        return pd.DataFrame()
