#!/usr/bin/env python3
"""
Convert PDB structures to PDBQT format for AutoDock Vina.

This script:
1. Converts the receptor PDB to PDBQT
2. Converts all peptide variant PDBs to PDBQT (if they exist)

Requires: OpenBabel installed (conda install -c conda-forge openbabel)
"""

import os
import subprocess
import logging
from pathlib import Path

from fetch_data import DATA_DIR

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

RECEPTOR_PDB = os.path.join(DATA_DIR, "P43250_receptor_clean.pdb")
RECEPTOR_PDBQT = os.path.join(DATA_DIR, "P43250_receptor_clean.pdbqt")
DOCKING_ROOT = os.path.join(DATA_DIR, "docking")


def check_openbabel():
    """Check if OpenBabel is installed."""
    try:
        result = subprocess.run(
            ["obabel", "-V"],
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def convert_pdb_to_pdbqt(pdb_path, pdbqt_path):
    """
    Convert PDB file to PDBQT format using OpenBabel.
    
    Args:
        pdb_path: Input PDB file path
        pdbqt_path: Output PDBQT file path
    """
    if not os.path.exists(pdb_path):
        logging.warning(f"Input file not found: {pdb_path}")
        return False
    
    try:
        # -xr flag: remove non-polar hydrogens (needed for Vina)
        cmd = ["obabel", pdb_path, "-O", pdbqt_path, "-xr"]
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode == 0 and os.path.exists(pdbqt_path):
            logging.info(f"Converted: {pdb_path} -> {pdbqt_path}")
            return True
        else:
            logging.error(f"Conversion failed: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        logging.error(f"Conversion timed out for {pdb_path}")
        return False
    except Exception as e:
        logging.error(f"Error converting {pdb_path}: {e}")
        return False


def convert_receptor():
    """Convert receptor PDB to PDBQT."""
    if not os.path.exists(RECEPTOR_PDB):
        logging.error(f"Receptor PDB not found: {RECEPTOR_PDB}")
        logging.info("Run prepare_structures.py first to generate it.")
        return False
    
    return convert_pdb_to_pdbqt(RECEPTOR_PDB, RECEPTOR_PDBQT)


def convert_peptide_variants():
    """Convert all peptide variant PDBs to PDBQT."""
    if not os.path.exists(DOCKING_ROOT):
        logging.warning(f"Docking directory not found: {DOCKING_ROOT}")
        logging.info("Run docking_prep.py first to create variant directories.")
        return False
    
    converted = 0
    skipped = 0
    
    for variant_dir in Path(DOCKING_ROOT).glob("GRK6_variant_*"):
        if not variant_dir.is_dir():
            continue
            
        ligand_pdb = variant_dir / "ligand.pdb"
        ligand_pdbqt = variant_dir / "ligand.pdbqt"
        
        if not ligand_pdb.exists():
            logging.debug(f"No PDB file found for {variant_dir.name}, skipping")
            skipped += 1
            continue
        
        if convert_pdb_to_pdbqt(str(ligand_pdb), str(ligand_pdbqt)):
            converted += 1
        else:
            skipped += 1
    
    logging.info(f"Converted {converted} peptide variants, skipped {skipped}")
    return converted > 0


def main():
    """Main conversion workflow."""
    logging.info("=" * 60)
    logging.info("PDB to PDBQT Conversion")
    logging.info("=" * 60)
    
    # Check for OpenBabel
    if not check_openbabel():
        logging.error("OpenBabel not found!")
        logging.info("Install with: conda install -c conda-forge openbabel")
        logging.info("Or download from: https://openbabel.org")
        return False
    
    logging.info("OpenBabel found ✓")
    
    # Convert receptor
    logging.info("\n[Step 1/2] Converting receptor...")
    receptor_success = convert_receptor()
    
    # Convert peptides
    logging.info("\n[Step 2/2] Converting peptide variants...")
    peptides_success = convert_peptide_variants()
    
    if receptor_success:
        logging.info(f"\n✓ Receptor converted: {RECEPTOR_PDBQT}")
    else:
        logging.warning("\n✗ Receptor conversion failed or skipped")
    
    if peptides_success:
        logging.info("✓ Peptide variants converted")
    else:
        logging.warning("✗ No peptide PDB files found to convert")
        logging.info("  Generate 3D structures first (see NEXT_STEPS.md)")
    
    logging.info("\n" + "=" * 60)
    if receptor_success:
        logging.info("Conversion complete! Ready for docking.")
    else:
        logging.warning("Conversion incomplete. Check errors above.")
    logging.info("=" * 60)
    
    return receptor_success


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)






