#!/usr/bin/env python3
"""
Run AutoDock Vina docking for all peptide variants.

This script:
1. Checks that all required files exist
2. Runs Vina docking for each variant
3. Saves results to variant directories
"""

import os
import subprocess
import logging
from pathlib import Path

from fetch_data import DATA_DIR

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

RECEPTOR_PDBQT = os.path.join(DATA_DIR, "P43250_receptor_clean.pdbqt")
CONFIG_FILE = os.path.join(DATA_DIR, "docking_config.txt")
DOCKING_ROOT = os.path.join(DATA_DIR, "docking")


def check_vina():
    """Check if Vina is installed."""
    try:
        result = subprocess.run(
            ["vina", "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0 or "vina" in result.stdout.lower() or "vina" in result.stderr.lower()
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def create_default_config():
    """Create a default docking configuration file."""
    if os.path.exists(CONFIG_FILE):
        return
    
    logging.info(f"Creating default config file: {CONFIG_FILE}")
    with open(CONFIG_FILE, "w") as f:
        f.write(f"""receptor = {RECEPTOR_PDBQT}
center_x = 0.0
center_y = 0.0
center_z = 0.0
size_x = 20.0
size_y = 20.0
size_z = 20.0
exhaustiveness = 8
num_modes = 10
""")
    logging.warning("Default config created with center at (0,0,0). Update with actual binding site coordinates!")


def run_docking():
    """Run docking for all variants."""
    logging.info("=" * 70)
    logging.info("AUTODOCK VINA DOCKING")
    logging.info("=" * 70)
    
    # Check Vina
    if not check_vina():
        logging.error("AutoDock Vina not found!")
        logging.info("Install with: conda install -c conda-forge autodock-vina")
        return False
    
    logging.info("✓ AutoDock Vina found")
    
    # Check receptor
    if not os.path.exists(RECEPTOR_PDBQT):
        logging.error(f"Receptor PDBQT not found: {RECEPTOR_PDBQT}")
        logging.info("Run: python src/convert_to_pdbqt.py")
        return False
    
    logging.info(f"✓ Receptor found: {RECEPTOR_PDBQT}")
    
    # Create config if needed
    create_default_config()
    
    # Find variant directories
    if not os.path.exists(DOCKING_ROOT):
        logging.error(f"Docking directory not found: {DOCKING_ROOT}")
        logging.info("Run: python src/docking_prep.py")
        return False
    
    variant_dirs = sorted(Path(DOCKING_ROOT).glob("GRK6_variant_*"))
    if not variant_dirs:
        logging.error("No variant directories found")
        logging.info("Run: python src/docking_prep.py")
        return False
    
    logging.info(f"Found {len(variant_dirs)} variant directories")
    print()
    
    # Run docking for each variant
    success_count = 0
    fail_count = 0
    
    for variant_dir in variant_dirs:
        variant_id = variant_dir.name
        ligand_pdbqt = variant_dir / "ligand.pdbqt"
        output_pdbqt = variant_dir / "docked.pdbqt"
        log_file = variant_dir / "log.txt"
        
        if not ligand_pdbqt.exists():
            logging.warning(f"Skipping {variant_id}: ligand.pdbqt not found")
            fail_count += 1
            continue
        
        logging.info(f"Docking {variant_id}...")
        
        try:
            cmd = [
                "vina",
                "--receptor", RECEPTOR_PDBQT,
                "--ligand", str(ligand_pdbqt),
                "--config", CONFIG_FILE,
                "--out", str(output_pdbqt),
                "--log", str(log_file),
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout per variant
            )
            
            if result.returncode == 0:
                if log_file.exists() and output_pdbqt.exists():
                    logging.info(f"✓ {variant_id} complete")
                    success_count += 1
                else:
                    logging.warning(f"⚠ {variant_id} completed but output files missing")
                    fail_count += 1
            else:
                logging.error(f"✗ {variant_id} failed")
                logging.debug(f"Error: {result.stderr}")
                fail_count += 1
                
        except subprocess.TimeoutExpired:
            logging.error(f"✗ {variant_id} timed out (>10 minutes)")
            fail_count += 1
        except Exception as e:
            logging.error(f"✗ {variant_id} error: {e}")
            fail_count += 1
        
        print()
    
    # Summary
    logging.info("=" * 70)
    logging.info(f"Docking complete: {success_count} succeeded, {fail_count} failed")
    logging.info("=" * 70)
    
    if success_count > 0:
        logging.info("Next step: python src/parse_docking.py")
    
    return success_count > 0


if __name__ == "__main__":
    success = run_docking()
    exit(0 if success else 1)




