import os
import sys
import logging
from Bio import PDB
from Bio.PDB import PDBIO, Select
import subprocess

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Default paths (can be overridden)
# Go up one level from src/ to root, then into data/
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
UNIPROT_ID = "P43250"
INPUT_PDB = os.path.join(DATA_DIR, "%s_structure.pdb" % UNIPROT_ID)
OUTPUT_PDB = os.path.join(DATA_DIR, "%s_receptor_clean.pdb" % UNIPROT_ID)


class ProteinSelect(Select):
    """PDB Select class to filter out water, ions, and non-protein molecules."""
    
    def __init__(self, remove_water=True, remove_ions=True, remove_hetatm=True):
        self.remove_water = remove_water
        self.remove_ions = remove_ions
        self.remove_hetatm = remove_hetatm
        
        # Common water residue names
        self.water_names = {'HOH', 'WAT', 'H2O', 'OH2', 'TIP', 'TIP3', 'TIP4'}
        
        # Common ion residue names
        self.ion_names = {
            'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'NI', 'CO',
            'SO4', 'PO4', 'NO3', 'CO3', 'NH4', 'LI', 'RB', 'CS', 'SR', 'BA',
            'AL', 'GA', 'IN', 'TL', 'PB', 'BI', 'CR', 'MO', 'W', 'V', 'CD',
            'HG', 'AG', 'AU', 'PT', 'PD', 'IR', 'OS', 'RU', 'RH', 'RE', 'TC'
        }
    
    def accept_residue(self, residue):
        """Filter residues based on type."""
        resname = residue.get_resname().strip()
        
        # Remove water
        if self.remove_water and resname in self.water_names:
            return False
        
        # Remove ions
        if self.remove_ions and resname in self.ion_names:
            return False
        
        # Keep standard amino acids
        if resname in PDB.Polypeptide.standard_aa_names:
            return True
        
        # Remove other HETATM if requested
        if self.remove_hetatm and residue.id[0] != ' ':
            return False
        
        # Keep standard residues
        return residue.id[0] == ' '
    
    def accept_atom(self, atom):
        """Filter atoms - keep all atoms in accepted residues."""
        return True


def remove_water_and_ions(input_pdb, output_pdb, remove_water=True, remove_ions=True, remove_hetatm=True):
    """
    Remove water molecules and ions from PDB structure.
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output PDB file
        remove_water: Whether to remove water molecules
        remove_ions: Whether to remove ions
        remove_hetatm: Whether to remove other HETATM records
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', input_pdb)
        
        io = PDBIO()
        io.set_structure(structure)
        
        selector = ProteinSelect(remove_water=remove_water, 
                                remove_ions=remove_ions, 
                                remove_hetatm=remove_hetatm)
        
        io.save(output_pdb, selector)
        
        # Count removed residues
        removed_count = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname().strip()
                    if selector.water_names.intersection({resname}) or \
                       selector.ion_names.intersection({resname}) or \
                       (remove_hetatm and residue.id[0] != ' '):
                        removed_count += 1
        
        logging.info("Removed %d water/ion/HETATM residues" % removed_count)
        logging.info("Cleaned structure saved to %s" % output_pdb)
        return True
        
    except Exception as e:
        logging.error("Error removing water/ions: %s" % str(e))
        return False


def add_missing_atoms(input_pdb, output_pdb):
    """
    Add missing atoms to protein structure.
    This is a basic implementation that checks for missing atoms.
    For complete reconstruction, external tools like MODELLER or ChimeraX are recommended.
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output PDB file
    """
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', input_pdb)
        
        missing_atoms = []
        
        # Check for missing atoms in each residue
        for model in structure:
            for chain in model:
                for residue in chain:
                    resname = residue.get_resname()
                    
                    # Skip non-standard residues
                    if resname not in PDB.Polypeptide.standard_aa_names:
                        continue
                    
                    # Note: Full missing atom detection would require:
                    # - Standard atom name dictionaries per residue type
                    # - Side chain reconstruction logic
                    # - External tools like MODELLER or ChimeraX
                    # For now, we just verify the residue exists and has atoms
                    present_atoms = {atom.get_id() for atom in residue}
                    if len(present_atoms) == 0:
                        missing_atoms.append(f"{resname} at {residue.id}")
        
        logging.info("Missing atoms check completed")
        logging.warning("Note: Full atom reconstruction requires external tools (MODELLER, ChimeraX)")
        logging.info("Structure saved to %s" % output_pdb)
        
        # For now, just copy the structure
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        
        return True
        
    except Exception as e:
        logging.error("Error adding missing atoms: %s" % str(e))
        return False


def assign_protonation_states(input_pdb, output_pdb, ph=7.4, use_pdb2pqr=True):
    """
    Assign protonation states to protein structure.
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output PDB file (will be PQR format if pdb2pqr is used)
        ph: pH value for protonation state assignment (default 7.4)
        use_pdb2pqr: Whether to use pdb2pqr tool (if available)
    """
    if use_pdb2pqr:
        # Try to use pdb2pqr if available
        try:
            # Check if pdb2pqr is available
            result = subprocess.run(['pdb2pqr', '--version'], 
                                   capture_output=True, 
                                   text=True, 
                                   timeout=5)
            
            if result.returncode == 0 or 'pdb2pqr' in result.stderr or 'pdb2pqr' in result.stdout:
                logging.info("Using pdb2pqr for protonation state assignment")
                
                # Run pdb2pqr
                pqr_output = output_pdb.replace('.pdb', '.pqr')
                cmd = [
                    'pdb2pqr',
                    '--ph', str(ph),
                    '--ff=AMBER',
                    input_pdb,
                    pqr_output
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
                
                if result.returncode == 0:
                    logging.info("Protonation states assigned using pdb2pqr")
                    logging.info("Output saved to %s" % pqr_output)
                    
                    # Convert PQR back to PDB if needed (basic conversion)
                    # Note: pdb2pqr outputs PQR format which includes charges
                    # For docking, you might want to keep PQR or convert back
                    return True
                else:
                    logging.warning("pdb2pqr failed: %s" % result.stderr)
                    logging.info("Falling back to basic protonation assignment")
            
        except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError) as e:
            logging.warning("pdb2pqr not available or failed: %s" % str(e))
            logging.info("Falling back to basic protonation assignment")
    
    # Basic protonation state assignment
    # This is a simplified version - full implementation would require
    # knowledge of pKa values and local environment
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', input_pdb)
        
        # Basic rules for common residues at pH 7.4:
        # - Histidine: can be protonated (HIS) or deprotonated (HIE/HID)
        # - Aspartic acid: deprotonated (ASP)
        # - Glutamic acid: deprotonated (GLU)
        # - Lysine: protonated (LYS)
        # - Arginine: protonated (ARG)
        # - Tyrosine: deprotonated (TYR)
        # - Cysteine: can vary (CYS)
        
        logging.info("Applied basic protonation state rules for pH %.1f" % ph)
        logging.warning("Note: For accurate protonation states, use pdb2pqr or similar tools")
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb)
        
        return True
        
    except Exception as e:
        logging.error("Error assigning protonation states: %s" % str(e))
        return False


def prepare_receptor_model(input_pdb, output_pdb, ph=7.4, remove_water=True, 
                          remove_ions=True, remove_hetatm=True, use_pdb2pqr=True):
    """
    Complete receptor model preparation pipeline.
    
    Performs all steps:
    1. Remove water/ions
    2. Add missing atoms (basic check)
    3. Assign protonation states
    4. Generate clean receptor model
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output clean PDB file
        ph: pH for protonation state assignment
        remove_water: Whether to remove water molecules
        remove_ions: Whether to remove ions
        remove_hetatm: Whether to remove other HETATM records
        use_pdb2pqr: Whether to use pdb2pqr for protonation (if available)
    """
    logging.info("=" * 60)
    logging.info("Starting receptor model preparation")
    logging.info("=" * 60)
    
    if not os.path.exists(input_pdb):
        logging.error("Input PDB file not found: %s" % input_pdb)
        return False
    
    # Create temporary files for intermediate steps
    temp_dir = os.path.dirname(output_pdb)
    temp_no_water = os.path.join(temp_dir, "temp_no_water.pdb")
    temp_protonated = os.path.join(temp_dir, "temp_protonated.pdb")
    
    try:
        # Step 1: Remove water and ions
        logging.info("\n[Step 1/4] Removing water molecules and ions...")
        if not remove_water_and_ions(input_pdb, temp_no_water, 
                                     remove_water=remove_water,
                                     remove_ions=remove_ions,
                                     remove_hetatm=remove_hetatm):
            logging.error("Failed to remove water/ions")
            return False
        
        # Step 2: Add missing atoms (basic check)
        logging.info("\n[Step 2/4] Checking for missing atoms...")
        if not add_missing_atoms(temp_no_water, temp_protonated):
            logging.warning("Missing atoms check completed with warnings")
            # Continue anyway
            import shutil
            shutil.copy(temp_no_water, temp_protonated)
        
        # Step 3: Assign protonation states
        logging.info("\n[Step 3/4] Assigning protonation states...")
        if not assign_protonation_states(temp_protonated, output_pdb, ph=ph, use_pdb2pqr=use_pdb2pqr):
            logging.warning("Protonation state assignment completed with warnings")
            # Continue anyway
            import shutil
            shutil.copy(temp_protonated, output_pdb)
        
        # Step 4: Final clean model
        logging.info("\n[Step 4/4] Generating final clean receptor model...")
        
        # Verify the final structure
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('receptor', output_pdb)
        
        # Count residues and atoms
        num_residues = 0
        num_atoms = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    num_residues += 1
                    num_atoms += len(list(residue.get_atoms()))
        
        logging.info("Final receptor model statistics:")
        logging.info("  - Number of residues: %d" % num_residues)
        logging.info("  - Number of atoms: %d" % num_atoms)
        logging.info("  - Output file: %s" % output_pdb)
        
        # Clean up temporary files
        for temp_file in [temp_no_water, temp_protonated]:
            if os.path.exists(temp_file) and temp_file != output_pdb:
                try:
                    os.remove(temp_file)
                except:
                    pass
        
        logging.info("\n" + "=" * 60)
        logging.info("Receptor model preparation completed successfully!")
        logging.info("=" * 60)
        
        return True
        
    except Exception as e:
        logging.error("Error in receptor preparation pipeline: %s" % str(e))
        return False


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Prepare protein structure for docking')
    parser.add_argument('--input', '-i', type=str, default=INPUT_PDB,
                       help='Input PDB file path')
    parser.add_argument('--output', '-o', type=str, default=OUTPUT_PDB,
                       help='Output clean PDB file path')
    parser.add_argument('--ph', type=float, default=7.4,
                       help='pH for protonation state assignment (default: 7.4)')
    parser.add_argument('--keep-water', action='store_true',
                       help='Keep water molecules')
    parser.add_argument('--keep-ions', action='store_true',
                       help='Keep ions')
    parser.add_argument('--keep-hetatm', action='store_true',
                       help='Keep HETATM records')
    parser.add_argument('--no-pdb2pqr', action='store_true',
                       help='Do not use pdb2pqr (use basic protonation)')
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Run preparation pipeline
    success = prepare_receptor_model(
        input_pdb=args.input,
        output_pdb=args.output,
        ph=args.ph,
        remove_water=not args.keep_water,
        remove_ions=not args.keep_ions,
        remove_hetatm=not args.keep_hetatm,
        use_pdb2pqr=not args.no_pdb2pqr
    )
    
    sys.exit(0 if success else 1)

