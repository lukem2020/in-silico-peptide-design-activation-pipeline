import os
import csv
import logging
from typing import List, Tuple

from design_library import load_first_fasta_sequence
from fetch_data import DATA_DIR, UNIPROT_ID


"""
Docking preparation utilities.

For this demo pipeline we:
1. Take the prepared receptor structure (from prepare_structures.py).
2. Read the designed peptide library FASTA.
3. Create a simple directory structure for each variant.
4. Write small metadata/config files that a real docking engine
   (e.g., AutoDock Vina) could consume.

We intentionally avoid calling heavy external docking tools here so the
repository remains self-contained. An employer can see clearly where
to plug in their preferred docking backend.
"""


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

RECEPTOR_PDB = os.path.join(DATA_DIR, f"{UNIPROT_ID}_receptor_clean.pdb")
LIBRARY_FASTA = os.path.join(DATA_DIR, "library.fasta")
DOCKING_ROOT = os.path.join(DATA_DIR, "docking")


def read_library_fasta(path: str) -> List[Tuple[str, str]]:
    """
    Read a multi-FASTA file and return (canonical_id, sequence) pairs.

    Canonical IDs are derived from the header up to the first '|' character,
    which keeps them filesystem-friendly while ignoring annotation fields.
    """
    records: List[Tuple[str, str]] = []
    current_id = None
    seq_chunks: List[str] = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(seq_chunks)))
                raw_id = line[1:].split()[0]
                current_id = raw_id.split("|")[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if current_id is not None:
        records.append((current_id, "".join(seq_chunks)))

    return records


def prepare_docking_inputs() -> None:
    """
    Create docking directories and simple config files for each variant.

    Directory layout (under data/docking/):
        data/docking/
            GRK6_variant_0001/
                ligand.fasta
                docking_meta.csv
            GRK6_variant_0002/
                ...

    The receptor PDB is not copied per variant; we reference the shared path.
    """
    if not os.path.exists(RECEPTOR_PDB):
        raise FileNotFoundError(
            f"Receptor PDB not found at {RECEPTOR_PDB}. "
            "Run prepare_structures.py first to generate it."
        )

    if not os.path.exists(LIBRARY_FASTA):
        raise FileNotFoundError(
            f"Library FASTA not found at {LIBRARY_FASTA}. "
            "Run design_library.py first to generate it."
        )

    os.makedirs(DOCKING_ROOT, exist_ok=True)

    variants = read_library_fasta(LIBRARY_FASTA)
    logging.info("Preparing docking inputs for %d variants", len(variants))

    for variant_id, seq in variants:
        variant_dir = os.path.join(DOCKING_ROOT, variant_id)
        os.makedirs(variant_dir, exist_ok=True)

        # Save ligand sequence (placeholder – a real workflow would include 3D structures)
        ligand_fasta_path = os.path.join(variant_dir, "ligand.fasta")
        with open(ligand_fasta_path, "w") as f:
            f.write(f">{variant_id}\n{seq}\n")

        # Minimal metadata/config that an external docking script could use
        meta_path = os.path.join(variant_dir, "docking_meta.csv")
        with open(meta_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["variant_id", "receptor_pdb", "ligand_fasta"])
            writer.writerow([variant_id, RECEPTOR_PDB, ligand_fasta_path])

    logging.info("Docking input preparation complete. Root directory: %s", DOCKING_ROOT)


if __name__ == "__main__":
    prepare_docking_inputs()


