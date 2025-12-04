import os
import csv
import logging
from typing import List, Dict

from fetch_data import DATA_DIR

"""
Selection utilities for choosing peptide variants for synthesis/testing.

This module:
1. Reads scored_variants.csv produced by scoring.py.
2. Selects the top-N variants by composite score.
3. Writes:
   - selected_variants.csv with full scoring info for the chosen set.
   - selected_variants.fasta for downstream experimental workflows.
"""

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

SCORED_CSV = os.path.join(DATA_DIR, "scored_variants.csv")
SELECTED_CSV = os.path.join(DATA_DIR, "selected_variants.csv")
SELECTED_FASTA = os.path.join(DATA_DIR, "selected_variants.fasta")


def read_scored_variants(path: str) -> List[Dict[str, str]]:
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Scored variants file not found at {path}. "
            "Run scoring.py first."
        )

    records: List[Dict[str, str]] = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append(row)
    return records


def select_top_variants(records: List[Dict[str, str]], top_n: int) -> List[Dict[str, str]]:
    if not records:
        return []

    # Assume records are already ranked by composite_score (as in scoring.py)
    return records[:top_n]


def write_selected_csv(records: List[Dict[str, str]], out_path: str) -> None:
    if not records:
        logging.warning("No variants selected; nothing to write to %s", out_path)
        return

    fieldnames = list(records[0].keys())
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)

    logging.info("Selected variants written to %s", out_path)


def write_selected_fasta(records: List[Dict[str, str]], out_path: str) -> None:
    if not records:
        logging.warning("No variants selected; nothing to write to %s", out_path)
        return

    lines: List[str] = []
    for rec in records:
        vid = rec["variant_id"]
        seq = rec["sequence"]
        header = f">{vid}|composite={rec['composite_score']}|dock={rec['docking_score']}"
        lines.append(header)
        lines.append(seq)

    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    logging.info("Selected variants FASTA written to %s", out_path)


def main(top_n: int = 10) -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    records = read_scored_variants(SCORED_CSV)
    selected = select_top_variants(records, top_n=top_n)
    write_selected_csv(selected, SELECTED_CSV)
    write_selected_fasta(selected, SELECTED_FASTA)
    logging.info("Total selected variants: %d", len(selected))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Select top peptide variants for synthesis/testing.")
    parser.add_argument("--top", type=int, default=10, help="Number of top variants to select (default: 10)")
    args = parser.parse_args()

    main(top_n=args.top)


