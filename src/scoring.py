import os
import csv
import logging
from typing import Dict, Tuple, List

from fetch_data import DATA_DIR
from design_library import compute_properties

"""
Scoring utilities for peptide variants.

This module combines:
1. Docking scores (if available) from data/docking_results.csv.
2. Simple sequence-derived properties (charge, hydrophobicity, length).

into a single composite score suitable for ranking variants. The scoring
formula is intentionally simple and transparent so that employers can
see how to swap in their own weighting or more sophisticated models.
"""

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

LIBRARY_FASTA = os.path.join(DATA_DIR, "library.fasta")
DOCKING_RESULTS = os.path.join(DATA_DIR, "docking_results.csv")
SCORED_CSV = os.path.join(DATA_DIR, "scored_variants.csv")


def read_library_fasta(path: str) -> Dict[str, str]:
    """
    Read multi-FASTA into a dict: {variant_id: sequence}.
    
    Variant IDs are extracted as the canonical ID (before the first '|')
    to match the format used in docking_results.csv.
    """
    sequences: Dict[str, str] = {}
    current_id = None
    seq_chunks: List[str] = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(seq_chunks)
                raw_id = line[1:].split()[0]
                # Extract canonical ID (before first '|') to match docking_results.csv
                current_id = raw_id.split("|")[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

    if current_id is not None:
        sequences[current_id] = "".join(seq_chunks)

    return sequences


def read_docking_results(path: str) -> Dict[str, float]:
    """
    Read docking_results.csv into a dict: {variant_id: docking_score}.
    Missing file -> empty dict.
    """
    if not os.path.exists(path):
        logging.warning("Docking results file not found at %s", path)
        return {}

    scores: Dict[str, float] = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                scores[row["variant_id"]] = float(row["docking_score"])
            except (KeyError, ValueError):
                continue
    return scores


def compute_composite_scores() -> List[Dict[str, str]]:
    """
    Combine docking scores and sequence properties into a composite score.

    Scoring heuristic (lower is better):
        composite = docking_score
                    + 0.1 * abs(net_charge)
                    + 0.2 * abs(avg_hydrophobicity)

    If docking_score is missing, we treat it as 0.0 (neutral) so that
    sequence properties still drive differentiation.
    """
    if not os.path.exists(LIBRARY_FASTA):
        raise FileNotFoundError(
            f"Library FASTA not found at {LIBRARY_FASTA}. "
            "Run design_library.py first."
        )

    sequences = read_library_fasta(LIBRARY_FASTA)
    docking_scores = read_docking_results(DOCKING_RESULTS)

    records: List[Dict[str, str]] = []
    for vid, seq in sequences.items():
        props = compute_properties(seq)
        docking_score = docking_scores.get(vid, 0.0)
        composite = docking_score + 0.1 * abs(props.net_charge) + 0.2 * abs(props.avg_hydrophobicity)

        records.append(
            {
                "variant_id": vid,
                "sequence": seq,
                "docking_score": f"{docking_score:.3f}",
                "net_charge": f"{props.net_charge:.3f}",
                "avg_hydrophobicity": f"{props.avg_hydrophobicity:.3f}",
                "length": str(props.length),
                "composite_score": f"{composite:.3f}",
            }
        )

    # Sort by composite score ascending (lower = better)
    records.sort(key=lambda r: float(r["composite_score"]))
    return records


def write_scored_csv(records: List[Dict[str, str]], out_path: str) -> None:
    if not records:
        logging.warning("No records to score")
        return

    fieldnames = [
        "rank",
        "variant_id",
        "sequence",
        "docking_score",
        "net_charge",
        "avg_hydrophobicity",
        "length",
        "composite_score",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, rec in enumerate(records, start=1):
            row = {"rank": i}
            row.update(rec)
            writer.writerow(row)

    logging.info("Scored variants written to %s", out_path)


if __name__ == "__main__":
    os.makedirs(DATA_DIR, exist_ok=True)
    recs = compute_composite_scores()
    write_scored_csv(recs, SCORED_CSV)


