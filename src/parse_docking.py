import os
import csv
import glob
import logging
from typing import List, Tuple

from fetch_data import DATA_DIR

"""
Parsing utilities for docking outputs.

In a real setup, this module would parse log files from a docking engine
such as AutoDock Vina. For this demo pipeline we:

1. Look for Vina-like log files under data/docking/*/log.txt.
2. Parse the best binding score from each log (if present).
3. If no real logs are found, synthesize a simple score placeholder so the
   rest of the pipeline still runs end-to-end.

This keeps the code honest (no hard-coded winners) while avoiding a
hard dependency on an external docking installation.
"""

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

DOCKING_ROOT = os.path.join(DATA_DIR, "docking")
RESULTS_CSV = os.path.join(DATA_DIR, "docking_results.csv")


def find_variant_dirs(root: str) -> List[str]:
    """
    Return a list of variant directories under the docking root.
    """
    if not os.path.exists(root):
        return []
    return [p for p in glob.glob(os.path.join(root, "*")) if os.path.isdir(p)]


def parse_vina_log(path: str) -> float | None:
    """
    Very minimal parser for AutoDock Vina log files.

    We look for a line starting with 'REMARK VINA RESULT:' and parse the
    first numeric token as the binding score (kcal/mol).
    """
    if not os.path.exists(path):
        return None

    try:
        with open(path) as f:
            for line in f:
                if "REMARK VINA RESULT:" in line:
                    parts = line.split()
                    for token in parts:
                        try:
                            return float(token)
                        except ValueError:
                            continue
    except Exception:
        return None

    return None


def collect_docking_results() -> List[Tuple[str, float]]:
    """
    Collect best docking scores for each variant.

    If no real logs are present, we assign a neutral placeholder score (0.0)
    and let later stages (e.g., scoring.py) differentiate variants using
    sequence-derived properties.
    """
    variant_dirs = find_variant_dirs(DOCKING_ROOT)
    if not variant_dirs:
        logging.warning("No docking variant directories found under %s", DOCKING_ROOT)
        return []

    results: List[Tuple[str, float]] = []
    any_real_score = False

    for vdir in sorted(variant_dirs):
        variant_id = os.path.basename(vdir)
        log_path = os.path.join(vdir, "log.txt")
        score = parse_vina_log(log_path)
        if score is not None:
            any_real_score = True
        else:
            score = 0.0  # neutral placeholder; refined in scoring.py
        results.append((variant_id, score))

    if not any_real_score:
        logging.info(
            "No real docking logs detected. Using placeholder scores (0.0) for all variants. "
            "Subsequent scoring will rely on sequence properties."
        )

    return results


def write_results_csv(results: List[Tuple[str, float]], out_path: str) -> None:
    """
    Write (variant_id, score) pairs to CSV.
    """
    if not results:
        logging.warning("No docking results to write")
        return

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["variant_id", "docking_score"])
        for vid, score in results:
            writer.writerow([vid, score])

    logging.info("Docking results written to %s", out_path)


if __name__ == "__main__":
    os.makedirs(DATA_DIR, exist_ok=True)
    res = collect_docking_results()
    write_results_csv(res, RESULTS_CSV)


