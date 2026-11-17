"""Emit the dauer contactome (membrane-adjacency) edge list from Yim et al. 2024.

The "contactome" is the set of physical membrane contacts (appositions) between neurons — the
substrate that bounds the synaptic connectome (a synapse can only form where membranes touch, so
synapses are a subset of contacts). Yim et al. reconstructed BOTH the dauer synaptic connectome
(already in CIRCE as ``yim_2024_dauer``) and the dauer contactome from the SAME EM volume, so the
synapse-subset-of-contacts property is a clean within-dataset check (holds for ~99% of dauer
synapse pairs). This reads the dauer contact-area matrix and emits the committed
``dauer_contactome.csv`` the build ingests as ``adjacency`` connections.

Source (not vendored — large third-party file):
  Yim H, Walker DS, ... Cook SJ. "Comparative connectomics of dauer reveals developmental
  plasticity." Nat Commun 15:1546 (2024). PMID 38413604; doi:10.1038/s41467-024-45943-3.
  Supplementary Data 6 (``41467_2024_45943_MOESM9_ESM.xlsx``), sheet "Dauer": a symmetric
  neuron x neuron **contact-area** matrix (raw EM-measured area; 0 = no contact). The file also
  carries a "Dauer_normalized" sheet and other stages (L1-L3, Adult) -- we read the raw "Dauer"
  sheet only. The matrix header is cell names (row 0) / indices (row 1); column 0/1 are the row
  cell / index; data starts at row 2, column 2.

    python -I build_contactome.py 41467_2024_45943_MOESM9_ESM.xlsx   # or set YIM_CONTACTOME_XLSX
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parent


def xlsx_path() -> Path:
    if len(sys.argv) > 1:
        return Path(sys.argv[1]).expanduser()
    env = os.environ.get("YIM_CONTACTOME_XLSX")
    if env:
        return Path(env).expanduser()
    sys.exit(
        "usage: python -I build_contactome.py <...MOESM9_ESM.xlsx>  (or set YIM_CONTACTOME_XLSX)"
    )


def main() -> None:
    wb = openpyxl.load_workbook(xlsx_path(), read_only=True, data_only=True)
    rows = list(wb["Dauer"].iter_rows(values_only=True))
    wb.close()
    col_names = [str(c).strip() if c is not None else "" for c in rows[0][2:]]
    edges: list[tuple[str, str, float]] = []
    for r in rows[2:]:
        a = str(r[0]).strip() if r[0] is not None else ""
        if not a:
            continue
        for j, v in enumerate(r[2:]):
            b = col_names[j]
            if a < b and isinstance(v, (int, float)) and v > 0:  # upper triangle, nonzero
                edges.append((a, b, float(v)))
    edges.sort()
    with open(HERE / "dauer_contactome.csv", "w", newline="") as f:
        wri = csv.writer(f)
        wri.writerow(["cell_a", "cell_b", "contact_area"])
        for a, b, w in edges:
            wri.writerow([a, b, int(w) if float(w).is_integer() else w])
    cells = {c for a, b, _ in edges for c in (a, b)}
    print(f"wrote dauer_contactome.csv: {len(edges)} contacts over {len(cells)} dauer neurons")


if __name__ == "__main__":
    main()
