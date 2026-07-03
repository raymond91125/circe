#!/usr/bin/env python
"""Parse the Wang et al. neurotransmitter atlas (Table S2) into per-neuron assignments.

Source: TableS2_expression.xlsx — supplement of Wang, Vidal, Sural et al., "A neurotransmitter
atlas of C. elegans males and hermaphrodites", eLife 2024 (doi:10.7554/eLife.95402). See
MANIFEST.md. The sheet has per-neuron reporter-expression columns (eat-4, unc-17, unc-25,
unc-47, cat-1, tph-1, cat-2, …) plus a consolidated "Neurotransmitter(s)" column; we take the
consolidated call (the reporter columns are color-coded and not machine-readable via values).

Run:  uv run --with 'openpyxl' python parse_atlas.py
Output: assignments.csv (class, neuron, neurotransmitters) next to this script.
"""

from __future__ import annotations

import csv
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parent
XLSX = HERE / "TableS2_expression.xlsx"
HEADER_ROW = 4  # 1-based; columns: Class | Neuron | Lineage | reporters… | Neurotransmitter(s)
COL_CLASS, COL_NEURON, COL_NT = 1, 2, 20  # 0-based indices


def parse(path: Path = XLSX) -> list[dict]:
    ws = openpyxl.load_workbook(path, read_only=True, data_only=True).worksheets[0]
    rows = list(ws.iter_rows(values_only=True))
    out: list[dict] = []
    current_class = ""
    for r in rows[HEADER_ROW:]:
        neuron = str(r[COL_NEURON]).strip() if len(r) > COL_NEURON and r[COL_NEURON] else ""
        if not neuron:
            continue
        cls = str(r[COL_CLASS]).strip() if len(r) > COL_CLASS and r[COL_CLASS] else ""
        current_class = cls or current_class  # class cell is merged over its members
        nt = str(r[COL_NT]).strip() if len(r) > COL_NT and r[COL_NT] else ""
        if nt:
            out.append({"class": current_class, "neuron": neuron, "neurotransmitters": nt})
    return out


def main() -> None:
    records = parse()
    with (HERE / "assignments.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["class", "neuron", "neurotransmitters"])
        w.writeheader()
        w.writerows(records)
    print(f"parsed {len(records)} neuron assignments -> assignments.csv")


if __name__ == "__main__":
    main()
