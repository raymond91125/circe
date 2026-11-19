"""Project CeNGEN threshold-4 expression for the Fenyves ionotropic receptor genes onto
CIRCE's individual neurons (template: ``../cengen-monoamine-ionotropic/build_expression.py``).

Reads:
  * ``ionotropic_receptor_genes.csv`` (this dir; the committed gene table, 64 rows / 62 WBGene).
  * the CeNGEN **published** threshold-4 matrix ``021821_stringent_threshold4.csv`` (gene x 128
    neuron CLASSES, keyed by ``Wormbase_ID``; value > 0 = expressed). NOT vendored (9.8 MB);
    give its path via env ``CENGEN_T4`` (sha256 in ``MANIFEST.md``).
  * ``../neuron-graph/neurons.json`` (the CIRCE cell universe + per-cell class token).

Writes (next to this file):
  * ``cengen_class_map.csv`` -- every CeNGEN class -> its CIRCE member cells (auditable).
  * ``ionotropic_receptor_expression.csv`` -- ``cell,gene_symbol,wbgene`` for expressed only.

CeNGEN's 128 classes collapse neurons: 105 match a CIRCE class token, 4 a cell name, and the 19
aggregates below are mapped explicitly. Two collapses are irreducibly ambiguous and resolved
**conservatively (union)**, noted in MANIFEST:
  * ``VD_DD`` -> all 13 VD + 6 DD (CeNGEN cannot separate VD from DD motor neurons).
  * ``AWC_ON`` / ``AWC_OFF`` -> both AWCL and AWCR (ON/OFF identity is stochastic, not L/R-fixed).

Usage:  CENGEN_T4=/path/to/021821_stringent_threshold4.csv python -I build_expression.py
"""

from __future__ import annotations

import csv
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
GENE_CSV = HERE / "ionotropic_receptor_genes.csv"
NEURONS_JSON = HERE.parent / "neuron-graph" / "neurons.json"

# CeNGEN aggregate class -> explicit CIRCE cell members (the 19 that match no token/name).
SPECIAL: dict[str, list[str]] = {
    "AS": [f"AS{i}" for i in range(1, 12)],  # AS1..AS11
    "DA": [f"DA{i}" for i in range(1, 9)],  # DA1..DA8 (DA9 is its own CeNGEN class)
    "DB": [f"DB{i}" for i in range(2, 8)],  # DB2..DB7 (DB01 separate)
    "DB01": ["DB1"],
    "VA": [f"VA{i}" for i in range(1, 12)],  # VA1..VA11 (VA12 separate)
    "VB": [f"VB{i}" for i in range(3, 12)],  # VB3..VB11 (VB01/VB02 separate)
    "VB01": ["VB1"],
    "VB02": ["VB2"],
    "VC": ["VC1", "VC2", "VC3", "VC6"],  # VC_4_5 separate
    "VC_4_5": ["VC4", "VC5"],
    "VD_DD": [f"VD{i}" for i in range(1, 14)] + [f"DD{i}" for i in range(1, 7)],
    "IL2_DV": ["IL2DL", "IL2DR", "IL2VL", "IL2VR"],
    "IL2_LR": ["IL2L", "IL2R"],
    "RMD_DV": ["RMDDL", "RMDDR", "RMDVL", "RMDVR"],
    "RMD_LR": ["RMDL", "RMDR"],
    "RME_DV": ["RMED", "RMEV"],
    "RME_LR": ["RMEL", "RMER"],
    "AWC_ON": ["AWCL", "AWCR"],
    "AWC_OFF": ["AWCL", "AWCR"],
}


def _pos(v: str) -> bool:
    try:
        return float(v) > 0
    except (TypeError, ValueError):
        return False


def main() -> None:
    cengen_csv = os.environ.get("CENGEN_T4")
    if not cengen_csv:
        sys.exit("set CENGEN_T4 to the path of 021821_stringent_threshold4.csv (not vendored)")

    cells = json.loads(NEURONS_JSON.read_text())
    names = {c["name"] for c in cells}
    by_class: dict[str, list[str]] = defaultdict(list)
    for c in cells:
        by_class[c["classes"]].append(c["name"])

    def resolve(cengen_class: str) -> tuple[str, list[str]]:
        if cengen_class in SPECIAL:
            kind, members = "special", SPECIAL[cengen_class]
        elif cengen_class in by_class:
            kind, members = "class-token", sorted(by_class[cengen_class])
        elif cengen_class in names:
            kind, members = "cell-name", [cengen_class]
        else:
            return "UNRESOLVED", []
        return kind, [m for m in members if m in names]

    gene_rows = list(csv.DictReader(GENE_CSV.open(newline="")))
    wb_of = {r["symbol"].strip(): r["wbgene"].strip().removeprefix("WB:") for r in gene_rows}
    wanted = set(wb_of.values())

    with open(cengen_csv, newline="") as f:
        rd = csv.reader(f)
        classes = [h.strip().strip('"') for h in next(rd)[3:]]
        expr_by_wb: dict[str, set[str]] = {}
        for row in rd:
            wb = row[2].strip().strip('"')
            if wb in wanted:
                expr_by_wb[wb] = {classes[i] for i, v in enumerate(row[3:]) if _pos(v)}

    unresolved = []
    with open(HERE / "cengen_class_map.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cengen_class", "resolution", "n_cells", "cells"])
        for cl in classes:
            kind, members = resolve(cl)
            if kind == "UNRESOLVED":
                unresolved.append(cl)
            w.writerow([cl, kind, len(members), ";".join(members)])

    rows: set[tuple[str, str, str]] = set()
    missing = []
    for sym, wb in wb_of.items():
        if wb not in expr_by_wb:
            missing.append(sym)
            continue
        for cl in expr_by_wb[wb]:
            for cell in resolve(cl)[1]:
                rows.add((cell, sym, f"WB:{wb}"))

    with open(HERE / "ionotropic_receptor_expression.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cell", "gene_symbol", "wbgene"])
        w.writerows(sorted(rows))

    print(f"classes={len(classes)} unresolved={unresolved}")
    print(
        f"genes with neuronal expression: {len(expr_by_wb)}/{len(wanted)} "
        f"(not in CeNGEN neuronal set: {sorted(set(missing))})"
    )
    print(f"expression rows={len(rows)} cells={len({r[0] for r in rows})}")


if __name__ == "__main__":
    main()
