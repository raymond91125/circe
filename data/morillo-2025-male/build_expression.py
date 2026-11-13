"""Emit per-neuron neuropeptide-ligand expression from the Morillo et al. 2025 male/herm atlas.

Phase A of the Morillo 2025 ingestion (see ``docs/proposals/morillo-2025-male-expression.md``):
the sex-shared FLP/NLP/INS neuropeptide **ligands**, for both sexes, from supplement **mmc8**
("Herm/Male - Sexshared FlpNlpIns SCT"). mmc8 holds SCT-normalized mean expression per neuron
column x gene row, continuous (not binary). We binarize at a CeNGEN-calibrated cutoff and emit
the per-neuron ``cell,gene_symbol`` lists the build ingests (one per sex).

Binarization cutoff (calibrated, not arbitrary)
-----------------------------------------------
The Morillo-herm SCT value is thresholded at **T = 0.2526**, the cutoff that best reproduces
(max F1) CeNGEN threshold-4 hermaphrodite calls over the 42 neuropeptide-ligand genes shared with
``np_genes.csv`` and the mmc8 columns that map 1:1 to a CeNGEN neuron/class (AUC = 0.907,
precision 0.80 / recall 0.64; CeNGEN prevalence 12.2% vs Morillo@T 9.7%). The *same* T is applied
to the male sheet: male and herm mmc8 share one platform, pipeline, and normalization, so T is a
property of the assay, calibrated where a gold standard exists (herm) and transferred to the sex
where none does (male). See the proposal doc for the full validation.

Label handling (reconciliation to CIRCE cells)
----------------------------------------------
mmc8 columns are a mix of individual neurons (``HOA``), neuron classes (``AIY``, ``CEM``), and
clusters the snRNA could not resolve (``PHA/PHB/ASH``) or generic groupings (``RnA``,
``CP Neurons``, ``Male 1``). Policy for Phase A (per the proposal's resolved curation choices):

  * **Drop merged clusters** -- any label naming more than one class (``/`` or ``,``): the value is
    for the pooled group, not per neuron, so it is not ingested as a per-cell call.
  * **Drop non-neuronal / ungroundable groupings** (``Intestine``, ``Sperm``, ``RnA``,
    ``CP Neurons``, ``Male N``, ...): not individual CIRCE neurons.
  * **Alias ventral-cord classes** to CIRCE's ``-n`` class names (``VB`` -> ``VBn`` etc.).
  * Everything kept is an individual cell name or a neuron class; the build (assemble) expands a
    class label to its member cells present in that sex.

The gene symbol -> WBGene table (``morillo_genes.csv``) is a committed, pre-resolved curation
artifact (MyGene.info, species 6239; cross-checked against ``np_genes.csv`` with zero mismatches),
so this build step is offline and deterministic. nlp-73 and nlp-77 (not in MyGene under those
names) were resolved by WormBase curation to their canonical synonyms ctdh-1 / srlf-1.

    MORILLO_MMC8=/path/to/mmc8.xlsx python -I build_expression.py
    # or: python -I build_expression.py /path/to/mmc8.xlsx
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parent
GENES_CSV = HERE / "morillo_genes.csv"

#: Calibrated binarization cutoff on the mmc8 SCT value (see module docstring).
CUTOFF_T = 0.2526

#: mmc8 sheet -> (sex, output filename).
SHEETS = {
    "Herm - Sexshared FlpNlpIns SCT": ("hermaphrodite", "morillo_herm_expression.csv"),
    "Male - Sexshared FlpNlpIns SCT": ("male", "morillo_male_expression.csv"),
}

#: Morillo ventral-cord class label -> CIRCE ``-n`` class name.
ALIAS = {"VB": "VBn", "VC": "VCn", "VD": "VDn", "DA": "DAn", "DB": "DBn", "AS": "ASn"}

#: Non-neuronal and ungroundable grouping labels dropped outright (not individual CIRCE neurons).
DROP_LABELS = {
    # non-neuronal
    "Intestine",
    "Int.",
    "Sperm",
    "Spermatheca",
    "Hypodermis",
    "Hyp.2",
    "Pharynx",
    "Excretory Cell",
    "Ex. C.",
    "Ex. C",
    "Musc.",
    # functional super-groups / unresolved generic clusters
    "Mechanosensory Neurons",
    "Motor Neurons",
    "M neurons",
    "I Neurons",
    "I neurons",
    "CP Neurons",
    "RnA",
    "RnB",
    "IL2DV",
    "IL2LT",
    "PVS",
}


def keep(label: str) -> bool:
    """True if a label is a single groundable neuron/class (not merged/pooled/non-neuronal)."""
    if "/" in label or "," in label:
        return False  # merged cluster of >1 class
    if label.startswith("Male "):
        return False  # generic unidentified male cluster (Male 1..12, Male Inter 8)
    return label not in DROP_LABELS


def mmc8_path() -> Path:
    if len(sys.argv) > 1:
        return Path(sys.argv[1]).expanduser()
    env = os.environ.get("MORILLO_MMC8")
    if env:
        return Path(env).expanduser()
    sys.exit("usage: python -I build_expression.py <mmc8.xlsx>  (or set MORILLO_MMC8)")


def load_gene_symbols() -> set[str]:
    with open(GENES_CSV, newline="") as f:
        return {r["symbol"].strip() for r in csv.DictReader(f)}


def read_sheet(wb, sheet: str, symbols: set[str]) -> list[tuple[str, str]]:
    """(label, gene_symbol) for every kept column x gene whose SCT value exceeds the cutoff."""
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    labels = [str(c).strip() if c is not None else "" for c in next(it)[1:]]
    out: set[tuple[str, str]] = set()
    for row in it:
        if row[0] is None:
            continue
        sym = str(row[0]).strip()
        if sym not in symbols:
            continue
        for i, raw in enumerate(labels):
            if not raw or not keep(raw):
                continue
            v = row[i + 1]
            try:
                val = float(v) if v is not None else 0.0
            except (TypeError, ValueError):
                val = 0.0
            if val > CUTOFF_T:
                out.add((ALIAS.get(raw, raw), sym))
    return sorted(out)


def main() -> None:
    path = mmc8_path()
    symbols = load_gene_symbols()
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for sheet, (sex, fname) in SHEETS.items():
        rows = read_sheet(wb, sheet, symbols)
        with open(HERE / fname, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["cell", "gene_symbol"])
            w.writerows(rows)
        cells = sorted({c for c, _ in rows})
        genes = sorted({g for _, g in rows})
        print(f"wrote {fname}: {len(rows)} rows, {len(cells)} labels, {len(genes)} genes ({sex})")
    wb.close()


if __name__ == "__main__":
    main()
