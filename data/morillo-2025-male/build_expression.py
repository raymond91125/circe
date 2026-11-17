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
mmc8 columns are a mix of individual neurons (``HOA``), neuron classes (``AIY``, ``CEM``), clusters
the snRNA could not resolve (``PHA/PHB/ASH``), serial-homolog groupings (``RnA``, ``CP Neurons``,
``R5,7,9A``), and non-neuronal / generic clusters (``Intestine``, ``Male 1``). Policy:

  * **Keep serial-homolog groupings** -- labels that are one ontology class or serial homologs of
    one (rays ``RnA``/``RnB``, ``IL2DV``/``IL2LT``, ``CP Neurons``, ``VC4/5``, see
    ``morillo_groupings.csv``): kept even with ``/``/``,``; the build expands them to member cells
    with a distinct ``grouped`` confidence (likely-similar expression across homologs).
  * **Keep ``PVS``** -- an individual male-specific neuron (WBbt:0007837), a per-cell call.
  * **Drop merged clusters of DISTINCT types** (``PHA/PHB``, ``AVE/AVA``, ``CEM/Ray1,5,7B``): a
    pooled value over different classes, not groundable per neuron.
  * **Drop non-neuronal / heterogeneous super-groups** (``Intestine``, ``Sperm``, ``Motor
    Neurons``, ``Male N``, ...).
  * **Alias ventral-cord classes** to CIRCE's ``-n`` class names (``VB`` -> ``VBn`` etc.).
  * Everything else kept is an individual cell or a neuron class; the build expands a class label to
    its member cells present in that sex.

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
GROUPINGS_CSV = HERE / "morillo_groupings.csv"


def load_grouping_labels() -> set[str]:
    """Serial-homolog grouping labels to KEEP (expanded to member cells by the build, see
    morillo_groupings.csv). These are ontology-confirmed single classes (rays, IL2, CP, VC), so
    unlike merged clusters of distinct types they are not dropped."""
    with open(GROUPINGS_CSV, newline="") as f:
        return {r["label"].strip() for r in csv.DictReader(f)}


GROUPINGS = load_grouping_labels()

#: Calibrated binarization cutoff on the mmc8 SCT value (see module docstring).
CUTOFF_T = 0.2526

#: mmc8 sheet -> (sex, output filename).
SHEETS = {
    "Herm - Sexshared FlpNlpIns SCT": ("hermaphrodite", "morillo_herm_expression.csv"),
    "Male - Sexshared FlpNlpIns SCT": ("male", "morillo_male_expression.csv"),
}

#: Morillo ventral-cord class label -> CIRCE ``-n`` class name.
ALIAS = {"VB": "VBn", "VC": "VCn", "VD": "VDn", "DA": "DAn", "DB": "DBn", "AS": "ASn"}

#: Non-neuronal and heterogeneous super-group labels dropped outright (not a single groundable
#: neuron/class). Serial-homolog groupings (rays, IL2, CP, VC) are NOT here -- they are kept and
#: expanded (see GROUPINGS / morillo_groupings.csv). PVS is NOT here either -- it is an individual
#: male-specific neuron (WBbt:0007837), kept as a per-cell call.
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
    # heterogeneous functional super-groups (span many distinct classes, not one serial-homolog set)
    "Mechanosensory Neurons",
    "Motor Neurons",
    "M neurons",
    "I Neurons",
    "I neurons",
    "Touch Neurons",
}


def keep(label: str) -> bool:
    """True if a label is a single groundable neuron/class OR a serial-homolog grouping.

    Serial-homolog groupings (ontology-confirmed one class; rays/IL2/CP/VC) are kept even when they
    contain ``/`` or ``,`` (e.g. ``R5,7,9A``), to be expanded to member cells downstream. Merged
    clusters of *distinct* types (``PHA/PHB``, ``AVE/AVA``) are still dropped.
    """
    if label in GROUPINGS:
        return True
    if "/" in label or "," in label:
        return False  # merged cluster of >1 distinct class
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
