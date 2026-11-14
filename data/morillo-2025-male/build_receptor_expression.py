"""Emit per-neuron neuropeptide-receptor (NPP-GPCR) expression from the Morillo 2025 atlas (mmc12).

Phase B of the Morillo 2025 ingestion (see ``docs/proposals/morillo-2025-male-expression.md``): the
**neuropeptide-receptor GPCRs** for both sexes, from supplement **mmc12** (the genome-wide
single-nucleus atlas), sheets ``HermThr1(0.5%)`` / ``MaleThr1(0.5%)`` -- SCT-normalized mean
expression per **cluster** column x gene row. This is the receptor side the connectome needs (the
ligands are already ingested per-neuron from mmc8 in Phase A, ``build_expression.py``); together
they let a male neuropeptide connectome be derived in Phase C, paralleling the hermaphrodite Ripoll
layer.

Binarization cutoff (calibrated, receptor-specific)
---------------------------------------------------
Receptors are lower-expressed and noisier than ligands, so they get their **own** cutoff. The
Morillo-herm SCT value is thresholded at **T = 0.0491**, the cutoff that best reproduces (max F1)
CeNGEN threshold-4 hermaphrodite receptor calls over the 45 NPP-GPCRs and the mmc12 clusters that
map 1:1 to a CeNGEN neuron/class (AUC = 0.834; precision 0.47 / recall 0.66; CeNGEN prevalence
16.0% vs Morillo@T 22.8%). The *same* T is transferred to the male sheet (shared platform /
pipeline / normalization). The modest precision reflects snRNA detecting low-level receptor
transcripts below CeNGEN's stringent threshold-4 -- a caveat for any edges derived from this layer
(Phase C), not for the expression records themselves. See the proposal doc.

Resolution & labels
-------------------
mmc12 is **cluster-level** (145 clusters; row 1 = cluster number, row 2 = annotation), coarser than
mmc8's per-neuron columns. The same reconciliation policy as Phase A applies: merged clusters
(``/`` or ``,``) and non-neuronal / ungroundable grouping columns are dropped; ventral-cord classes
are aliased to CIRCE's ``-n`` names; kept labels are an individual cell or a neuron class, expanded
to member cells by the build. Consequence: neurons that appear only inside a merged cluster get no
receptor call from this layer.

Six NPP-GPCRs of the Beets/Ripoll pair set (npr-37/39/40/41/42/43) fall below mmc12's 0.5% gene-
inclusion threshold and are absent from the atlas, so they are not in ``morillo_receptor_genes.csv``.

    MORILLO_MMC12=/path/to/mmc12.xlsx python -I build_receptor_expression.py
    # or: python -I build_receptor_expression.py /path/to/mmc12.xlsx
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import openpyxl

HERE = Path(__file__).resolve().parent
GENES_CSV = HERE / "morillo_receptor_genes.csv"

#: Calibrated receptor binarization cutoff on the mmc12 SCT value (see module docstring).
CUTOFF_T = 0.0491

#: mmc12 sheet -> (sex, output filename).
SHEETS = {
    "HermThr1(0.5%)": ("hermaphrodite", "morillo_herm_receptor_expression.csv"),
    "MaleThr1(0.5%)": ("male", "morillo_male_receptor_expression.csv"),
}

#: Morillo ventral-cord class label -> CIRCE ``-n`` class name.
ALIAS = {"VB": "VBn", "VC": "VCn", "VD": "VDn", "DA": "DAn", "DB": "DBn", "AS": "ASn"}

#: Non-neuronal and ungroundable grouping labels dropped outright (matched case-insensitively).
DROP_LABELS = {
    s.lower()
    for s in (
        "Intestine",
        "Int.",
        "Sperm",
        "Spermatheca",
        "Hypodermis",
        "Hyp.2",
        "Pharynx",
        "Excretory Cell",
        "Excretory cell",
        "Ex. C.",
        "Ex. C",
        "Musc.",
        "Mechanosensory Neurons",
        "Motor Neurons",
        "M neurons",
        "I Neurons",
        "I neurons",
        "Touch Neurons",
        "CP Neurons",
        "RnA",
        "RnB",
        "IL2DV",
        "IL2LT",
        "PVS",
    )
}


def keep(label: str) -> bool:
    """True if a label is a single groundable neuron/class (not merged/pooled/non-neuronal)."""
    if "/" in label or "," in label:
        return False
    if label.startswith("Male "):
        return False
    return label.lower() not in DROP_LABELS


def mmc12_path() -> Path:
    if len(sys.argv) > 1:
        return Path(sys.argv[1]).expanduser()
    env = os.environ.get("MORILLO_MMC12")
    if env:
        return Path(env).expanduser()
    sys.exit("usage: python -I build_receptor_expression.py <mmc12.xlsx>  (or set MORILLO_MMC12)")


def load_gene_symbols() -> set[str]:
    with open(GENES_CSV, newline="") as f:
        return {r["symbol"].strip() for r in csv.DictReader(f)}


def read_sheet(wb, sheet: str, symbols: set[str]) -> list[tuple[str, str]]:
    """(label, gene_symbol) for every kept cluster x gene whose SCT value exceeds the cutoff."""
    ws = wb[sheet]
    it = ws.iter_rows(values_only=True)
    next(it)  # row 1: cluster numbers
    labels = [str(c).strip() if c is not None else "" for c in next(it)[1:]]  # row 2: annotations
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
            val = float(v) if isinstance(v, (int, float)) else 0.0
            if val > CUTOFF_T:
                out.add((ALIAS.get(raw, raw), sym))
    return sorted(out)


def main() -> None:
    path = mmc12_path()
    symbols = load_gene_symbols()
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    for sheet, (sex, fname) in SHEETS.items():
        rows = read_sheet(wb, sheet, symbols)
        with open(HERE / fname, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["cell", "gene_symbol"])
            w.writerows(rows)
        labels = sorted({c for c, _ in rows})
        genes = sorted({g for _, g in rows})
        print(f"wrote {fname}: {len(rows)} rows, {len(labels)} labels, {len(genes)} genes ({sex})")
    wb.close()


if __name__ == "__main__":
    main()
