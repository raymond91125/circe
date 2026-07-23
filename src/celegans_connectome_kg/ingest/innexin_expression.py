"""Ingest the Bhattacharya et al. 2019 (Cell 176:1174-1189) neuronal innexin expression map.

The data is Figure 1B (a color-coded matrix, no machine-readable supplement), extracted from the
high-resolution figure and manually verified. Two vendored CSVs:

- ``innexin_expression.csv`` (``neuron, innexin, expression, developmental_plasticity``): per
  neuron-class × innexin, with expression scored ``both`` (non-dauer & dauer), ``non-dauer_only``,
  or ``dauer_only``.
- ``innexin_genes.csv`` (``fig_label, symbol, isoform, wbgene, category, systematic_name``): the 17
  neuronally-expressed innexins keyed to persistent WormBase gene ids (isoforms a/b as qualifiers).

This module only parses; the build maps neuron-class labels to member cells and mints reified
GeneExpression records across a non-dauer and a dauer dataset (so the plasticity is explicit).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class InnexinGene:
    wbgene: str
    symbol: str
    category: str
    systematic_name: str


@dataclass(frozen=True)
class InnexinExpr:
    neuron_label: str  # Bhattacharya class label, e.g. "ADA", "IL1-L/R", "DD", "VC"
    wbgene: str
    isoform: str  # "" when gene-level
    state: str  # "both" / "non-dauer_only" / "dauer_only"
    plasticity: bool  # Fig 1B "*": developmental plasticity within non-dauer


@dataclass(frozen=True)
class InnexinData:
    genes: list[InnexinGene]
    expressions: list[InnexinExpr]


def read_innexin_expression(expr_csv: Path, gene_map_csv: Path) -> InnexinData:
    """Parse the verified Fig 1B extraction + its resolved gene map."""
    gmap = {r["fig_label"]: r for r in csv.DictReader(open(gene_map_csv, newline=""))}
    genes = {
        r["wbgene"]: InnexinGene(r["wbgene"], r["symbol"], r["category"], r["systematic_name"])
        for r in gmap.values()
        if r["wbgene"]
    }
    exprs: list[InnexinExpr] = []
    for row in csv.DictReader(open(expr_csv, newline="")):
        g = gmap[row["innexin"]]
        exprs.append(
            InnexinExpr(
                neuron_label=row["neuron"].strip(),
                wbgene=g["wbgene"],
                isoform=g["isoform"],
                state=row["expression"].strip(),
                plasticity=row.get("developmental_plasticity", "").strip() == "yes",
            )
        )
    return InnexinData(genes=list(genes.values()), expressions=exprs)
