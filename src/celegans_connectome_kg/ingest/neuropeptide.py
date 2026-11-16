"""Ingest the Ripoll-Sánchez et al. 2023 neuropeptidergic connectome (Neuron 111:3570).

A predicted extrasynaptic ("wireless") signaling network: a directed edge source->target where the
source expresses a neuropeptide precursor and the target a cognate GPCR, for a biochemically
validated NPP-GPCR pair (CeNGEN threshold-4 expression + EC50 <= 500 nM). Three range models
(short/mid/long) reflect assumed peptide diffusion distance.

Each vendored file is a labeled 302x302 directed weighted adjacency matrix (row = source, column =
target, weight = number of NPP-GPCR pathways). Zero-padded motor-neuron names (``DA01``) are
normalized to CIRCE names (``DA1``). This module only parses; the build mints Connection records
with connection_type ``neuropeptidergic``.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from celegans_connectome_kg.ingest.neuron_graph import ConnectionRecord

#: Zero-padded serial motor-neuron names in the matrix (DA01, VD09) -> CIRCE form (DA1, VD9).
_PAD = re.compile(r"^([A-Za-z]+)0(\d)$")


def _norm(name: str) -> str:
    n = name.strip()
    m = _PAD.match(n)
    return f"{m.group(1)}{m.group(2)}" if m else n


@dataclass(frozen=True)
class NeuropeptideNetwork:
    connections: list[ConnectionRecord]
    dataset_id: str
    dataset_name: str
    dataset_description: str
    sex: str


def read_neuropeptide_network(
    csv_path: Path, dataset_id: str, dataset_name: str, dataset_description: str
) -> NeuropeptideNetwork:
    """Read one range-model adjacency matrix into directed weighted ConnectionRecords."""
    rows = list(csv.reader(open(csv_path, newline="")))
    targets = [_norm(x) for x in rows[0][1:]]
    conns: list[ConnectionRecord] = []
    for r in rows[1:]:
        src = _norm(r[0])
        for j, val in enumerate(r[1:]):
            w = int(val) if val not in ("", None) else 0
            if w > 0:
                conns.append(
                    ConnectionRecord(
                        dataset_id=dataset_id,
                        pre=src,
                        post=targets[j],
                        connection_type="neuropeptidergic",
                        weight=float(w),
                        syn=(),
                        ids=None,
                        pre_tid=None,
                        post_tid=None,
                    )
                )
    return NeuropeptideNetwork(
        conns, dataset_id, dataset_name, dataset_description, "hermaphrodite"
    )


@dataclass(frozen=True)
class NeuropeptidePair:
    """One deorphanized NPP-ligand -> GPCR-receptor pair (mechanistic layer)."""

    index: int
    ligand: str
    gpcr: str
    ec50_nm: float | None
    npp_family: str
    gpcr_class: str


def read_neuropeptide_pairs(csv_path: Path) -> list[NeuropeptidePair]:
    """Read the 92 validated NPP-GPCR pairs (mechanistic/npp_gpcr_pairs.csv)."""
    out: list[NeuropeptidePair] = []
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            ec50 = r["ec50_nm"].strip()
            out.append(
                NeuropeptidePair(
                    index=int(r["pair_index"]),
                    ligand=r["ligand"].strip(),
                    gpcr=r["gpcr"].strip(),
                    ec50_nm=float(ec50) if ec50 else None,
                    npp_family=r["npp_family"].strip(),
                    gpcr_class=r["gpcr_class"].strip(),
                )
            )
    return out


def read_edge_pairs(csv_path: Path) -> dict[tuple[str, str], list[int]]:
    """Read edge_pairs.csv into {(source, target): [pair_index, ...]} (names already CIRCE-normed)."""
    out: dict[tuple[str, str], list[int]] = {}
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            out.setdefault((r["source"], r["target"]), []).append(int(r["pair_index"]))
    return out


@dataclass(frozen=True)
class NeuropeptideGene:
    """An NPP or GPCR gene of the network, resolved to WormBase."""

    symbol: str
    wbgene: str
    systematic_name: str
    category: str


def read_neuropeptide_genes(csv_path: Path) -> dict[str, NeuropeptideGene]:
    """Read np_genes.csv into {symbol: NeuropeptideGene} (WBGene as a WB:… curie)."""
    out: dict[str, NeuropeptideGene] = {}
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            out[r["symbol"].strip()] = NeuropeptideGene(
                symbol=r["symbol"].strip(),
                wbgene=r["wbgene"].strip(),
                systematic_name=r["systematic_name"].strip(),
                category=r["category"].strip(),
            )
    return out


def read_neuropeptide_expression(csv_path: Path) -> list[tuple[str, str]]:
    """Read np_gene_expression.csv into a list of (cell, gene_symbol) records."""
    with open(csv_path, newline="") as f:
        return [(r["cell"].strip(), r["gene_symbol"].strip()) for r in csv.DictReader(f)]


def derive_neuropeptide_edges(
    ligand_cells: dict[str, set[str]],
    receptor_cells: dict[str, set[str]],
    pairs: list[tuple[str, str]],
) -> dict[tuple[str, str], int]:
    """Derive a long-range neuropeptide connectome from expression + deorphanized pairs.

    The exact method underlying the Ripoll-Sánchez et al. 2023 *long-range* model: a directed edge
    source->target exists when, for some NPP-GPCR pair, the source expresses the ligand and the
    target the cognate GPCR; the edge weight is the number of such pairs. "Long-range" = no
    anatomical gating (every expressing source reaches every expressing target). Keys are whatever
    identifies a gene consistently across the three inputs (this codebase uses the bare WBGene id).

    ``ligand_cells``:  ligand-gene -> set of source cells expressing it.
    ``receptor_cells``: receptor-gene -> set of target cells expressing it.
    ``pairs``:          (ligand_gene, receptor_gene) for each deorphanized pair; a pair whose
                        ligand or receptor is absent from the expression maps contributes nothing.

    Summing CIRCE's own Ripoll expression over the 92 pairs with this function reproduces the
    published 53,558-edge hermaphrodite long-range network exactly (identical per-edge weights),
    so the same engine applied to another sex's expression is the validated method, not a new one.
    """
    edges: dict[tuple[str, str], int] = {}
    for lig, gpcr in pairs:
        srcs = ligand_cells.get(lig)
        tgts = receptor_cells.get(gpcr)
        if not srcs or not tgts:
            continue
        for s in srcs:
            for t in tgts:
                key = (s, t)
                edges[key] = edges.get(key, 0) + 1
    return edges
