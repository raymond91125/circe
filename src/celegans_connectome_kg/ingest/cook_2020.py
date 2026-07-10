"""Ingest the Cook et al. 2020 C. elegans pharyngeal connectome.

Reads the vendored ``data/cook-2020-pharynx/edges.csv`` — the weighted edge list aggregated from
Supplemental Data 3 (the raw combined synapse list) by ``aggregate_edges.py``; see that directory's
MANIFEST for source, sha256, and why SI 3 is used over the SI 4 export / PDF. One sex-tagged
dataset, ``cook_2020_pharynx``.

The pharynx is shared between the sexes; the reconstruction is from hermaphrodite (N2) specimens,
so the dataset is tagged ``hermaphrodite`` (we make no sex-agnostic assumption).

Returns the same ``CookData`` shape as :mod:`cook_2019`, so ``assemble`` folds it in identically
(name reconciliation via ``cook_name_aliases``, sex-presence derivation, weight aggregation).

IMPORTANT — weight semantics: weights are **total EM serial sections** summed per pair across all
synapses/series — the SAME definition as Cook 2019 (synapse number x size), so the pharynx is
directly comparable to the 2019 pharyngeal subset. (Distinct from neuron-graph synapse counts; keep
tagged to its dataset.)
"""

from __future__ import annotations

import csv
from pathlib import Path

from celegans_connectome_kg.ingest.cook_2019 import CookData, CookDataset
from celegans_connectome_kg.ingest.neuron_graph import ConnectionRecord

#: Vendored aggregated pharyngeal edge list (relative to the repo root).
DEFAULT_COOK_2020_EDGES = (
    Path(__file__).resolve().parents[3] / "data" / "cook-2020-pharynx" / "edges.csv"
)

DATASET_ID = "cook_2020_pharynx"
_SEX = "hermaphrodite"
#: edge-list ``type`` -> our CellType-enum connection label.
_TYPE = {"chemical": "chemical", "electrical": "gap_junction"}


def read_cook_2020(edges_path: Path = DEFAULT_COOK_2020_EDGES) -> CookData:
    """Parse edges.csv into a one-dataset ``CookData`` (chemical directed, gap junctions undirected).

    Gap junctions are emitted once per unordered pair in the canonical (alphabetically-sorted)
    orientation, matching the KG's single-record-per-gap-junction convention; the source list
    already carries one row per pair (no reciprocals).
    """
    connections: list[ConnectionRecord] = []
    cells: set[str] = set()
    with Path(edges_path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            src, tgt = row["source"].strip(), row["target"].strip()
            ctype = _TYPE[row["type"].strip().lower()]
            weight = float(row["weight"])
            cells.add(src)
            cells.add(tgt)
            if ctype == "gap_junction":
                pre, post = sorted((src, tgt))  # canonical undirected orientation
            else:
                pre, post = src, tgt
            connections.append(
                ConnectionRecord(
                    dataset_id=DATASET_ID,
                    pre=pre,
                    post=post,
                    connection_type=ctype,
                    weight=weight,
                    syn=(),
                    ids=None,
                    pre_tid=None,
                    post_tid=None,
                )
            )
    dataset = CookDataset(
        id=DATASET_ID,
        sex=_SEX,
        name="Cook et al. 2020 (pharynx)",
        description=(
            "Pharyngeal connectome, Cook et al. 2020 (J Comp Neurol 528:2767-2784; PMC7601127); "
            "aggregated from Supplemental Data 3 (combined synapse list). Weights are total EM "
            "serial sections per pair across all synapses/series (same metric as Cook 2019). "
            "Reconstructed from hermaphrodite (N2) specimens."
        ),
    )
    return CookData(datasets=[dataset], connections=connections, cells_by_sex={_SEX: cells})
