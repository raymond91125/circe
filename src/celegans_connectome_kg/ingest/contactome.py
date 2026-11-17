"""Ingest a contactome (membrane-adjacency) edge list.

The contactome is the set of physical membrane contacts between neurons — observed appositions, not
synapses (a synapse can only form where membranes touch, so synapses are a subset of contacts). The
vendored edge list is undirected (``cell_a,cell_b,contact_size`` with ``cell_a < cell_b``); the build
mints ``adjacency`` Connection records in both orientations. See data/yim-2024-dauer.
"""

from __future__ import annotations

import csv
from pathlib import Path


def read_contactome(path: Path) -> list[tuple[str, str, float]]:
    """Read a contactome edge list into undirected (cell_a, cell_b, weight).

    Columns: ``cell_a,cell_b`` plus a third weight column named either ``contact_area`` (Yim 2024,
    EM-measured area) or ``contact_size`` (generic). Each row is one undirected contact.
    """
    out: list[tuple[str, str, float]] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        wcol = "contact_area" if "contact_area" in (reader.fieldnames or []) else "contact_size"
        for r in reader:
            out.append((r["cell_a"].strip(), r["cell_b"].strip(), float(r[wcol])))
    return out
