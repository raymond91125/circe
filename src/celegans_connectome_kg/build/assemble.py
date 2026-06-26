"""Phase 3 build: assemble the LinkML Connectome from ingest + match.

Orchestrates the earlier stages in-process — ingest the pinned neuron-graph snapshot, match
cell names to WBBT — and instantiates schema-conformant data model objects. Instance ids are
minted as ``cckg:`` CURIEs so the RDF export lands cleanly in the schema namespace.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType

from celegans_connectome_kg.build.datamodel import datamodel
from celegans_connectome_kg.ingest.neuron_graph import (
    CellRecord,
    ConnectionRecord,
    load_neuron_graph,
)
from celegans_connectome_kg.match.matcher import match_cells
from celegans_connectome_kg.match.wbbt import WBBTIndex

CELL_PREFIX = "cckg:cell/"
DATASET_PREFIX = "cckg:dataset/"
CONN_PREFIX = "cckg:conn/"

#: neuron-graph connection typ name → typ code, for stable connection ids.
_TYPE_CODE = {"chemical": "0", "gap_junction": "2", "functional": "4"}


def classify_cell_type(nemanode_type: str | None) -> str:
    """Map the NemaNode ``typ`` code to our CellType enum.

    Mirrors neuron-graph's ``getTypeDisplayNames``: ``b`` is muscle; any code carrying a
    neuron letter (s/i/m/n) is a neuron; everything else (e.g. ``u`` marginal cells, empty)
    is non-neuronal / non-muscle → ``other``.
    """
    code = nemanode_type or ""
    if code == "b":
        return "muscle"
    if any(letter in code for letter in "simn"):
        return "neuron"
    return "other"


def _cell_id(name: str) -> str:
    return CELL_PREFIX + name


def _dataset_id(dataset_id: str) -> str:
    return DATASET_PREFIX + dataset_id


def _connection_id(conn: ConnectionRecord) -> str:
    return (
        f"{CONN_PREFIX}{conn.dataset_id}.{conn.pre}.{conn.post}.{_TYPE_CODE[conn.connection_type]}"
    )


@dataclass(frozen=True)
class BuildStats:
    """Summary of one build, for the CLI and tests."""

    cells: int
    datasets: int
    connections: int
    cells_with_anatomy: int
    connections_by_type: dict[str, int]
    unknown_connection_cells: int


def assemble(data_dir: Path, wbbt_path: Path) -> tuple[object, BuildStats]:
    """Assemble a Connectome data-model object plus build stats."""
    dm: ModuleType = datamodel()
    data = load_neuron_graph(data_dir)
    index = WBBTIndex.from_obograph(wbbt_path)

    cell_records: dict[str, CellRecord] = {c.name: c for c in data.cells}
    anatomy_by_name = {m.cell_name: m.wbbt_id for m in match_cells(data.cells, index) if m.wbbt_id}

    cells = [
        dm.Cell(
            id=_cell_id(c.name),
            name=c.name,
            cell_type=classify_cell_type(c.nemanode_type),
            anatomy=anatomy_by_name.get(c.name),
        )
        for c in sorted(data.cells, key=lambda c: c.name)
    ]

    datasets = [
        dm.Dataset(id=_dataset_id(d.id), name=d.name, description=d.description)
        for d in sorted(data.datasets, key=lambda d: d.id)
    ]

    connections = []
    seen_ids: set[str] = set()
    unknown_cells = 0
    for conn in sorted(
        data.connections, key=lambda c: (c.dataset_id, c.pre, c.post, c.connection_type)
    ):
        conn_id = _connection_id(conn)
        if conn_id in seen_ids:
            # Same (dataset, pre, post, type) twice — should not happen; disambiguate.
            suffix = 1
            while f"{conn_id}#{suffix}" in seen_ids:
                suffix += 1
            conn_id = f"{conn_id}#{suffix}"
        seen_ids.add(conn_id)
        if conn.pre not in cell_records:
            unknown_cells += 1
        if conn.post not in cell_records:
            unknown_cells += 1
        connections.append(
            dm.Connection(
                id=conn_id,
                pre=_cell_id(conn.pre),
                post=_cell_id(conn.post),
                connection_type=conn.connection_type,
                weight=conn.weight,
                dataset=_dataset_id(conn.dataset_id),
                evidence=[
                    dm.Evidence(
                        evidence_type="aggregated_weight",
                        dataset=_dataset_id(conn.dataset_id),
                        value=str(conn.weight),
                    )
                ],
            )
        )

    connectome = dm.Connectome(cells=cells, datasets=datasets, connections=connections)
    stats = BuildStats(
        cells=len(cells),
        datasets=len(datasets),
        connections=len(connections),
        cells_with_anatomy=sum(1 for c in cells if c.anatomy),
        connections_by_type=dict(Counter(str(c.connection_type) for c in connections)),
        unknown_connection_cells=unknown_cells,
    )
    return connectome, stats
