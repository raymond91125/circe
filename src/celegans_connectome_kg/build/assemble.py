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
from celegans_connectome_kg.match.curation import load_curation, load_endpoint_cells
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


def assemble(
    data_dir: Path,
    wbbt_path: Path,
    curation_path: Path | None = None,
    endpoint_cells_path: Path | None = None,
) -> tuple[object, BuildStats]:
    """Assemble a Connectome data-model object plus build stats.

    If ``curation_path`` is given, manual anatomy resolutions take precedence over lexical
    matches. If ``endpoint_cells_path`` is given, stub cells are minted for class-level /
    aggregate connection endpoints that are not in neuron-graph's cell list (see
    :mod:`celegans_connectome_kg.match.curation`).
    """
    dm: ModuleType = datamodel()
    data = load_neuron_graph(data_dir)
    index = WBBTIndex.from_obograph(wbbt_path)
    curation = load_curation(curation_path) if curation_path else None
    endpoint_cells = load_endpoint_cells(endpoint_cells_path) if endpoint_cells_path else []

    cell_records: dict[str, CellRecord] = {c.name: c for c in data.cells}
    anatomy_by_name = {
        m.cell_name: m.wbbt_id for m in match_cells(data.cells, index, curation) if m.wbbt_id
    }

    cells = [
        dm.Cell(
            id=_cell_id(c.name),
            name=c.name,
            cell_type=classify_cell_type(c.nemanode_type),
            anatomy=anatomy_by_name.get(c.name),
            cell_class=c.cell_class,
            neurotransmitter=c.neurotransmitter,
            nemanode_type=c.nemanode_type,
            embryonic=c.embryonic,
            in_head=c.in_head,
            in_tail=c.in_tail,
        )
        for c in sorted(data.cells, key=lambda c: c.name)
    ]

    # Stub cells for class-level / aggregate connection endpoints (KG-only; no neuron-graph
    # attributes, so the viz cells projection excludes them).
    stub_names = {ec.name for ec in endpoint_cells}
    cells += [
        dm.Cell(id=_cell_id(ec.name), name=ec.name, cell_type=ec.cell_type, anatomy=ec.wbbt_id)
        for ec in sorted(endpoint_cells, key=lambda ec: ec.name)
    ]

    datasets = [
        dm.Dataset(id=_dataset_id(d.id), name=d.name, description=d.description)
        for d in sorted(data.datasets, key=lambda d: d.id)
    ]

    # Aggregate by (dataset, pre, post, type), summing weight — mirrors neuron-graph's
    # populate-connections (`synapseCount += ...`), which collapses duplicate listings of the
    # same edge in one dataset into a single connection. One connection per key results.
    summed_weight: dict[ConnectionRecord, float] = {}
    by_key: dict[tuple[str, str, str, str], ConnectionRecord] = {}
    for conn in data.connections:
        key = (conn.dataset_id, conn.pre, conn.post, conn.connection_type)
        if key not in by_key:
            by_key[key] = conn
            summed_weight[conn] = 0.0
        summed_weight[by_key[key]] += conn.weight

    connections = []
    for conn in sorted(
        by_key.values(), key=lambda c: (c.dataset_id, c.pre, c.post, c.connection_type)
    ):
        weight = summed_weight[conn]
        connections.append(
            dm.Connection(
                id=_connection_id(conn),
                pre=_cell_id(conn.pre),
                post=_cell_id(conn.post),
                connection_type=conn.connection_type,
                weight=weight,
                dataset=_dataset_id(conn.dataset_id),
                evidence=[
                    dm.Evidence(
                        evidence_type="aggregated_weight",
                        dataset=_dataset_id(conn.dataset_id),
                        value=str(weight),
                    )
                ],
            )
        )

    referenced = {c.pre for c in by_key.values()} | {c.post for c in by_key.values()}
    known = set(cell_records) | stub_names
    unknown_cells = sorted(name for name in referenced if name not in known)

    connectome = dm.Connectome(cells=cells, datasets=datasets, connections=connections)
    stats = BuildStats(
        cells=len(cells),
        datasets=len(datasets),
        connections=len(connections),
        cells_with_anatomy=sum(1 for c in cells if c.anatomy),
        connections_by_type=dict(Counter(str(c.connection_type) for c in connections)),
        unknown_connection_cells=len(unknown_cells),
    )
    return connectome, stats
