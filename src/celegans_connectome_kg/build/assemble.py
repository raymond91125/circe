"""Phase 3 build: assemble the LinkML Connectome from ingest + match.

Orchestrates the earlier stages in-process — ingest the pinned neuron-graph snapshot, match
cell names to WBBT — and instantiates schema-conformant data model objects. Instance ids are
minted as ``cckg:`` CURIEs so the RDF export lands cleanly in the schema namespace.

Sex extension (M5): when Cook et al. 2019 inputs are supplied, the male and hermaphrodite
whole-animal connectomes are folded into the same graph as additional datasets tagged by
``Dataset.sex``. Cook cell names are reconciled to canonical KG names (``cook_name_aliases``);
Cook-only (male-specific / finer-grained) cells are minted and grounded via
``cook_anatomy_curation``. Every cell's ``sexes`` is derived from the sex(es) it is present in.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from types import ModuleType

from celegans_connectome_kg.build.datamodel import datamodel
from celegans_connectome_kg.ingest.cook_2019 import read_cook
from celegans_connectome_kg.ingest.cook_2020 import read_cook_2020
from celegans_connectome_kg.ingest.neuron_graph import (
    CellRecord,
    ConnectionRecord,
    load_neuron_graph,
)
from celegans_connectome_kg.match.curation import (
    load_cook_aliases,
    load_curation,
    load_endpoint_cells,
    load_nt_curation,
)
from celegans_connectome_kg.match.matcher import match_cells
from celegans_connectome_kg.match.wbbt import WBBTIndex

CELL_PREFIX = "cckg:cell/"
DATASET_PREFIX = "cckg:dataset/"
CONN_PREFIX = "cckg:conn/"

#: neuron-graph connection typ name → typ code, for stable connection ids.
_TYPE_CODE = {"chemical": "0", "gap_junction": "2", "functional": "4"}

#: Redundant neuron-graph datasets dropped at build. `randi_funconn_wildcp` is byte-identical to
#: `randi_funconn_wildty` — the same Randi wild-type functional recording that neuron-graph loads
#: twice (into its "complete" and "head" database collections). Keeping both would double-count
#: every wild-type functional weight in KG aggregates, so we keep one (`wildty`). The pinned
#: source snapshot still contains both; this is a build-time de-duplication.
_REDUNDANT_NG_DATASETS = frozenset({"randi_funconn_wildcp"})

_HERMAPHRODITE = "hermaphrodite"

#: neuron-graph endpoint name normalization. G1/G2 are uppercase misnomers of the pharyngeal
#: gland cells g1/g2 (gland cells are conventionally lowercase); rename so they unify with the
#: g1/g2 gland cells shared by the Cook datasets.
_NEURON_GRAPH_ALIAS = {"G1": "g1", "G2": "g2"}


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


def _wbbt_ancestry(wbbt_path: Path) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Return (curie → label, curie → is_a parents) from the OBO-graph JSON."""
    graph = json.loads(Path(wbbt_path).read_text())["graphs"][0]
    label: dict[str, str] = {}
    parents: dict[str, set[str]] = defaultdict(set)
    for node in graph["nodes"]:
        if "WBbt_" in node.get("id", ""):
            label["WBbt:" + node["id"].split("WBbt_")[1]] = node.get("lbl") or ""
    for edge in graph.get("edges", []):
        if (
            edge.get("pred") == "is_a"
            and "WBbt_" in edge.get("sub", "")
            and "WBbt_" in edge.get("obj", "")
        ):
            sub = "WBbt:" + edge["sub"].split("WBbt_")[1]
            obj = "WBbt:" + edge["obj"].split("WBbt_")[1]
            parents[sub].add(obj)
    return label, parents


def _cell_kind_from_wbbt(curie: str, label: dict[str, str], parents: dict[str, set[str]]) -> str:
    """Classify a cell as neuron / muscle / other from its WBBT term's ancestry labels."""
    seen: set[str] = set()
    stack = [curie]
    while stack:
        c = stack.pop()
        if c in seen:
            continue
        seen.add(c)
        stack.extend(parents.get(c, ()))
    text = " ".join(label.get(c, "") for c in seen).lower()
    if "neuron" in text:
        return "neuron"
    if "muscle" in text:
        return "muscle"
    return "other"


def _cell_id(name: str) -> str:
    return CELL_PREFIX + name


def _dataset_id(dataset_id: str) -> str:
    return DATASET_PREFIX + dataset_id


def _connection_id(dataset_id: str, pre: str, post: str, connection_type: str) -> str:
    return f"{CONN_PREFIX}{dataset_id}.{pre}.{post}.{_TYPE_CODE[connection_type]}"


@dataclass(frozen=True)
class BuildStats:
    """Summary of one build, for the CLI and tests."""

    cells: int
    datasets: int
    connections: int
    cells_with_anatomy: int
    connections_by_type: dict[str, int]
    unknown_connection_cells: int
    cells_by_sex: dict[str, int]
    datasets_by_sex: dict[str, int]


def _aggregate(records: list[ConnectionRecord], alias: dict[str, str]):
    """Sum weight by (dataset, pre, post, type); apply name aliases to endpoints.

    Mirrors neuron-graph's populate-connections (``synapseCount += ...``): duplicate listings
    of one edge in one dataset collapse into a single connection.
    """
    summed: dict[tuple[str, str, str, str], float] = {}
    for c in records:
        pre, post = alias.get(c.pre, c.pre), alias.get(c.post, c.post)
        key = (c.dataset_id, pre, post, c.connection_type)
        summed[key] = summed.get(key, 0.0) + c.weight
    return summed


def assemble(
    data_dir: Path,
    wbbt_path: Path,
    curation_path: Path | None = None,
    endpoint_cells_path: Path | None = None,
    nt_curation_path: Path | None = None,
    cook_xlsx_path: Path | None = None,
    cook_aliases_path: Path | None = None,
    cook_anatomy_path: Path | None = None,
    cook_2020_edges_path: Path | None = None,
) -> tuple[object, BuildStats]:
    """Assemble a Connectome data-model object plus build stats.

    Curation/endpoint/nt args behave as before. When ``cook_xlsx_path`` is given (with its
    ``cook_aliases_path`` and ``cook_anatomy_path``), the Cook et al. 2019 male + hermaphrodite
    connectomes are merged in as sex-tagged datasets, producing a unified sex-aware graph. When
    ``cook_2020_edges_path`` is given, the Cook et al. 2020 pharyngeal connectome is merged in as
    an additional (hermaphrodite) dataset over the same cells.
    """
    dm: ModuleType = datamodel()
    data = load_neuron_graph(data_dir)
    # Drop byte-identical duplicate datasets (see _REDUNDANT_NG_DATASETS) before assembly.
    data = replace(
        data,
        datasets=[d for d in data.datasets if d.id not in _REDUNDANT_NG_DATASETS],
        connections=[c for c in data.connections if c.dataset_id not in _REDUNDANT_NG_DATASETS],
    )
    index = WBBTIndex.from_obograph(wbbt_path)
    curation = load_curation(curation_path) if curation_path else None
    endpoint_cells = load_endpoint_cells(endpoint_cells_path) if endpoint_cells_path else []
    nt_curation = load_nt_curation(nt_curation_path) if nt_curation_path else {}

    cell_records: dict[str, CellRecord] = {c.name: c for c in data.cells}
    stub_names = {ec.name for ec in endpoint_cells}
    anatomy_by_name = {
        m.cell_name: m.wbbt_id for m in match_cells(data.cells, index, curation) if m.wbbt_id
    }

    # sex-presence per cell name (derived from which sexes each cell is present in)
    sexes: dict[str, set[str]] = defaultdict(set)
    for c in data.cells:
        sexes[c.name].add(_HERMAPHRODITE)  # neuron-graph is hermaphrodite
    for ec in endpoint_cells:
        sexes[ec.name].add(_HERMAPHRODITE)

    # --- Cook (optional): reconcile names, mint new cells, tag datasets by sex ---
    # Two bundles, folded identically: 2019 whole-animal (male + hermaphrodite) and 2020 pharynx.
    cook = read_cook(cook_xlsx_path) if cook_xlsx_path else None
    cook_2020 = read_cook_2020(cook_2020_edges_path) if cook_2020_edges_path else None
    cook_bundles = [b for b in (cook, cook_2020) if b]
    cook_alias = load_cook_aliases(cook_aliases_path) if cook_aliases_path else {}
    cook_anatomy = load_curation(cook_anatomy_path) if cook_anatomy_path else {}
    dataset_sex: dict[str, str] = {d.id: _HERMAPHRODITE for d in data.datasets}
    cook_new_cells: dict[str, str] = {}  # canonical name -> wbbt id (cells not in neuron-graph)
    for bundle in cook_bundles:
        for d in bundle.datasets:
            dataset_sex[d.id] = d.sex
        for sex, names in bundle.cells_by_sex.items():
            for n in names:
                canonical = cook_alias.get(n, n)
                sexes[canonical].add(sex)
                if canonical not in cell_records and canonical not in stub_names:
                    cook_new_cells[canonical] = cook_anatomy.get(canonical, "")

    kind_label, kind_parents = _wbbt_ancestry(wbbt_path)

    def sexes_of(name: str) -> list[str]:
        return sorted(sexes.get(name, set()))

    # --- Cells: neuron-graph cells, endpoint stubs, then Cook-only cells ---
    cells = [
        dm.Cell(
            id=_cell_id(c.name),
            name=c.name,
            cell_type=classify_cell_type(c.nemanode_type),
            anatomy=anatomy_by_name.get(c.name),
            cell_class=c.cell_class,
            neurotransmitter=nt_curation.get(c.name, c.neurotransmitter),
            nemanode_type=c.nemanode_type,
            embryonic=c.embryonic,
            in_head=c.in_head,
            in_tail=c.in_tail,
            sexes=sexes_of(c.name),
        )
        for c in sorted(data.cells, key=lambda c: c.name)
    ]
    cells += [
        dm.Cell(
            id=_cell_id(ec.name),
            name=ec.name,
            cell_type=ec.cell_type,
            anatomy=ec.wbbt_id,
            sexes=sexes_of(ec.name),
        )
        for ec in sorted(endpoint_cells, key=lambda ec: ec.name)
    ]
    cells += [
        dm.Cell(
            id=_cell_id(name),
            name=name,
            cell_type=_cell_kind_from_wbbt(wbbt, kind_label, kind_parents) if wbbt else "other",
            anatomy=wbbt or None,
            sexes=sexes_of(name),
        )
        for name, wbbt in sorted(cook_new_cells.items())
    ]

    # --- Datasets (with sex) ---
    dataset_defs = [(d.id, d.name, d.description) for d in data.datasets]
    for bundle in cook_bundles:
        dataset_defs += [(d.id, d.name, d.description) for d in bundle.datasets]
    datasets = [
        dm.Dataset(id=_dataset_id(did), name=name, description=desc, sex=dataset_sex[did])
        for did, name, desc in sorted(dataset_defs)
    ]

    # --- Connections: neuron-graph (G1/G2 gland-misnomer rename) + Cook bundles (aliased) ---
    summed = _aggregate(data.connections, _NEURON_GRAPH_ALIAS)
    for bundle in cook_bundles:
        for key, w in _aggregate(bundle.connections, cook_alias).items():
            summed[key] = summed.get(key, 0.0) + w

    connections = []
    for (did, pre, post, ctype), weight in sorted(summed.items()):
        connections.append(
            dm.Connection(
                id=_connection_id(did, pre, post, ctype),
                pre=_cell_id(pre),
                post=_cell_id(post),
                connection_type=ctype,
                weight=weight,
                dataset=_dataset_id(did),
                evidence=[
                    dm.Evidence(
                        evidence_type="aggregated_weight",
                        dataset=_dataset_id(did),
                        value=str(weight),
                    )
                ],
            )
        )

    referenced = {pre for _, pre, _, _ in summed} | {post for _, _, post, _ in summed}
    known = set(cell_records) | stub_names | set(cook_new_cells)
    unknown_cells = sorted(name for name in referenced if name not in known)

    connectome = dm.Connectome(cells=cells, datasets=datasets, connections=connections)
    stats = BuildStats(
        cells=len(cells),
        datasets=len(datasets),
        connections=len(connections),
        cells_with_anatomy=sum(1 for c in cells if c.anatomy),
        connections_by_type=dict(Counter(str(c.connection_type) for c in connections)),
        unknown_connection_cells=len(unknown_cells),
        cells_by_sex=dict(Counter("+".join(str(s) for s in c.sexes) or "none" for c in cells)),
        datasets_by_sex=dict(Counter(str(d.sex) for d in datasets)),
    )
    return connectome, stats
