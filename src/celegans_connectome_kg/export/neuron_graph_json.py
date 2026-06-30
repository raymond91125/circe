"""Phase 3 export: neuron-graph JSON projection (consumer 2 — the viz).

Projects the LinkML Connectome back into the exact shapes neuron-graph's HTTP API serves,
so the enriched, ontology-linked graph can feed the existing visualization:

  /api/cells       → bare array of
                     {name, class, type, neurotransmitter, embryonic, inhead, intail}
  /api/connections → bare array of {pre, post, type, annotations, synapses}
                     where synapses is {dataset_id: count} aggregated across datasets

The KG stores one Connection per (dataset, pre, post, type); the viz groups by
(pre, post, type) with a per-dataset synapse map. We reproduce neuron-graph's grouping,
its gap_junction→"electrical" label, and its additive gap-junction merge (flipped pre/post
pairs combined under the alphabetically-sorted orientation). This projection is lossless and
unfiltered — no per-type thresholds and no annotations (annotations are always []), matching
the API called over all cells/datasets.
"""

from __future__ import annotations

from collections import defaultdict

from celegans_connectome_kg.build.assemble import CELL_PREFIX, DATASET_PREFIX
from celegans_connectome_kg.match.wbbt import STRONG_KINDS, WBBTIndex

#: our CellType-enum connection label → neuron-graph API `type` string.
_API_TYPE = {"chemical": "chemical", "gap_junction": "electrical", "functional": "functional"}


def _strip(curie: str, prefix: str) -> str:
    return curie[len(prefix) :] if curie.startswith(prefix) else curie


def cells_projection(connectome: object) -> list[dict]:
    """Project Cells into the /api/cells shape (ordered, booleans as true/false).

    Stub cells minted for class-level connection endpoints are not neuron-graph cells (no
    NemaNode type) and are not in /api/cells, so they are excluded from this projection.
    """
    out = []
    for c in connectome.cells:
        if c.nemanode_type is None:
            continue
        out.append(
            {
                "name": c.name,
                "class": c.cell_class,
                "type": c.nemanode_type,
                "neurotransmitter": c.neurotransmitter,
                "embryonic": bool(c.embryonic),
                "inhead": bool(c.in_head),
                "intail": bool(c.in_tail),
            }
        )
    return out


def anatomy_terms_map(
    connectome: object,
    index: WBBTIndex,
    class_curation: dict[str, str] | None = None,
) -> dict[str, str]:
    """Build a node-name → WBbt map for the viz cell-info WormBase link.

    The cell-info panel links by whatever ``DataService.cellClass()`` returns for the clicked
    node — a cell *class* (e.g. ``AVA``) when grouped, or the cell/endpoint *name* (e.g.
    ``pm2D``) when shown individually or when the node isn't a neuron-graph cell. So the map
    must cover both: every cell name (incl. minted endpoint stubs, all grounded in the KG)
    and every class.

    Class resolution: manual class curation → unique strong (label/exact) WBBT class-name
    match → single-cell-class reuse of that cell's anatomy. Keys are upper-cased; the viz
    looks up case-insensitively. Entities with no WBbt term are omitted (the viz then renders
    no link rather than a broken one).
    """
    class_curation = class_curation or {}
    members: dict[str, list[str]] = defaultdict(list)
    cell_anatomy: dict[str, str] = {}
    for cell in connectome.cells:
        if cell.cell_class:
            members[cell.cell_class].append(cell.name)
        if cell.anatomy:
            cell_anatomy[cell.name] = str(cell.anatomy)

    out: dict[str, str] = {}
    # Individual cell / endpoint-stub names (handles individual display + non-neuron nodes).
    for name, wbbt in cell_anatomy.items():
        out[name.upper()] = wbbt
    # Cell classes (the default grouped view).
    for cls, mem in members.items():
        if cls in class_curation:
            wbbt = class_curation[cls]
        else:
            strong = {h.curie for h in index.lookup(cls) if h.kind in STRONG_KINDS}
            if len(strong) == 1:
                wbbt = next(iter(strong))
            elif len(mem) == 1 and mem[0] in cell_anatomy:
                wbbt = cell_anatomy[mem[0]]
            else:
                continue
        out[cls.upper()] = wbbt
    return dict(sorted(out.items()))


def _merge_gap_junctions(gap_junctions: list[dict]) -> list[dict]:
    """Merge flipped pre/post electrical pairs additively (neuron-graph mergeGapJunctions)."""
    by_key: dict[str, dict] = {}
    for gj in gap_junctions:
        key = "$".join(sorted([gj["pre"], gj["post"]]))
        if key not in by_key:
            # Copy so accumulation does not mutate the source dict's synapses.
            by_key[key] = {**gj, "synapses": dict(gj["synapses"])}
        else:
            for dataset, count in gj["synapses"].items():
                by_key[key]["synapses"][dataset] = by_key[key]["synapses"].get(dataset, 0) + count

    merged = []
    for key, gj in by_key.items():
        pre, post = key.split("$")
        merged.append(
            {
                "pre": pre,
                "post": post,
                "type": gj["type"],
                "annotations": [],
                "synapses": gj["synapses"],
            }
        )
    return merged


def connections_projection(connectome: object) -> list[dict]:
    """Project Connections into the /api/connections shape.

    Group by (pre, post, api_type) into a {dataset: count} synapse map, then merge gap
    junctions. Final order: merged gap junctions, then chemical, then functional.
    """
    grouped: dict[tuple[str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for conn in connectome.connections:
        pre = _strip(conn.pre, CELL_PREFIX)
        post = _strip(conn.post, CELL_PREFIX)
        api_type = _API_TYPE[str(conn.connection_type)]
        dataset = _strip(str(conn.dataset), DATASET_PREFIX)
        grouped[(pre, post, api_type)][dataset] += int(conn.weight)

    chemical, electrical, functional = [], [], []
    for (pre, post, api_type), synapses in grouped.items():
        obj = {
            "pre": pre,
            "post": post,
            "type": api_type,
            "annotations": [],
            "synapses": dict(synapses),
        }
        if api_type == "electrical":
            electrical.append(obj)
        elif api_type == "chemical":
            chemical.append(obj)
        else:
            functional.append(obj)

    return [*_merge_gap_junctions(electrical), *chemical, *functional]
