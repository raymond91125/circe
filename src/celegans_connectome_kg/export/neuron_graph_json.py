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

from celegans_connectome_kg.build.assemble import CELL_PREFIX, DATASET_PREFIX, _wbbt_ancestry
from celegans_connectome_kg.match.wbbt import STRONG_KINDS, WBBTIndex

#: WBBT "pharyngeal cell" — all pharyngeal neurons/muscles/marginal/epithelial/gland cells is_a this.
_PHARYNGEAL_CELL = "WBbt:0005460"

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


def pharyngeal_cells(connectome: object, wbbt_path: object) -> list[str]:
    """Upper-cased cell names + classes whose WBbt term is_a "pharyngeal cell" (WBbt:0005460).

    The viz cell-info shows these as location "Pharynx": NemaNode's inhead/intail flags mark
    membership in the somatic head/tail *ganglia*, which excludes the pharyngeal nervous system,
    so pharyngeal cells otherwise fall through to the misleading "Body". Keys are upper-cased to
    match the viz's case-insensitive (uppercased) node lookup.
    """
    label, parents = _wbbt_ancestry(wbbt_path)

    def is_pharyngeal(curie: str) -> bool:
        seen: set[str] = set()
        stack = [curie]
        while stack:
            c = stack.pop()
            if c == _PHARYNGEAL_CELL:
                return True
            if c in seen:
                continue
            seen.add(c)
            stack.extend(parents.get(c, ()))
        return False

    members: dict[str, list[object]] = defaultdict(list)
    out: set[str] = set()
    for cell in connectome.cells:
        if cell.cell_class:
            members[cell.cell_class].append(cell)
        if cell.anatomy and is_pharyngeal(str(cell.anatomy)):
            out.add(cell.name.upper())
    for cls, mem in members.items():
        if any(m.anatomy and is_pharyngeal(str(m.anatomy)) for m in mem):
            out.add(cls.upper())
    return sorted(out)


def cell_sexes_map(connectome: object) -> dict[str, list[str]]:
    """Cell name → sexes present, for every cell the viz loads (hermaphrodite + male projections).

    Drives which viz database a cell is a valid node of (hermaphrodite views vs the additive male
    view) in cell-info.js. Restricted to *projected* cells — those with a NemaNode or synthesised
    type; class/endpoint stubs are not viz nodes. Keys and each cell's sexes are sorted for a
    stable output (``hermaphrodite`` sorts before ``male``).
    """
    viz_names = {c["name"] for c in cells_projection(connectome)} | {
        c["name"] for c in male_cells_projection(connectome)
    }
    out = {
        cell.name: sorted(str(s) for s in cell.sexes)
        for cell in connectome.cells
        if cell.name in viz_names and cell.sexes
    }
    return dict(sorted(out.items()))


def pharynx_database_cells(connectome: object, dataset_id: str = "cook_2020_pharynx") -> list[str]:
    """Upper-cased nodes of the Cook 2020 pharyngeal viz database (the additive "pharynx" database).

    Every cell appearing in the ``cook_2020_pharynx`` connections plus their classes, upper-cased to
    match the viz's case-insensitive node lookup. Populates ``validNodes['pharynx']`` so those cells
    are valid nodes there. Distinct from :func:`pharyngeal_cells` (WBbt "pharyngeal cell" ancestry,
    used for the *location* label): this is the concrete node set of one dataset's connectome.
    """
    class_of = {c.name: c.cell_class for c in connectome.cells}
    names: set[str] = set()
    for conn in connectome.connections:
        if _strip(str(conn.dataset), DATASET_PREFIX) != dataset_id:
            continue
        for end in (_strip(conn.pre, CELL_PREFIX), _strip(conn.post, CELL_PREFIX)):
            names.add(end.upper())
            cls = class_of.get(end)
            if cls:
                names.add(cls.upper())
    return sorted(names)


def anatomy_labels_map(terms: dict[str, str], index: WBBTIndex) -> dict[str, str]:
    """WBbt curie → human-readable term label, for every curie referenced in ``terms``.

    Lets the viz cell-info panel show the anatomy term name (e.g. "RIM" → "RIM neuron")
    beside the WBbt id, rather than only the opaque id. Curies absent from the ontology or
    without a label are omitted (the viz then shows just the id).
    """
    labels: dict[str, str] = {}
    for curie in sorted(set(terms.values())):
        term = index.terms.get(curie)
        if term and term.label:
            labels[curie] = term.label
    return labels


#: WormAtlas handbook pages are hand-curated (not ontology-derived). The neuron pattern is
#: keyed by neuron class; non-neuron categories need explicit handbook URLs.
_WORMATLAS_NEURON = "https://www.wormatlas.org/neurons/Individual%20Neurons/{}frameset.html"
_WORMATLAS_BODY_WALL_MUSCLE = (
    "https://wormatlas.org/hermaphrodite/musclesomatic/MusSomaticframeset.html"
)


def _is_body_wall_muscle(name: str) -> bool:
    u = name.upper()
    return u.startswith("BWM") or u in ("BODYWALLMUSCLES", "LEGACYBODYWALLMUSCLES")


def wormatlas_links_map(
    connectome: object, male_url_curation: dict[str, str] | None = None
) -> dict[str, str]:
    """Build a node-name → WormAtlas URL map for the viz cell-info link.

    WormAtlas links to handbook pages, which are not derivable from the cell name except for
    neurons. Neurons → the per-class neuron page; body wall muscles → the somatic-muscle
    handbook page; other non-neuron categories are omitted for now (the viz then renders no
    link rather than a broken neuron-pattern URL). Keys are upper-cased; lookup is
    case-insensitive in the viz.

    Non-neurons have no neuron-pattern URL; ``male_url_curation`` (cell name → exact WormAtlas URL)
    supplies it directly and takes precedence. For a curated cell that groups under a real
    ``cell_class`` (e.g. ``pm3D`` → class ``pm3``), the class key is emitted too, so the grouped
    node the viz actually clicks (the class) resolves as well as the individual cell. Cells absent
    from the curation fall back to the KG class / name neuron pattern as before.
    """
    male_url_curation = male_url_curation or {}
    members: dict[str, list[object]] = defaultdict(list)
    for cell in connectome.cells:
        if cell.cell_class:
            members[cell.cell_class].append(cell)

    def url_for(name: str, cell_type: str, cls: str | None) -> str | None:
        if _is_body_wall_muscle(name) or (cls and _is_body_wall_muscle(cls)):
            return _WORMATLAS_BODY_WALL_MUSCLE
        if cell_type == "neuron":
            return _WORMATLAS_NEURON.format(cls or name)
        return None

    out: dict[str, str] = {}
    for cell in connectome.cells:
        curated = male_url_curation.get(cell.name)
        url = curated or url_for(cell.name, str(cell.cell_type), cell.cell_class)
        if url:
            out[cell.name.upper()] = url
        # A curated grouped cell is clicked in the viz as its class; emit that key too.
        if curated and cell.cell_class:
            out.setdefault(cell.cell_class.upper(), curated)
    for cls, mem in members.items():
        url = url_for(cls, str(mem[0].cell_type), cls)
        if url:
            out.setdefault(cls.upper(), url)
    return dict(sorted(out.items()))


#: KG connection type → (out-relation, in-relation) codes for the class-level cell-info map.
#: Chemical & functional are directional; gap junctions are symmetric (single "e").
_KG_REL = {
    "chemical": ("o", "i"),
    "functional": ("fo", "fi"),
    "gap_junction": ("e", "e"),
}


def _kg_num(weight: float) -> float | int:
    """Integer EM section counts stay ints; functional (float) weights round to 3 places."""
    return int(weight) if float(weight).is_integer() else round(float(weight), 3)


def kg_connections_map(connectome: object) -> dict:
    """Class-level *complete* connectivity for the viz cell-info "Connections (knowledge graph)".

    The viz only draws connections for the currently-selected database at/above the per-type
    weight threshold, so a weak edge (below threshold) or one from a KG dataset not loaded into
    the viz DB looks like "no connections" — e.g. ``g2R``'s only input, ``M5→g2R`` (weight 1), is
    hidden by the default chemical threshold of 3. This map is the full KG connectivity so the
    panel can show what exists regardless of view/threshold, keyed by cell *class* (what the
    cell-info panel resolves the clicked node to).

    Shape (compact for bundling into the client):

        {"datasets": [dataset_id, ...],      # list index i ↔ code str(i)
         "conn": {class: {rel: {partner_class: {dataset_code: weight}}}}}

    ``rel`` is ``o``/``i`` (chemical out/in), ``e`` (gap junction, symmetric), or ``fo``/``fi``
    (functional out/in). Cells with no KG ``cell_class`` key by their own name (matching the male
    projection's individual display). Gap junctions are stored redundantly in both orientations in
    the KG (equal weights), so each unordered pair is counted once, then attributed to both
    endpoints' classes.
    """
    datasets = sorted({_strip(str(c.dataset), DATASET_PREFIX) for c in connectome.connections})
    code = {d: str(i) for i, d in enumerate(datasets)}
    cls = {c.name: (c.cell_class or c.name) for c in connectome.cells}

    def klass(curie: str) -> str:
        name = _strip(curie, CELL_PREFIX)
        return cls.get(name, name)

    # class -> rel -> partner_class -> dataset_code -> summed weight
    conn: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(float))))
    # Gap junctions first canonicalised per (dataset, unordered cell pair) to drop the redundant
    # reverse orientation, keyed at cell level so distinct member pairs still aggregate by class.
    gap_seen: set[tuple[str, str, str]] = set()

    for c in connectome.connections:
        ctype = str(c.connection_type)
        dset = code[_strip(str(c.dataset), DATASET_PREFIX)]
        pre_cls, post_cls = klass(c.pre), klass(c.post)
        if ctype == "gap_junction":
            a, b = _strip(c.pre, CELL_PREFIX), _strip(c.post, CELL_PREFIX)
            key = (dset, *sorted((a, b)))
            if key in gap_seen:
                continue
            gap_seen.add(key)
            conn[pre_cls]["e"][post_cls][dset] += c.weight
            if post_cls != pre_cls:
                conn[post_cls]["e"][pre_cls][dset] += c.weight
        else:
            rel_out, rel_in = _KG_REL[ctype]
            conn[pre_cls][rel_out][post_cls][dset] += c.weight
            conn[post_cls][rel_in][pre_cls][dset] += c.weight

    out = {
        cell: {
            rel: {
                p: {d: _kg_num(w) for d, w in sorted(ds.items())} for p, ds in sorted(parts.items())
            }
            for rel, parts in sorted(rels.items())
        }
        for cell, rels in sorted(conn.items())
    }
    return {"datasets": datasets, "conn": out}


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
        dataset = _strip(str(conn.dataset), DATASET_PREFIX)
        # The neuron-graph projection reproduces neuron-graph's own datasets only. Cook 2019
        # (the sex-aware KG datasets) are excluded here; a sex-aware / male viz projection is
        # a separate concern (M6), so the existing hermaphrodite viz feed is unchanged.
        if dataset.startswith("cook_"):
            continue
        pre = _strip(conn.pre, CELL_PREFIX)
        post = _strip(conn.post, CELL_PREFIX)
        api_type = _API_TYPE[str(conn.connection_type)]
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


# --- Male projection (M6): the sex-aware KG projected as a "male" viz database -----------------

#: KG cell_type → a neuron-graph `type` code for male-specific cells lacking a NemaNode type.
#: The subtype (sensory/inter/motor) isn't in the KG for Cook-only cells, so neurons get the
#: generic interneuron code; coloring by neurotransmitter is unaffected.
_MALE_TYPE_FALLBACK = {"neuron": "i", "muscle": "b", "other": ""}


def _viz_type(cell: object) -> str:
    if cell.nemanode_type:
        return str(cell.nemanode_type)
    return _MALE_TYPE_FALLBACK.get(str(cell.cell_type), "")


def male_cells_projection(connectome: object) -> list[dict]:
    """Project the cells present in the male (shared + male-specific) into the /api/cells shape.

    Unlike :func:`cells_projection` (neuron-graph only), this includes Cook male-specific cells,
    synthesizing a ``type`` from the KG ``cell_type`` and defaulting ``class`` to the cell name
    when the KG has no ``cell_class`` (so male-specific cells render individually).
    """
    out = []
    for c in connectome.cells:
        if "male" not in {str(s) for s in (c.sexes or [])}:
            continue
        out.append(
            {
                "name": c.name,
                "class": c.cell_class or c.name,
                "type": _viz_type(c),
                "neurotransmitter": c.neurotransmitter or "u",
                "embryonic": bool(c.embryonic),
                "inhead": bool(c.in_head),
                "intail": bool(c.in_tail),
            }
        )
    return out


def male_connections_projection(
    connectome: object, dataset_id: str = "cook_2019_male"
) -> list[dict]:
    """Project the male connectome (one dataset) into the /api/connections shape."""
    grouped: dict[tuple[str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for conn in connectome.connections:
        dataset = _strip(str(conn.dataset), DATASET_PREFIX)
        if dataset != dataset_id:
            continue
        pre = _strip(conn.pre, CELL_PREFIX)
        post = _strip(conn.post, CELL_PREFIX)
        api_type = _API_TYPE[str(conn.connection_type)]
        grouped[(pre, post, api_type)][dataset] += int(conn.weight)

    chemical, electrical = [], []
    for (pre, post, api_type), synapses in grouped.items():
        obj = {
            "pre": pre,
            "post": post,
            "type": api_type,
            "annotations": [],
            "synapses": dict(synapses),
        }
        (electrical if api_type == "electrical" else chemical).append(obj)
    return [*_merge_gap_junctions(electrical), *chemical]


def male_dataset(dataset_id: str = "cook_2019_male") -> dict:
    """The neuron-graph dataset entry for the male connectome (a 'male' viz database)."""
    return {
        "id": dataset_id,
        "name": "Cook et al. 2019 (male)",
        "type": "male",
        "time": None,
        "visualTime": None,
        "description": "Whole-animal male connectome, Cook et al. 2019 (Nature 571:63-71).",
        "datatypes": "cs,gj",
    }
