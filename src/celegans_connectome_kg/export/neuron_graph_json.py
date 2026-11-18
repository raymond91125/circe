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
_API_TYPE = {
    "chemical": "chemical",
    "gap_junction": "electrical",
    "functional": "functional",
    "neuropeptidergic": "neuropeptidergic",
    "monoaminergic": "monoaminergic",
    "adjacency": "adjacency",
}


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
#: Chemical, functional & neuropeptidergic are directional; gap junctions are symmetric ("e").
_KG_REL = {
    "chemical": ("o", "i"),
    "functional": ("fo", "fi"),
    "gap_junction": ("e", "e"),
    "neuropeptidergic": ("npo", "npi"),
    "monoaminergic": ("mao", "mai"),
}


def _kg_num(weight: float) -> float | int:
    """Integer EM section counts stay ints; functional (float) weights round to 3 places."""
    return int(weight) if float(weight).is_integer() else round(float(weight), 3)


def _signed_weight(conn: object) -> float:
    """Viz-facing signed weight. The KG stores a non-negative magnitude + an explicit ``sign``; the
    neuron-graph projections mirror the original convention of a *signed* functional weight
    (- for inhibitory), so returning markers/edges downstream are unchanged."""
    w = conn.weight
    if str(conn.connection_type) == "functional" and str(getattr(conn, "sign", None) or "") == (
        "inhibitory"
    ):
        return -w
    return w


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

    ``rel`` is ``o``/``i`` (chemical out/in), ``e`` (gap junction, symmetric), ``fo``/``fi``
    (functional out/in), ``npo``/``npi`` (predicted neuropeptide out/in), ``mao``/``mai`` (predicted
    monoamine out/in), or ``adj`` (membrane contact / contactome, symmetric; weight = contact area).
    Cells with no KG ``cell_class`` key by their own name (matching the male projection's individual
    display). The symmetric relations (gap junction ``e``, contact ``adj``) are stored redundantly in
    both orientations in the KG, so each unordered pair is counted once, then attributed to both
    endpoints' classes.

    The predicted neuropeptidergic network is included so the cell-info panel is complete ("all
    partners across every dataset"); it aggregates to class level here, so the bundled map stays
    far smaller than the ~125k raw edges.
    """
    conns = connectome.connections
    # Datasets that contribute a surfaced relation get a code. Membrane adjacency (the contactome) is
    # now surfaced as the symmetric ``adj`` relation, so the contactome datasets are included (and are
    # referenced by ``conn`` -- not orphan codes).
    kg_types = set(_KG_REL) | {"adjacency"}
    datasets = sorted(
        {
            _strip(str(c.dataset), DATASET_PREFIX)
            for c in conns
            if str(c.connection_type) in kg_types
        }
    )
    code = {d: str(i) for i, d in enumerate(datasets)}
    cls = {c.name: (c.cell_class or c.name) for c in connectome.cells}

    def klass(curie: str) -> str:
        name = _strip(curie, CELL_PREFIX)
        return cls.get(name, name)

    # class -> rel -> partner_class -> dataset_code -> summed weight
    conn: dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(float))))
    # Symmetric relations -- gap junction (``e``) and membrane contact (``adj``) -- are stored in both
    # orientations in the KG, so canonicalise per (rel, dataset, unordered cell pair) to drop the
    # redundant reverse, keyed at cell level so distinct member pairs still aggregate by class.
    symmetric = {"gap_junction": "e", "adjacency": "adj"}
    sym_seen: set[tuple[str, str, str, str]] = set()

    for c in conns:
        ctype = str(c.connection_type)
        if ctype not in kg_types:
            continue
        dset = code[_strip(str(c.dataset), DATASET_PREFIX)]
        pre_cls, post_cls = klass(c.pre), klass(c.post)
        if ctype in symmetric:
            rel = symmetric[ctype]
            a, b = _strip(c.pre, CELL_PREFIX), _strip(c.post, CELL_PREFIX)
            key = (rel, dset, *sorted((a, b)))
            if key in sym_seen:
                continue
            sym_seen.add(key)
            conn[pre_cls][rel][post_cls][dset] += c.weight
            if post_cls != pre_cls:
                conn[post_cls][rel][pre_cls][dset] += c.weight
        else:
            rel_out, rel_in = _KG_REL[ctype]
            w = _signed_weight(
                c
            )  # functional inhibitory edges stay negative in the cell-info panel
            conn[pre_cls][rel_out][post_cls][dset] += w
            conn[post_cls][rel_in][pre_cls][dset] += w

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
        # The hermaphrodite neuron-graph projection reproduces neuron-graph's own datasets only.
        # Allowlist the neuron-graph-native ids so KG-added datasets (Cook 2019/2020, Bhatla 2015,
        # and any future addition) never leak into the viz's complete/head/tail databases — where
        # their differing weight scales would be mixed with the White/Witvliet synapse counts.
        if not dataset.startswith(("white_1986_", "witvliet_2020_", "randi_funconn_")):
            continue
        pre = _strip(conn.pre, CELL_PREFIX)
        post = _strip(conn.post, CELL_PREFIX)
        api_type = _API_TYPE[str(conn.connection_type)]
        grouped[(pre, post, api_type)][dataset] += int(_signed_weight(conn))

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
    # Prefer the male-specific neurotransmitter call (Wang 2024) where the KG has one — this is
    # the only NT for male-specific neurons, and the male call for sexually-dimorphic neurons.
    male_nt = {
        _strip(str(a.cell), CELL_PREFIX): str(a.neurotransmitter)
        for a in (getattr(connectome, "neurotransmitter_assignments", None) or [])
        if str(a.sex) == "male" and a.neurotransmitter
    }
    out = []
    for c in connectome.cells:
        if "male" not in {str(s) for s in (c.sexes or [])}:
            continue
        out.append(
            {
                "name": c.name,
                "class": c.cell_class or c.name,
                "type": _viz_type(c),
                "neurotransmitter": male_nt.get(c.name) or c.neurotransmitter or "u",
                "embryonic": bool(c.embryonic),
                "inhead": bool(c.in_head),
                "intail": bool(c.in_tail),
            }
        )
    return out


def _single_dataset_connections(connectome: object, dataset_id: str) -> list[dict]:
    """Project one dataset's connections into the /api/connections shape (chemical + gap only)."""
    grouped: dict[tuple[str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for conn in connectome.connections:
        dataset = _strip(str(conn.dataset), DATASET_PREFIX)
        if dataset != dataset_id:
            continue
        pre = _strip(conn.pre, CELL_PREFIX)
        post = _strip(conn.post, CELL_PREFIX)
        api_type = _API_TYPE[str(conn.connection_type)]
        grouped[(pre, post, api_type)][dataset] += int(_signed_weight(conn))

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


def male_connections_projection(
    connectome: object, dataset_id: str = "cook_2019_male"
) -> list[dict]:
    """Project the male connectome (one dataset) into the /api/connections shape."""
    return _single_dataset_connections(connectome, dataset_id)


def male_dataset(dataset_id: str = "cook_2019_male") -> dict:
    """The neuron-graph dataset entry for the male connectome (a 'male' viz database)."""
    return {
        "id": dataset_id,
        "name": "Cook et al. 2019 (male)",
        "type": "male",
        "time": 50,
        "visualTime": 50,
        "description": "Whole-animal male connectome, Cook et al. 2019 (Nature 571:63-71).",
        "datatypes": "cs,gj",
    }


# --- Pharynx projection: the Cook 2020 pharyngeal connectome as a "pharynx" viz database --------


def pharynx_cells_projection(
    connectome: object, dataset_id: str = "cook_2020_pharynx"
) -> list[dict]:
    """Project the cells of the Cook 2020 pharyngeal connectome into the /api/cells shape.

    The cell set is every endpoint of a ``cook_2020_pharynx`` connection (the full-record analogue
    of :func:`pharynx_database_cells`), sharing :func:`_viz_type` and class defaulting with the
    male projection.
    """
    endpoints: set[str] = set()
    for conn in connectome.connections:
        if _strip(str(conn.dataset), DATASET_PREFIX) != dataset_id:
            continue
        endpoints.add(_strip(conn.pre, CELL_PREFIX))
        endpoints.add(_strip(conn.post, CELL_PREFIX))
    return [
        {
            "name": c.name,
            "class": c.cell_class or c.name,
            "type": _viz_type(c),
            "neurotransmitter": c.neurotransmitter or "u",
            "embryonic": bool(c.embryonic),
            "inhead": bool(c.in_head),
            "intail": bool(c.in_tail),
        }
        for c in connectome.cells
        if c.name in endpoints
    ]


def pharynx_connections_projection(
    connectome: object, dataset_id: str = "cook_2020_pharynx"
) -> list[dict]:
    """Project the pharyngeal connectome (one dataset) into the /api/connections shape."""
    return _single_dataset_connections(connectome, dataset_id)


def pharynx_dataset(dataset_id: str = "cook_2020_pharynx") -> dict:
    """The neuron-graph dataset entry for the pharyngeal connectome (a 'pharynx' viz database)."""
    return {
        "id": dataset_id,
        "name": "Cook et al. 2020 (pharynx)",
        "type": "pharynx",
        "time": 50,
        "visualTime": 50,
        "description": "Pharyngeal connectome, Cook et al. 2020 (J Comp Neurol 528:2767-2784).",
        "datatypes": "cs,gj",
    }


# --- Dauer projection: the Yim 2024 dauer nerve-ring connectome as a head life-stage dataset -----

#: Dauer is an alternative L3-stage larva, off the continuous L1->adult clock. It is folded into
#: the "head" developmental series (a nerve-ring reconstruction with synapse-count weights, like
#: the Witvliet head datasets) and positioned in the L3 region of the shared life-stage timeline.
_DAUER_VISUAL_TIME = 30


def dauer_cells_projection(connectome: object, dataset_id: str = "yim_2024_dauer") -> list[dict]:
    """Project the cells of the Yim 2024 dauer connectome into the /api/cells shape."""
    endpoints: set[str] = set()
    for conn in connectome.connections:
        if _strip(str(conn.dataset), DATASET_PREFIX) != dataset_id:
            continue
        endpoints.add(_strip(conn.pre, CELL_PREFIX))
        endpoints.add(_strip(conn.post, CELL_PREFIX))
    return [
        {
            "name": c.name,
            "class": c.cell_class or c.name,
            "type": _viz_type(c),
            "neurotransmitter": c.neurotransmitter or "u",
            "embryonic": bool(c.embryonic),
            "inhead": bool(c.in_head),
            "intail": bool(c.in_tail),
        }
        for c in connectome.cells
        if c.name in endpoints
    ]


def dauer_connections_projection(
    connectome: object, dataset_id: str = "yim_2024_dauer"
) -> list[dict]:
    """Project the dauer connectome (one dataset) into the /api/connections shape."""
    return _single_dataset_connections(connectome, dataset_id)


def dauer_dataset(dataset_id: str = "yim_2024_dauer") -> dict:
    """The neuron-graph dataset entry for the dauer connectome.

    Folded into the "head" collection so it is browsable as a life stage alongside the Witvliet
    L1->adult series; positioned in the L3 region of the shared timeline (dauer branches near L3).
    """
    return {
        "id": dataset_id,
        "name": "Yim et al. 2024 (dauer)",
        "type": "head",
        "time": _DAUER_VISUAL_TIME,
        "visualTime": _DAUER_VISUAL_TIME,
        "description": "Dauer nerve-ring connectome, Yim et al. 2024 (Nat Commun 15:1546).",
        "datatypes": "cs",
    }


# --- Neuropeptide projection: the Ripoll-Sánchez 2023 predicted neuropeptidergic network ---------
#
# The extrasynaptic neuropeptide connectome is a *predicted* signalling network (CeNGEN
# single-cell expression of NPP ligands + GPCR receptors, filtered by EC50 ≤ 500 nM), not observed
# synapses. It is surfaced in the viz as a distinct "np" edge datatype layered on the whole-animal
# hermaphrodite ("complete") database, deliberately kept off by default (see default-datasets.js)
# so its density never overwhelms the observed wiring. Three diffusion range models are projected
# as separate selectable datasets; short-range is the primary one.
#
# The viz `datasets.id` column is varchar(20); the KG dataset ids (``ripoll_2023_neuropeptide_sr``,
# 27 chars) are too long, so each projects under a short viz id (``ripoll_2023_np_sr``).
#
# Each range model is projected as its OWN viz database (collection np_sr/np_mr/np_lr), rather
# than three datasets crammed onto the "complete" timeline where their bookmarks and tooltips
# overlapped and were hard to read/select. Each database then holds a single dataset (one clean
# tick at the adult end), selected via the database dropdown — consistent with the male/pharynx/
# dauer databases. The circe-vis side renders the single tick as a distinct amber circle marker.
# Entries: (KG dataset id, viz dataset id, viz collection, viz name, one-line description).
_NEUROPEPTIDE_DATASETS = [
    (
        "ripoll_2023_neuropeptide_sr",
        "ripoll_2023_np_sr",
        "np_sr",
        "Ripoll-Sánchez 2023 (NP short-range, predicted)",
        "Predicted short-range neuropeptidergic connectome, Ripoll-Sánchez et al. 2023 "
        "(Neuron 111:3570-3589). Extrasynaptic — not observed synapses.",
    ),
    (
        "ripoll_2023_neuropeptide_mr",
        "ripoll_2023_np_mr",
        "np_mr",
        "Ripoll-Sánchez 2023 (NP mid-range, predicted)",
        "Predicted mid-range neuropeptidergic connectome, Ripoll-Sánchez et al. 2023 "
        "(Neuron 111:3570-3589). Extrasynaptic — not observed synapses.",
    ),
    (
        "ripoll_2023_neuropeptide_lr",
        "ripoll_2023_np_lr",
        "np_lr",
        "Ripoll-Sánchez 2023 (NP long-range, predicted)",
        "Predicted long-range neuropeptidergic connectome, Ripoll-Sánchez et al. 2023 "
        "(Neuron 111:3570-3589). Extrasynaptic — not observed synapses.",
    ),
]

#: All three predicted networks are modelled on the adult hermaphrodite → one adult tick each.
_NEUROPEPTIDE_VISUAL_TIME = 50

#: KG dataset id → short viz dataset id (the projected synapses map + dataset entries use the viz id).
_NEUROPEPTIDE_VIZ_ID = {kg_id: viz_id for kg_id, viz_id, *_rest in _NEUROPEPTIDE_DATASETS}


def _neuropeptide_dataset_ids() -> list[str]:
    return [d[0] for d in _NEUROPEPTIDE_DATASETS]


def neuropeptide_cells_projection(connectome: object) -> list[dict]:
    """Project the neurons of the neuropeptide network(s) into the /api/cells shape.

    The cell set is every endpoint of a ``ripoll_2023_neuropeptide_*`` connection (all standard
    hermaphrodite neurons); shares :func:`_viz_type` and class defaulting with the male projection.
    """
    ids = set(_neuropeptide_dataset_ids())
    endpoints: set[str] = set()
    for conn in connectome.connections:
        if _strip(str(conn.dataset), DATASET_PREFIX) not in ids:
            continue
        endpoints.add(_strip(conn.pre, CELL_PREFIX))
        endpoints.add(_strip(conn.post, CELL_PREFIX))
    return [
        {
            "name": c.name,
            "class": c.cell_class or c.name,
            "type": _viz_type(c),
            "neurotransmitter": c.neurotransmitter or "u",
            "embryonic": bool(c.embryonic),
            "inhead": bool(c.in_head),
            "intail": bool(c.in_tail),
        }
        for c in connectome.cells
        if c.name in endpoints
    ]


def neuropeptide_connections_projection(connectome: object) -> list[dict]:
    """Project the neuropeptide network(s) into the /api/connections shape.

    Neuropeptidergic edges are directed and predicted (no gap-junction merge, no per-synapse ids).
    Grouped by (pre, post) with a per-dataset weight map, keeping each range model's weight under
    its own dataset id so the viz can show/select them independently.
    """
    ids = set(_neuropeptide_dataset_ids())
    grouped: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for conn in connectome.connections:
        if str(conn.connection_type) != "neuropeptidergic":
            continue
        dataset = _strip(str(conn.dataset), DATASET_PREFIX)
        if dataset not in ids:
            continue
        pre = _strip(conn.pre, CELL_PREFIX)
        post = _strip(conn.post, CELL_PREFIX)
        grouped[(pre, post)][_NEUROPEPTIDE_VIZ_ID[dataset]] += int(conn.weight)
    return [
        {
            "pre": pre,
            "post": post,
            "type": "neuropeptidergic",
            "annotations": [],
            "synapses": dict(synapses),
        }
        for (pre, post), synapses in grouped.items()
    ]


def neuropeptide_database_cells(connectome: object) -> list[str]:
    """Upper-cased node set (cell names + classes) of the predicted neuropeptide network.

    Populates ``validNodes`` for each neuropeptide viz database (np_sr/np_mr/np_lr — all share the
    same 302 neurons), the analogue of :func:`pharynx_database_cells` for those databases. Upper-
    cased to match the viz's case-insensitive node lookup.
    """
    class_of = {c.name: c.cell_class for c in connectome.cells}
    ids = set(_neuropeptide_dataset_ids())
    names: set[str] = set()
    for conn in connectome.connections:
        if _strip(str(conn.dataset), DATASET_PREFIX) not in ids:
            continue
        for end in (_strip(conn.pre, CELL_PREFIX), _strip(conn.post, CELL_PREFIX)):
            names.add(end.upper())
            cls = class_of.get(end)
            if cls:
                names.add(cls.upper())
    return sorted(names)


def neuropeptide_pairs_map(connectome: object, edge_pairs_path: object) -> dict:
    """Class-level NPP-GPCR pair attribution for the viz cell-info neuropeptide sections.

    Answers "which peptide->receptor pairs mediate this predicted edge". The 92 pairs come from the
    KG (``connectome.neuropeptide_receptor_pairs``); the per-edge attribution is read from the
    vendored ``mechanistic/edge_pairs.csv`` (kept out of RDF — 145k rows). Cell-level edges are
    aggregated to the class level the cell-info panel resolves nodes to, taking the union of pair
    indices across member-cell edges.

    Shape (compact, lazy-loaded by the client):

        {"pairs": [[ligand, gpcr, ec50_nm, gpcr_class], ...],   # list index i <-> pair
         "conn":  {source_class: {target_class: [pair_i, ...]}}}

    Class keys are upper-cased (case-insensitive viz lookup); a cell with no KG ``cell_class`` keys
    by its own upper-cased name. ``pair_i`` is the 0-based index into ``pairs``.
    """
    import csv

    pairs = sorted(
        getattr(connectome, "neuropeptide_receptor_pairs", None) or [],
        key=lambda p: int(str(p.id).rsplit("/", 1)[-1]),
    )
    # 1-based pair_index (from edge_pairs.csv) -> 0-based position in the emitted list.
    index_of = {int(str(p.id).rsplit("/", 1)[-1]): i for i, p in enumerate(pairs)}
    pairs_out = [
        [
            p.ligand,
            p.gpcr,
            (float(p.ec50_nm) if p.ec50_nm is not None else None),
            str(p.gpcr_class or ""),
        ]
        for p in pairs
    ]

    cls = {c.name: (c.cell_class or c.name) for c in connectome.cells}

    def klass(name: str) -> str:
        return cls.get(name, name).upper()

    conn: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    with open(edge_pairs_path, newline="") as f:
        for r in csv.DictReader(f):
            i = index_of.get(int(r["pair_index"]))
            if i is None:
                continue
            conn[klass(r["source"])][klass(r["target"])].add(i)

    out = {
        s: {t: sorted(idxs) for t, idxs in sorted(tgts.items())} for s, tgts in sorted(conn.items())
    }
    return {"pairs": pairs_out, "conn": out}


def monoamine_pairs_map(connectome: object, edge_pairs_path: object) -> dict:
    """Class-level monoamine->receptor pair attribution for the viz cell-info monoamine sections.

    The monoamine analogue of ``neuropeptide_pairs_map``: answers "which monoamine->receptor pairs
    mediate this predicted edge". The 14 pairs come from the KG
    (``connectome.monoamine_receptor_pairs``); the per-edge attribution is read from the vendored
    ``mechanistic/edge_pairs.csv``, aggregated to the class level with a union of pair indices.

    Shape (compact, lazy-loaded by the client):

        {"pairs": [[monoamine, receptor], ...],   # list index i <-> pair
         "conn":  {source_class: {target_class: [pair_i, ...]}}}
    """
    import csv

    pairs = sorted(
        getattr(connectome, "monoamine_receptor_pairs", None) or [],
        key=lambda p: int(str(p.id).rsplit("/", 1)[-1]),
    )
    index_of = {int(str(p.id).rsplit("/", 1)[-1]): i for i, p in enumerate(pairs)}
    pairs_out = [[p.monoamine, p.receptor] for p in pairs]

    cls = {c.name: (c.cell_class or c.name) for c in connectome.cells}

    def klass(name: str) -> str:
        return cls.get(name, name).upper()

    conn: dict[str, dict[str, set]] = defaultdict(lambda: defaultdict(set))
    with open(edge_pairs_path, newline="") as f:
        for r in csv.DictReader(f):
            i = index_of.get(int(r["pair_index"]))
            if i is None:
                continue
            conn[klass(r["source"])][klass(r["target"])].add(i)

    out = {
        s: {t: sorted(idxs) for t, idxs in sorted(tgts.items())} for s, tgts in sorted(conn.items())
    }
    return {"pairs": pairs_out, "conn": out}


def neuropeptide_datasets() -> list[dict]:
    """The neuron-graph dataset entries for the neuropeptide range models.

    Each range model is its own viz database (collection np_sr/np_mr/np_lr) holding a single
    dataset, so it is selected via the database dropdown rather than crammed onto the "complete"
    timeline. Opt-in (not a default database) so the predicted network never crowds the observed one.
    """
    return [
        {
            "id": viz_id,
            "name": name,
            "type": collection,
            "time": _NEUROPEPTIDE_VISUAL_TIME,
            "visualTime": _NEUROPEPTIDE_VISUAL_TIME,
            "description": description,
            "datatypes": "np",
        }
        for _kg_id, viz_id, collection, name, description in _NEUROPEPTIDE_DATASETS
    ]


# --- Monoamine projection: the Ripoll-Sánchez 2023 predicted aminergic network ------------------
#
# The companion predicted network to the neuropeptide one. Its own viz database (collection
# "monoamine") with a single dataset, surfaced as an opt-in "ma" edge datatype; off by default so
# it never crowds the observed connectome. Modelled on the adult hermaphrodite (visual_time 50).
#: KG dataset id (RDF curie; no length limit) and the short viz dataset id (varchar(20) in the DB).
_MONOAMINE_DATASET_ID = "ripoll_2023_monoamine"
_MONOAMINE_VIZ_ID = "ripoll_2023_ma"


def monoamine_cells_projection(connectome: object) -> list[dict]:
    """Project the neurons of the monoamine network into the /api/cells shape."""
    endpoints: set[str] = set()
    for conn in connectome.connections:
        if str(conn.connection_type) != "monoaminergic":
            continue
        endpoints.add(_strip(conn.pre, CELL_PREFIX))
        endpoints.add(_strip(conn.post, CELL_PREFIX))
    return [
        {
            "name": c.name,
            "class": c.cell_class or c.name,
            "type": _viz_type(c),
            "neurotransmitter": c.neurotransmitter or "u",
            "embryonic": bool(c.embryonic),
            "inhead": bool(c.in_head),
            "intail": bool(c.in_tail),
        }
        for c in connectome.cells
        if c.name in endpoints
    ]


def monoamine_connections_projection(connectome: object) -> list[dict]:
    """Project the monoamine network into the /api/connections shape (directed, predicted)."""
    grouped: dict[tuple[str, str], int] = defaultdict(int)
    for conn in connectome.connections:
        if str(conn.connection_type) != "monoaminergic":
            continue
        pre = _strip(conn.pre, CELL_PREFIX)
        post = _strip(conn.post, CELL_PREFIX)
        grouped[(pre, post)] += int(conn.weight)
    return [
        {
            "pre": pre,
            "post": post,
            "type": "monoaminergic",
            "annotations": [],
            "synapses": {_MONOAMINE_VIZ_ID: w},
        }
        for (pre, post), w in grouped.items()
    ]


def monoamine_database_cells(connectome: object) -> list[str]:
    """Upper-cased node set (cell names + classes) of the monoamine network, for validNodes."""
    class_of = {c.name: c.cell_class for c in connectome.cells}
    names: set[str] = set()
    for conn in connectome.connections:
        if str(conn.connection_type) != "monoaminergic":
            continue
        for end in (_strip(conn.pre, CELL_PREFIX), _strip(conn.post, CELL_PREFIX)):
            names.add(end.upper())
            cls = class_of.get(end)
            if cls:
                names.add(cls.upper())
    return sorted(names)


def monoamine_datasets() -> list[dict]:
    """The neuron-graph dataset entry for the monoamine connectome (an "ma" datatype database)."""
    return [
        {
            "id": _MONOAMINE_VIZ_ID,
            "name": "Ripoll-Sánchez 2023 (monoamine, predicted)",
            "type": "monoamine",
            "time": _NEUROPEPTIDE_VISUAL_TIME,
            "visualTime": _NEUROPEPTIDE_VISUAL_TIME,
            "description": (
                "Predicted monoamine (aminergic) connectome, Ripoll-Sánchez et al. 2023 "
                "(Neuron 111:3570-3589); receptor pairs from Bentley et al. 2016. "
                "Extrasynaptic — not observed synapses."
            ),
            "datatypes": "ma",
        }
    ]


# --- Male neuropeptide projection: the DERIVED male long-range neuropeptide connectome -----------
#
# The male counterpart to the Ripoll-Sánchez hermaphrodite neuropeptide network (Phase C): derived
# in-build from the Morillo 2025 male expression layer. Surfaced the same way as the hermaphrodite
# network — its OWN viz database (collection "np_male") with a single dataset and the shared "np"
# edge datatype — but male-sexed, so it carries male-specific neurons (CEM, HOA/HOB, spicule/ray)
# that the hermaphrodite network cannot. Opt-in (not a default database); predicted, not observed.
#: KG dataset id (RDF curie; no length limit) and the short viz dataset id (varchar(20) in the DB).
_MORILLO_MALE_NP_DATASET_ID = "morillo_2025_male_neuropeptide_lr"
_MORILLO_MALE_NP_VIZ_ID = "morillo_np_male"


def _male_nt_map(connectome: object) -> dict[str, str]:
    """name -> male neurotransmitter (Wang 2024 male/dimorphic call), as in male_cells_projection."""
    return {
        _strip(str(a.cell), CELL_PREFIX): str(a.neurotransmitter)
        for a in (getattr(connectome, "neurotransmitter_assignments", None) or [])
        if str(a.sex) == "male" and a.neurotransmitter
    }


def male_neuropeptide_cells_projection(connectome: object) -> list[dict]:
    """Project the neurons of the derived male neuropeptide network into the /api/cells shape.

    Endpoints of the ``morillo_2025_male_neuropeptide_lr`` connections (shared + male-specific), with
    the male neurotransmitter call preferred — mirroring :func:`male_cells_projection`.
    """
    male_nt = _male_nt_map(connectome)
    endpoints: set[str] = set()
    for conn in connectome.connections:
        if _strip(str(conn.dataset), DATASET_PREFIX) != _MORILLO_MALE_NP_DATASET_ID:
            continue
        endpoints.add(_strip(conn.pre, CELL_PREFIX))
        endpoints.add(_strip(conn.post, CELL_PREFIX))
    return [
        {
            "name": c.name,
            "class": c.cell_class or c.name,
            "type": _viz_type(c),
            "neurotransmitter": male_nt.get(c.name) or c.neurotransmitter or "u",
            "embryonic": bool(c.embryonic),
            "inhead": bool(c.in_head),
            "intail": bool(c.in_tail),
        }
        for c in connectome.cells
        if c.name in endpoints
    ]


def male_neuropeptide_connections_projection(connectome: object) -> list[dict]:
    """Project the derived male neuropeptide network into the /api/connections shape.

    Directed, predicted edges (no gap-junction merge, no per-synapse ids); one weight under the male
    NP viz dataset id, so the viz shows/selects it like the hermaphrodite range models.
    """
    grouped: dict[tuple[str, str], int] = defaultdict(int)
    for conn in connectome.connections:
        if str(conn.connection_type) != "neuropeptidergic":
            continue
        if _strip(str(conn.dataset), DATASET_PREFIX) != _MORILLO_MALE_NP_DATASET_ID:
            continue
        pre = _strip(conn.pre, CELL_PREFIX)
        post = _strip(conn.post, CELL_PREFIX)
        grouped[(pre, post)] += int(conn.weight)
    return [
        {
            "pre": pre,
            "post": post,
            "type": "neuropeptidergic",
            "annotations": [],
            "synapses": {_MORILLO_MALE_NP_VIZ_ID: w},
        }
        for (pre, post), w in grouped.items()
    ]


def male_neuropeptide_database_cells(connectome: object) -> list[str]:
    """Upper-cased node set (cell names + classes) of the derived male NP network, for validNodes."""
    class_of = {c.name: c.cell_class for c in connectome.cells}
    names: set[str] = set()
    for conn in connectome.connections:
        if str(conn.connection_type) != "neuropeptidergic":
            continue
        if _strip(str(conn.dataset), DATASET_PREFIX) != _MORILLO_MALE_NP_DATASET_ID:
            continue
        for end in (_strip(conn.pre, CELL_PREFIX), _strip(conn.post, CELL_PREFIX)):
            names.add(end.upper())
            cls = class_of.get(end)
            if cls:
                names.add(cls.upper())
    return sorted(names)


def male_neuropeptide_datasets() -> list[dict]:
    """The neuron-graph dataset entry for the derived male NP connectome (an "np" datatype database)."""
    return [
        {
            "id": _MORILLO_MALE_NP_VIZ_ID,
            "name": "Morillo 2025 (male NP, predicted)",
            "type": "np_male",
            "time": _NEUROPEPTIDE_VISUAL_TIME,
            "visualTime": _NEUROPEPTIDE_VISUAL_TIME,
            "description": (
                "Predicted male long-range neuropeptidergic connectome — the male counterpart to the "
                "Ripoll-Sánchez 2023 hermaphrodite network, derived from the Morillo et al. 2025 male "
                "expression atlas (PMID 40682776) over the Beets et al. 2023 NPP-GPCR pairs. Carries "
                "male-specific neurons. Extrasynaptic — not observed synapses."
            ),
            "datatypes": "np",
        }
    ]
