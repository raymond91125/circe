"""Tests for the sex-aware assemble — Cook 2019 male + hermaphrodite merged [M5]."""

from pathlib import Path

import pytest

from celegans_connectome_kg.build.assemble import assemble

REPO = Path(__file__).resolve().parents[1]
NG = REPO / "data" / "neuron-graph"
WBBT = REPO / "data" / "wbbt" / "wbbt.json"
CUR = REPO / "data" / "curation"
COOK_XLSX = (
    REPO
    / "data"
    / "cook-2019-connectome"
    / ("SI5_connectome_adjacency_matrices_corrected_2020.xlsx")
)
COOK_2020_EDGES = REPO / "data" / "cook-2020-pharynx" / "edges.csv"
GENE_EXPR_XLSX = REPO / "data" / "cook-2020-pharynx" / "SI6_gene_expression.xlsx"
GENE_MAP = REPO / "data" / "cook-2020-pharynx" / "si6_genes.csv"
BHATLA_I2 = REPO / "data" / "bhatla-2015-i2" / "i2_synapses.csv"
DAUER = REPO / "data" / "yim-2024-dauer" / "dauer_connections.csv"
LIFE_STAGE = CUR / "dataset_life_stage.csv"
NEUROTRANSMITTER = REPO / "data" / "wang-neurotransmitter-atlas" / "sex_neurotransmitters.csv"
ATLAS_ONLY = CUR / "atlas_only_cells.csv"


@pytest.fixture(scope="module")
def built():
    connectome, stats = assemble(
        NG,
        WBBT,
        curation_path=CUR / "anatomy_curation.csv",
        endpoint_cells_path=CUR / "connection_endpoint_cells.csv",
        nt_curation_path=CUR / "neurotransmitter_curation.csv",
        cook_xlsx_path=COOK_XLSX,
        cook_aliases_path=CUR / "cook_name_aliases.csv",
        cook_anatomy_path=CUR / "cook_anatomy_curation.csv",
        cook_2020_edges_path=COOK_2020_EDGES,
        bhatla_i2_path=BHATLA_I2,
        dauer_path=DAUER,
        life_stage_path=LIFE_STAGE,
        gene_expr_xlsx_path=GENE_EXPR_XLSX,
        gene_map_path=GENE_MAP,
        neurotransmitter_path=NEUROTRANSMITTER,
        atlas_only_cells_path=ATLAS_ONLY,
        innexin_expr_path=REPO / "data" / "bhattacharya-2019-innexin" / "innexin_expression.csv",
        innexin_gene_map_path=REPO / "data" / "bhattacharya-2019-innexin" / "innexin_genes.csv",
    )
    return connectome, stats


def test_bhatla_i2_dataset(built) -> None:
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    sex = {strip(d.id): str(d.sex) for d in connectome.datasets}
    assert sex.get("bhatla_2015_i2") == "hermaphrodite"
    bh = [c for c in connectome.connections if strip(c.dataset) == "bhatla_2015_i2"]
    assert len(bh) == 26
    assert {strip(c.pre) for c in bh} == {"I2L", "I2R"}
    # weight = EM sections; the heavily-weighted I2 -> pharyngeal-muscle edges are present
    edge = {(strip(c.pre), strip(c.post)): c.weight for c in bh}
    assert edge[("I2L", "pm3VL")] == 133.0
    assert any(strip(c.post).startswith("pm") for c in bh)


def test_kg_added_datasets_excluded_from_herm_projection(built) -> None:
    """The hermaphrodite viz projection is neuron-graph-native only: KG-added datasets
    (Cook 2019/2020, Bhatla 2015, Yim 2024 dauer) must not leak in, or their differing weight
    scales would contaminate the viz's complete/head/tail databases."""
    from celegans_connectome_kg.export.neuron_graph_json import connections_projection

    connectome, _ = built
    datasets = {d for c in connections_projection(connectome) for d in c["synapses"]}
    assert datasets  # projection is non-empty
    assert not any(d.startswith(("cook_", "bhatla_", "yim_")) for d in datasets)
    assert all(d.startswith(("white_1986_", "witvliet_2020_", "randi_funconn_")) for d in datasets)


def test_dauer_dataset(built) -> None:
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    ds = {strip(d.id): d for d in connectome.datasets}
    dauer = ds["yim_2024_dauer"]
    assert str(dauer.sex) == "hermaphrodite"
    assert str(dauer.life_stage) == "dauer"
    conns = [c for c in connectome.connections if strip(c.dataset) == "yim_2024_dauer"]
    assert len(conns) == 2200
    # chemical only (the study did not reconstruct gap junctions); weight = synapse count
    assert all(str(c.connection_type) == "chemical" for c in conns)
    assert int(sum(c.weight for c in conns)) == 6371
    # the excretory duct cell (only non-neuron-graph partner) is a specific cell, not a placeholder
    exc = next(c for c in connectome.cells if c.name == "exc_duct")
    assert exc.unspecified is False and str(exc.anatomy) == "WBbt:0004540"


def test_dataset_life_stage_backfill(built) -> None:
    """Every dataset carries a curated life stage; the Witvliet developmental series and the
    dauer dataset are distinguishable as structured data."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    stage = {
        strip(d.id): (str(d.life_stage) if d.life_stage else None) for d in connectome.datasets
    }
    # every dataset is staged except the non-dauer innexin expression dataset, which spans stages
    assert {k for k, v in stage.items() if v is None} == {"bhattacharya_2019_innexin"}
    assert stage["white_1986_jsh"] == "L4" and stage["white_1986_n2u"] == "adult"
    assert stage["witvliet_2020_1"] == "L1" and stage["witvliet_2020_5"] == "L2"
    assert stage["witvliet_2020_6"] == "L3" and stage["witvliet_2020_7"] == "adult"
    assert stage["cook_2019_male"] == "adult" and stage["yim_2024_dauer"] == "dauer"


def test_placeholder_endpoint_cells_flagged_and_sexless(built) -> None:
    """Class-level endpoint placeholders (the unspecified VA-class 'VAn', pharyngeal muscle
    classes, glands, etc.) are flagged ``unspecified`` and carry no sex-presence, so they can't
    masquerade as sex-specific cells in per-cell analyses."""
    connectome, _ = built
    placeholders = [c for c in connectome.cells if c.unspecified]
    names = {c.name for c in placeholders}
    assert {"VAn", "pm3", "g1"} <= names  # a VA-class neuron endpoint + pharyngeal muscle/gland
    assert all(list(c.sexes) == [] for c in placeholders)  # no sex-presence
    # VAn stays a grounded neuron with its one edge preserved — only the spurious sex tag is gone
    van = next(c for c in connectome.cells if c.name == "VAn")
    assert van.unspecified and str(van.cell_type) == "neuron"
    assert str(van.anatomy) == "WBbt:0005339"  # "VA neuron" (the class)


def test_hermaphrodite_specific_neurons_are_canonical(built) -> None:
    """Only the true hermaphrodite-specific neurons (HSN + VC1-6) are hermaphrodite-only; the
    'VAn' placeholder no longer leaks in now that placeholders are sexless/unspecified."""
    connectome, _ = built
    herm_only = {
        c.name
        for c in connectome.cells
        if str(c.cell_type) == "neuron" and {str(s) for s in c.sexes} == {"hermaphrodite"}
    }
    assert herm_only == {"HSNL", "HSNR", "VC1", "VC2", "VC3", "VC4", "VC5", "VC6"}


def test_neurotransmitter_assignments_per_sex(built) -> None:
    """Wang 2024 (eLife 95402) male atlas: male-specific neurons get a call, and sexually-dimorphic
    sex-shared neurons carry distinct hermaphrodite vs male calls. Cell.neurotransmitter (the
    hermaphrodite/neuron-graph call) is left untouched."""
    connectome, stats = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    na = connectome.neurotransmitter_assignments
    assert len(na) == stats.neurotransmitter_assignments == 119
    by = {(strip(a.cell), str(a.sex)): str(a.neurotransmitter) for a in na}
    # male-specific neurons now have a neurotransmitter (were None on the cell)
    assert by[("CEMDL", "male")] == "a"  # cholinergic
    assert by[("R7AL", "male")] == "d"  # dopaminergic ray neuron
    assert by[("R3BL", "male")] == "s"  # serotonergic ray neuron
    # sexually dimorphic sex-shared neuron: Glu in hermaphrodite, ACh in male (AIM switch)
    assert by[("AIML", "hermaphrodite")] == "ls" and by[("AIML", "male")] == "as"
    assert by[("ADFL", "hermaphrodite")] == "as" and by[("ADFL", "male")] == "ags"  # male +GABA
    # provenance + confidence recorded; Cell.neurotransmitter untouched for male-specific cells
    assert all("95402" in str(a.source) for a in na)
    cemdl = next(c for c in connectome.cells if c.name == "CEMDL")
    assert cemdl.neurotransmitter is None  # the per-sex call lives on the assignment, not the cell


def test_atlas_only_cells_minted_with_neurotransmitter(built) -> None:
    """Neurons in the Wang atlas but absent from the Cook connectome (CP0, DX4, EF4) are minted as
    grounded male neurons so their atlas neurotransmitter attaches; they carry no connections."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    by_name = {c.name: c for c in connectome.cells}
    nt = {
        (strip(a.cell), str(a.sex)): str(a.neurotransmitter)
        for a in connectome.neurotransmitter_assignments
    }
    for name, wbbt, code in [
        ("CP0", "WBbt:0004903", "l"),
        ("DX4", "WBbt:0007845", "a"),
        ("EF4", "WBbt:0007841", "g"),
    ]:
        cell = by_name[name]
        assert str(cell.cell_type) == "neuron" and str(cell.anatomy) == wbbt
        assert {str(s) for s in cell.sexes} == {"male"}
        assert nt[(name, "male")] == code  # atlas call now attaches
    # they are connectivity-free (present in no connection)
    endpoints = {strip(c.pre) for c in connectome.connections} | {
        strip(c.post) for c in connectome.connections
    }
    assert {"CP0", "DX4", "EF4"}.isdisjoint(endpoints)


def test_male_projection_uses_male_neurotransmitter(built) -> None:
    from celegans_connectome_kg.export.neuron_graph_json import male_cells_projection

    connectome, _ = built
    by_name = {c["name"]: c for c in male_cells_projection(connectome)}
    # male-specific neurons carry their Wang-atlas call; dimorphic AIM shows its male call (ACh)
    assert by_name["CEMDL"]["neurotransmitter"] == "a"
    assert by_name["R7AL"]["neurotransmitter"] == "d"
    assert by_name["AIML"]["neurotransmitter"] == "as"  # male call, not the herm "ls"
    # a sex-shared, non-dimorphic neuron falls back to its (shared) neurotransmitter
    assert by_name["AVAL"]["neurotransmitter"] == "a"


def test_datasets_tagged_by_sex(built) -> None:
    connectome, _ = built
    sex = {d.id.split("/")[-1]: str(d.sex) for d in connectome.datasets}
    assert sex["cook_2019_male"] == "male"
    assert sex["cook_2019_hermaphrodite"] == "hermaphrodite"
    assert sex["white_1986_whole"] == "hermaphrodite"  # neuron-graph = hermaphrodite


def test_cell_sex_presence(built) -> None:
    connectome, _ = built
    sexes = {c.name: {str(s) for s in c.sexes} for c in connectome.cells}
    assert sexes["AVAL"] == {"hermaphrodite", "male"}  # shared core neuron
    assert sexes["HSNL"] == {"hermaphrodite"}  # hermaphrodite-specific
    assert sexes["CEMDL"] == {"male"}  # male-specific sensory neuron
    assert sexes["R1AL"] == {"male"}  # male ray neuron


def test_cook_specific_cells_present_and_grounded(built) -> None:
    connectome, _ = built
    by_name = {c.name: c for c in connectome.cells}
    # male-specific cells were minted and grounded to WBbt
    for name in ("CEMDL", "R1AL", "CA9", "ailL", "dglL1"):
        assert name in by_name and by_name[name].anatomy
    # WBBT-ancestry cell typing
    assert str(by_name["CEMDL"].cell_type) == "neuron"
    assert str(by_name["ailL"].cell_type) == "muscle"  # anterior inner longitudinal muscle


def test_male_specific_cell_classes_from_wbbt(built) -> None:
    """Cook-only cells derive a class from their WBbt ``is_a`` parent (bare class token only).

    Bilateral male-specific cells group to their class; serially-repeated neurons and
    pharyngeal endpoints stay classless, matching neuron-graph.
    """
    connectome, _ = built
    cls = {c.name: (str(c.cell_class) if c.cell_class else None) for c in connectome.cells}
    # bilateral pairs group to their WBbt class term
    assert cls["R5AL"] == cls["R5AR"] == "R5A"
    assert cls["CEMDL"] == cls["CEMDR"] == cls["CEMVL"] == cls["CEMVR"] == "CEM"
    assert cls["PCAL"] == "PCA" and cls["MCML"] == "MCM"
    # serial neurons (parent "CA neuron"/"CP neuron") and pharyngeal endpoints keep no class
    assert cls["CA1"] is None and cls["CA9"] is None and cls["CP6"] is None
    assert cls["pm3"] is None and cls["g1"] is None


def test_gene_expression_ingest(built) -> None:
    connectome, stats = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    genes = connectome.genes
    assert len(genes) == stats.genes
    assert all(g.id.startswith("WB:WBGene") for g in genes)
    assert len({g.id for g in genes}) == len(genes)  # genes deduped across sources by WBGene id
    assert {str(g.category) for g in genes} == {
        "metabotropic_receptor",
        "ionotropic_receptor",
        "innexin",
        "neuropeptide",
    }
    # Cook 2020 SI6 contribution is unchanged (309 records) by the added innexin datasets
    cook = [
        e for e in connectome.gene_expressions if strip(e.dataset) == "cook_2020_pharynx_expression"
    ]
    assert len(cook) == 309
    assert len(connectome.gene_expressions) == stats.gene_expressions
    assert any(d.id.endswith("cook_2020_pharynx_expression") for d in connectome.datasets)
    assert {str(e.confidence) for e in connectome.gene_expressions} <= {"reported", "putative"}


def test_innexin_expression_dauer_split(built) -> None:
    """Bhattacharya 2019 Fig 1B innexin expression, split into non-dauer + dauer datasets so the
    dauer plasticity is explicit; class labels expand to member cells; genes reuse Cook's."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    ge = connectome.gene_expressions
    nd = [e for e in ge if strip(e.dataset) == "bhattacharya_2019_innexin"]
    da = [e for e in ge if strip(e.dataset) == "bhattacharya_2019_innexin_dauer"]
    assert len(nd) == 2118 and len(da) == 1885
    ds = {strip(d.id): d for d in connectome.datasets}
    assert str(ds["bhattacharya_2019_innexin_dauer"].life_stage) == "dauer"
    assert ds["bhattacharya_2019_innexin"].life_stage is None  # non-dauer spans stages
    # the 6 innexins not already in Cook SI6 were added (keyed to WBGene)
    innex = {g.symbol for g in connectome.genes if str(g.category) == "innexin"}
    assert {"che-7", "inx-5", "inx-6", "inx-11", "inx-13", "eat-5"} <= innex
    # "both" -> record in both datasets (ADA/inx-1a); "non-dauer only" -> non-dauer only (HSN/inx-1)
    ada_nd = {
        (strip(e.cell), e.isoform)
        for e in nd
        if strip(e.cell) == "ADAL" and str(e.gene) == "WB:WBGene00002123"
    }
    ada_da = {
        (strip(e.cell), e.isoform)
        for e in da
        if strip(e.cell) == "ADAL" and str(e.gene) == "WB:WBGene00002123"
    }
    assert ("ADAL", "a") in ada_nd and ("ADAL", "a") in ada_da  # both
    assert any(str(e.cell).endswith("/HSNL") and str(e.gene) == "WB:WBGene00002123" for e in nd)
    assert not any(str(e.cell).endswith("/HSNL") and str(e.gene) == "WB:WBGene00002123" for e in da)


def test_gene_expression_per_cell_and_isoform(built) -> None:
    connectome, _ = built
    ge = connectome.gene_expressions

    def genes_of(cell):
        return {str(e.gene) for e in ge if str(e.cell).endswith("/" + cell)}

    # SI6's I1L/R class row expands to both member cells; gar-2 expressed, gar-1 not.
    for c in ("I1L", "I1R"):
        assert "WB:WBGene00001518" in genes_of(c)  # gar-2
        assert "WB:WBGene00001517" not in genes_of(c)  # gar-1 (blank in SI6)
    # inx-1 isoforms a and b are kept distinct (transcript qualifier), same persistent gene id.
    inx1 = [e for e in ge if str(e.cell).endswith("/I1L") and str(e.gene) == "WB:WBGene00002123"]
    assert {str(e.isoform) for e in inx1} == {"a", "b"}


def test_build_stats_partition(built) -> None:
    _, stats = built
    assert stats.datasets_by_sex["male"] == 1  # cook_2019_male
    assert stats.datasets_by_sex["hermaphrodite"] >= 15
    # most cells are shared; a substantial male-only and herm-only tail exist
    assert stats.cells_by_sex["hermaphrodite+male"] > 400
    assert stats.cells_by_sex.get("male", 0) > 100


def test_male_viz_projection(built) -> None:
    from celegans_connectome_kg.export.neuron_graph_json import (
        male_cells_projection,
        male_connections_projection,
        male_dataset,
    )

    connectome, _ = built
    cells = male_cells_projection(connectome)
    by_name = {c["name"]: c for c in cells}
    # shared + male-specific present; hermaphrodite-only excluded
    assert {"AVAL", "CEMDL", "R1AL"} <= set(by_name)
    assert "HSNL" not in by_name
    # type synthesis for male-specific cells (subtype absent -> neuron=i, muscle=b)
    assert by_name["CEMDL"]["type"] == "i" and by_name["ailL"]["type"] == "b"
    # shared cells keep their real NemaNode type + class
    assert by_name["AVAL"]["type"] == "i" and by_name["AVAL"]["class"] == "AVA"
    # male-specific bilateral pairs group by their WBbt class term (R5AL/R5AR -> R5A)
    assert by_name["R5AL"]["class"] == "R5A" and by_name["R5AR"]["class"] == "R5A"

    conns = male_connections_projection(connectome)
    assert conns and all("cook_2019_male" in c["synapses"] for c in conns)

    ds = male_dataset()
    assert ds["type"] == "male" and ds["datatypes"] == "cs,gj"


def test_pharynx_viz_projection(built) -> None:
    from celegans_connectome_kg.export.neuron_graph_json import (
        pharynx_cells_projection,
        pharynx_connections_projection,
        pharynx_dataset,
    )

    connectome, _ = built
    cells = pharynx_cells_projection(connectome)
    by_name = {c["name"]: c for c in cells}
    # pharyngeal cells present, with class grouping (M3L -> M3); non-pharyngeal excluded
    assert {"M3L", "M3R", "I1L", "MCL"} <= set(by_name)
    assert by_name["M3L"]["class"] == "M3" and by_name["M3R"]["class"] == "M3"
    assert "AVAL" not in by_name  # somatic neuron, not a pharynx-dataset endpoint

    conns = pharynx_connections_projection(connectome)
    assert conns and all("cook_2020_pharynx" in c["synapses"] for c in conns)

    ds = pharynx_dataset()
    assert ds["type"] == "pharynx" and ds["datatypes"] == "cs,gj"


def test_dauer_viz_projection(built) -> None:
    from celegans_connectome_kg.export.neuron_graph_json import (
        dauer_cells_projection,
        dauer_connections_projection,
        dauer_dataset,
    )

    connectome, _ = built
    cells = dauer_cells_projection(connectome)
    by_name = {c["name"]: c for c in cells}
    assert len(cells) == 221  # 181 neurons + muscle/other partners
    assert "exc_duct" in by_name  # the curated excretory-duct endpoint is projected

    conns = dauer_connections_projection(connectome)
    assert len(conns) == 2200
    # chemical only (no gap junctions reconstructed), weighted by synapse count
    assert all(c["type"] == "chemical" for c in conns)
    assert all("yim_2024_dauer" in c["synapses"] for c in conns)

    ds = dauer_dataset()
    # folded into the "head" life-stage series, positioned in the L3 region; chemical only
    assert ds["type"] == "head" and ds["datatypes"] == "cs"
    assert 25 <= ds["visualTime"] <= 34
