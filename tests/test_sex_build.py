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
        neuropeptide_dir=REPO / "data" / "ripoll-2023-neuropeptide",
        neuropeptide_pairs_path=REPO
        / "data"
        / "ripoll-2023-neuropeptide"
        / "mechanistic"
        / "npp_gpcr_pairs.csv",
        neuropeptide_genes_path=REPO
        / "data"
        / "ripoll-2023-neuropeptide"
        / "mechanistic"
        / "np_genes.csv",
        neuropeptide_expression_path=REPO
        / "data"
        / "ripoll-2023-neuropeptide"
        / "mechanistic"
        / "np_gene_expression.csv",
        monoamine_network_path=REPO / "data" / "ripoll-2023-monoamine" / "monoamine_network.csv",
        monoamine_pairs_path=REPO / "data" / "ripoll-2023-monoamine" / "monoamine_pairs.csv",
        monoamine_genes_path=REPO
        / "data"
        / "ripoll-2023-monoamine"
        / "mechanistic"
        / "monoamine_receptor_genes.csv",
        monoamine_expression_path=REPO
        / "data"
        / "ripoll-2023-monoamine"
        / "mechanistic"
        / "monoamine_receptor_expression.csv",
        monoamine_ionotropic_genes_path=REPO
        / "data"
        / "cengen-monoamine-ionotropic"
        / "monoamine_ionotropic_genes.csv",
        monoamine_ionotropic_expression_path=REPO
        / "data"
        / "cengen-monoamine-ionotropic"
        / "monoamine_ionotropic_expression.csv",
        morillo_genes_path=REPO / "data" / "morillo-2025-male" / "morillo_genes.csv",
        morillo_herm_expression_path=REPO
        / "data"
        / "morillo-2025-male"
        / "morillo_herm_expression.csv",
        morillo_male_expression_path=REPO
        / "data"
        / "morillo-2025-male"
        / "morillo_male_expression.csv",
        morillo_receptor_genes_path=REPO
        / "data"
        / "morillo-2025-male"
        / "morillo_receptor_genes.csv",
        morillo_herm_receptor_expression_path=REPO
        / "data"
        / "morillo-2025-male"
        / "morillo_herm_receptor_expression.csv",
        morillo_male_receptor_expression_path=REPO
        / "data"
        / "morillo-2025-male"
        / "morillo_male_receptor_expression.csv",
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
    assert not any(d.startswith(("cook_", "bhatla_", "yim_", "ripoll_")) for d in datasets)
    assert all(d.startswith(("white_1986_", "witvliet_2020_", "randi_funconn_")) for d in datasets)


def test_neuropeptide_connectome(built) -> None:
    """Ripoll-Sanchez 2023 predicted neuropeptidergic connectome: three range-model datasets with a
    new connection type; directed, weighted by NPP-GPCR pathway count; over existing cells."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    np = [c for c in connectome.connections if str(c.connection_type) == "neuropeptidergic"]
    from collections import Counter

    by_ds = Counter(strip(c.dataset) for c in np)
    assert by_ds == {
        "ripoll_2023_neuropeptide_sr": 31417,
        "ripoll_2023_neuropeptide_mr": 40425,
        "ripoll_2023_neuropeptide_lr": 53558,
    }
    assert 1 <= min(c.weight for c in np) and max(c.weight for c in np) == 20
    # datasets tagged hermaphrodite / adult; zero-padded matrix names normalized (DA01 -> DA1)
    ds = {strip(d.id): d for d in connectome.datasets}
    for did in by_ds:
        assert str(ds[did].sex) == "hermaphrodite" and str(ds[did].life_stage) == "adult"
    cells = {strip(c.pre) for c in np} | {strip(c.post) for c in np}
    assert "DA1" in cells and {"DA01", "VD09", "VC06"}.isdisjoint(cells)  # padded names normalized


def test_neuropeptide_in_viz_kg_connections(built) -> None:
    """The cell-info KG connectivity map is comprehensive ("all partners across every dataset"),
    so the predicted neuropeptide network is included, under directional npo/npi relations."""
    from celegans_connectome_kg.export.neuron_graph_json import kg_connections_map

    connectome, _ = built
    m = kg_connections_map(connectome)
    assert any(d.startswith("ripoll_") for d in m["datasets"])
    # At least one class carries a predicted-neuropeptide out or in relation.
    assert any("npo" in rels or "npi" in rels for rels in m["conn"].values())


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
        "neuropeptide_receptor",
        "monoamine_receptor",
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
    # cook_2019_male + morillo_2025_male (ligands) + morillo_2025_male_receptors
    assert stats.datasets_by_sex["male"] == 3
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


def test_neuropeptide_viz_projection(built) -> None:
    from celegans_connectome_kg.export.neuron_graph_json import (
        neuropeptide_cells_projection,
        neuropeptide_connections_projection,
        neuropeptide_database_cells,
        neuropeptide_datasets,
    )

    connectome, _ = built

    # Each range model is its own viz database (collection np_sr/np_mr/np_lr), one "np" dataset each.
    datasets = neuropeptide_datasets()
    ids = [d["id"] for d in datasets]
    assert ids == ["ripoll_2023_np_sr", "ripoll_2023_np_mr", "ripoll_2023_np_lr"]
    assert all(len(d["id"]) <= 20 for d in datasets)  # datasets.id is varchar(20)
    assert all(len(d["name"]) <= 50 for d in datasets)  # datasets.name is varchar(50)
    assert [d["type"] for d in datasets] == ["np_sr", "np_mr", "np_lr"]  # own collections
    assert all(len(d["type"]) <= 20 for d in datasets)  # datasets.collection is varchar(20)
    assert all(d["datatypes"] == "np" for d in datasets)
    assert all("predicted" in d["description"].lower() for d in datasets)

    # All 302 hermaphrodite neurons of the network are projected as cells.
    cells = neuropeptide_cells_projection(connectome)
    assert len(cells) == 302

    # The database node set (upper-cased names + classes) covers those neurons.
    db_nodes = neuropeptide_database_cells(connectome)
    assert db_nodes == sorted(db_nodes) and len(db_nodes) == len(set(db_nodes))
    assert "AVAL" in db_nodes and "AVA" in db_nodes  # a cell and its class

    # Directed neuropeptidergic edges, keyed by the short viz dataset ids; union across ranges.
    conns = neuropeptide_connections_projection(connectome)
    assert len(conns) == 53558  # long-range is the superset (union of all edges)
    assert all(c["type"] == "neuropeptidergic" for c in conns)
    assert all(c["annotations"] == [] for c in conns)
    per_dataset = {i: 0 for i in ids}
    for c in conns:
        for dset in c["synapses"]:
            assert dset in per_dataset  # only short viz ids leak into the projection
            per_dataset[dset] += 1
    assert per_dataset == {
        "ripoll_2023_np_sr": 31417,
        "ripoll_2023_np_mr": 40425,
        "ripoll_2023_np_lr": 53558,
    }


def test_monoamine_connectome(built) -> None:
    """Ripoll-Sanchez 2023 predicted monoamine (aminergic) connectome: a single hermaphrodite/adult
    dataset with a new monoaminergic connection type; directed, weight = monoamine-receptor pathway
    count; 19 canonical aminergic source neurons."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731
    ma = [c for c in connectome.connections if str(c.connection_type) == "monoaminergic"]
    assert len(ma) == 2881
    assert {int(c.weight) for c in ma} == {1, 2, 3}
    assert all(strip(c.dataset) == "ripoll_2023_monoamine" for c in ma)
    ds = {strip(d.id): d for d in connectome.datasets}["ripoll_2023_monoamine"]
    assert str(ds.sex) == "hermaphrodite"
    # the 19 sources are the canonical aminergic neurons (dopamine/serotonin/tyramine/octopamine)
    sources = {strip(c.pre) for c in ma}
    assert len(sources) == 19
    assert {"ADEL", "NSML", "RIML", "RICL", "CEPDL", "PDEL"} <= sources


def test_monoamine_receptor_pairs_ingested(built) -> None:
    """The 14 monoamine-receptor pairs are ingested (Bentley et al. 2016), 4 monoamines."""
    connectome, _ = built
    pairs = connectome.monoamine_receptor_pairs
    assert len(pairs) == 14
    assert {p.monoamine for p in pairs} == {"Ser", "DA", "Oct", "Tyr"}
    assert all(p.receptor for p in pairs)
    # serotonin has receptors ser-1/4/5/7
    ser = {p.receptor for p in pairs if p.monoamine == "Ser"}
    assert ser == {"ser-1", "ser-4", "ser-5", "ser-7"}


def test_monoamine_viz_projection(built) -> None:
    from celegans_connectome_kg.export.neuron_graph_json import (
        monoamine_connections_projection,
        monoamine_database_cells,
        monoamine_datasets,
    )

    connectome, _ = built
    conns = monoamine_connections_projection(connectome)
    assert len(conns) == 2881
    assert all(c["type"] == "monoaminergic" for c in conns)
    assert all("ripoll_2023_ma" in c["synapses"] for c in conns)  # short viz id (varchar(20))
    ds = monoamine_datasets()
    assert len(ds) == 1
    assert ds[0]["type"] == "monoamine" and ds[0]["datatypes"] == "ma"
    assert len(ds[0]["id"]) <= 20 and len(ds[0]["name"]) <= 50
    db = monoamine_database_cells(connectome)
    assert "AVAL" in db or "ADEL" in db


def test_monoamine_receptor_expression_ingested(built) -> None:
    """The mechanistic layer resolves the 14 receptors to WBGene and ingests their per-neuron CeNGEN
    expression under its own dataset; every pair carries a resolved receptor_gene link."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731

    ma_expr = [
        e
        for e in connectome.gene_expressions
        if strip(e.dataset) == "ripoll_2023_monoamine_expression"
    ]
    assert len(ma_expr) == 787
    assert all(str(e.confidence) == "reported" for e in ma_expr)

    # All 14 pairs link to a resolved receptor Gene node (WBGene), the stable expression join key.
    pairs = connectome.monoamine_receptor_pairs
    assert len(pairs) == 14
    gene_ids = {g.id for g in connectome.genes}
    assert all(p.receptor_gene and p.receptor_gene in gene_ids for p in pairs)

    ds = {strip(d.id): d for d in connectome.datasets}["ripoll_2023_monoamine_expression"]
    assert str(ds.sex) == "hermaphrodite" and str(ds.life_stage) == "adult"


def test_monoamine_ionotropic_receptors_ingested(built) -> None:
    """The eight ionotropic (ligand-gated) monoamine receptors -- the amine-gated Cys-loop channels
    mod-1/lgc-40/50/52/53/54/55/56 -- enter the KG as expression-only genes: category
    ionotropic_receptor, per-neuron CeNGEN threshold-4 expression under their own dataset, and --
    unlike the GPCR receptors -- no monoaminergic edges/pairs. lgc-51 (accessory subunit, not itself
    amine-gated) is excluded."""
    from collections import Counter

    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731

    # symbol -> expected (WBGene id, expressing-neuron count)
    expected = {
        "mod-1": ("WB:WBGene00003386", 28),
        "lgc-40": ("WB:WBGene00020767", 29),
        "lgc-53": ("WB:WBGene00020657", 68),
        "lgc-55": ("WB:WBGene00013746", 41),
        "lgc-50": ("WB:WBGene00020605", 60),
        "lgc-52": ("WB:WBGene00013517", 14),
        "lgc-54": ("WB:WBGene00020528", 17),
        "lgc-56": ("WB:WBGene00001588", 41),
    }
    genes = {g.symbol: g for g in connectome.genes}
    assert set(expected) <= set(genes)
    assert all(genes[s].id == wb for s, (wb, _) in expected.items())
    assert all(str(genes[s].category) == "ionotropic_receptor" for s in expected)
    assert "lgc-51" not in genes  # accessory subunit, deliberately excluded

    mi_expr = [
        e
        for e in connectome.gene_expressions
        if strip(e.dataset) == "cengen_2021_monoamine_ionotropic_expression"
    ]
    assert len(mi_expr) == sum(n for _, n in expected.values()) == 298
    by_gene = Counter(strip(e.gene) for e in mi_expr)
    assert all(by_gene[wb] == n for _, (wb, n) in expected.items())
    assert all(str(e.confidence) == "reported" for e in mi_expr)

    ds = {strip(d.id): d for d in connectome.datasets}[
        "cengen_2021_monoamine_ionotropic_expression"
    ]
    assert str(ds.sex) == "hermaphrodite" and str(ds.life_stage) == "adult"

    # Ionotropic receptors are not GPCRs: they must not appear in the monoamine (GPCR) pairs.
    assert not any(p.receptor in expected for p in connectome.monoamine_receptor_pairs)


def test_monoamine_pairs_map_matches_weight(built) -> None:
    """The reconstructed per-edge attribution reproduces the published monoamine weights exactly
    (0 mismatches), and the class-level viz map surfaces the mediating receptors by symbol."""
    import csv
    from collections import defaultdict

    from celegans_connectome_kg.export.neuron_graph_json import monoamine_pairs_map

    connectome, _ = built
    edge_pairs = REPO / "data" / "ripoll-2023-monoamine" / "mechanistic" / "edge_pairs.csv"

    # Cell-level check: for every edge, #mediating pairs == monoaminergic weight (0 mismatches).
    per_edge = defaultdict(int)
    with open(edge_pairs, newline="") as f:
        for r in csv.DictReader(f):
            per_edge[(r["source"], r["target"])] += 1
    assert sum(per_edge.values()) == 4127  # = sum of all 2,881 edge weights

    published = {}
    for c in connectome.connections:
        if str(c.connection_type) != "monoaminergic":
            continue
        pre = str(c.pre).rsplit("/", 1)[-1]
        post = str(c.post).rsplit("/", 1)[-1]
        published[(pre, post)] = int(c.weight)
    assert dict(per_edge) == published  # exact decomposition of the published weights
    assert len(published) == 2881

    # Class-level viz map: 14 pairs, and NSM->AIZ is mediated by all three serotonin receptors.
    m = monoamine_pairs_map(connectome, edge_pairs)
    assert len(m["pairs"]) == 14
    nsm_aiz = m["conn"]["NSM"]["AIZ"]
    recs = {m["pairs"][i][1] for i in nsm_aiz}
    assert recs == {"ser-1", "ser-4", "ser-5"}
    assert all(m["pairs"][i][0] == "Ser" for i in nsm_aiz)


def test_neuropeptide_receptor_pairs_ingested(built) -> None:
    """The 92 deorphanized NPP-GPCR pairs are ingested as KG entities, with EC50 + GPCR class."""
    connectome, _ = built
    pairs = connectome.neuropeptide_receptor_pairs
    assert len(pairs) == 92
    assert all(p.ligand and p.gpcr for p in pairs)
    assert all(p.ec50_nm is not None for p in pairs)  # all 92 matched an EC50 in Table S5
    # The first pair (canonical order) is nlp-40 -> aex-2.
    first = min(pairs, key=lambda p: int(str(p.id).rsplit("/", 1)[-1]))
    assert (first.ligand, first.gpcr) == ("nlp-40", "aex-2")


def test_neuropeptide_expression_ingested(built) -> None:
    """Per-neuron NPP/GPCR expression is ingested as GeneExpression under its own dataset, with the
    genes categorized and the 12 genes shared with existing datasets not duplicated."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731

    np_expr = [
        e for e in connectome.gene_expressions if strip(e.dataset) == "ripoll_2023_expression"
    ]
    assert len(np_expr) == 3514
    assert all(str(e.confidence) == "reported" for e in np_expr)

    # 49 NPP + 51 GPCR genes, categorized; genes keyed by WBGene (shared, not duplicated).
    npp = [g for g in connectome.genes if str(g.category) == "neuropeptide"]
    npr = [g for g in connectome.genes if str(g.category) == "neuropeptide_receptor"]
    assert len(npr) == 51 and len(npp) >= 49
    assert len({g.id for g in connectome.genes}) == len(connectome.genes)  # no dup gene ids

    # The expression dataset is a hermaphrodite / adult CeNGEN source.
    ds = {strip(d.id): d for d in connectome.datasets}["ripoll_2023_expression"]
    assert str(ds.sex) == "hermaphrodite" and str(ds.life_stage) == "adult"


def test_neuropeptide_network_derivable_from_expression(built) -> None:
    """The KG can DERIVE the neuropeptide network from first principles: joining per-neuron
    expression with the NPP-GPCR pairs (source expresses ligand, target expresses cognate GPCR)
    reproduces the published long-range connectome exactly."""
    from collections import defaultdict

    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731

    # expression keyed by WBGene (the stable join key — some ligand/gpcr symbols are synonyms of a
    # different canonical gene symbol, e.g. nlp-54 == trh-1, so a symbol join would be wrong).
    expr = defaultdict(set)
    for e in connectome.gene_expressions:
        if strip(e.dataset) == "ripoll_2023_expression":
            expr[e.gene].add(strip(e.cell))

    # derive edges by joining pairs to expression through the pairs' WBGene gene links
    derived = defaultdict(int)
    for p in connectome.neuropeptide_receptor_pairs:
        for s in expr.get(p.ligand_gene, ()):
            for t in expr.get(p.gpcr_gene, ()):
                derived[(s, t)] += 1

    # published long-range network from the KG connections
    published = {}
    for c in connectome.connections:
        if strip(c.dataset) == "ripoll_2023_neuropeptide_lr":
            published[(strip(c.pre), strip(c.post))] = int(c.weight)

    assert dict(derived) == published  # exact first-principles reconstruction, 53,558 edges
    assert len(published) == 53558


def test_morillo_male_expression_ingested(built) -> None:
    """Morillo et al. 2025 (PMID 40682776) male snRNA-seq atlas, Phase A: sex-shared FLP/NLP/INS
    neuropeptide-ligand expression as two sex-tagged GeneExpression datasets (male + herm),
    binarized at a CeNGEN-calibrated cutoff. Class columns expand to member cells; merged clusters
    and non-neuronal columns are dropped; male-specific neurons are included."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731

    # Two sex-tagged adult datasets.
    ds = {strip(d.id): d for d in connectome.datasets}
    assert str(ds["morillo_2025_male"].sex) == "male"
    assert str(ds["morillo_2025_hermaphrodite"].sex) == "hermaphrodite"
    assert all(
        str(ds[k].life_stage) == "adult"
        for k in ("morillo_2025_male", "morillo_2025_hermaphrodite")
    )

    herm = [
        e for e in connectome.gene_expressions if strip(e.dataset) == "morillo_2025_hermaphrodite"
    ]
    male = [e for e in connectome.gene_expressions if strip(e.dataset) == "morillo_2025_male"]
    assert len(herm) == 1556 and len(male) == 809
    assert all(str(e.confidence) == "reported" for e in herm + male)

    # 133 curated neuropeptide-ligand genes, keyed by WBGene so they share with the Ripoll
    # neuropeptide gene nodes (no duplication). nlp-73/nlp-77 resolve to their canonical
    # WormBase synonyms ctdh-1 (WBGene00013964) / srlf-1 (WBGene00015514).
    mo_genes = {strip(e.gene) for e in herm + male}
    assert len(mo_genes) <= 133
    # nlp-73 (ctdh-1), nlp-77 (srlf-1)
    assert {"WB:WBGene00013964", "WB:WBGene00015514"} <= mo_genes
    assert len({g.id for g in connectome.genes}) == len(connectome.genes)

    herm_cells = {strip(e.cell) for e in herm}
    male_cells = {strip(e.cell) for e in male}
    # class columns fan out to member cells (herm "ADE" -> ADEL + ADER)
    assert {"ADEL", "ADER"} <= herm_cells
    # male-specific neurons are included (CEM class + individually-named male neurons)
    assert {"CEMDL", "CEMDR", "CEMVL", "CEMVR", "HOA", "HOB"} <= male_cells
    male_only = {c.name for c in connectome.cells if {str(s) for s in c.sexes} == {"male"}}
    assert male_cells & male_only  # genuine male-specific coverage
    # merged clusters ("PHA/PHB/ASH", "PHA/PHB") are dropped, not fanned out
    assert "PHAL" not in herm_cells and "PHAL" not in male_cells


def test_morillo_receptor_expression_ingested(built) -> None:
    """Morillo 2025 Phase B: neuropeptide-receptor (NPP-GPCR) expression from the genome-wide atlas
    (mmc12, cluster-level), both sexes, binarized at a receptor-specific CeNGEN-calibrated cutoff.
    The receptor side the male neuropeptide connectome needs (ligands came from mmc8 in Phase A)."""
    connectome, _ = built
    strip = lambda s: str(s).split("/")[-1]  # noqa: E731

    ds = {strip(d.id): d for d in connectome.datasets}
    assert str(ds["morillo_2025_male_receptors"].sex) == "male"
    assert str(ds["morillo_2025_hermaphrodite_receptors"].sex) == "hermaphrodite"
    assert all(
        str(ds[k].life_stage) == "adult"
        for k in ("morillo_2025_male_receptors", "morillo_2025_hermaphrodite_receptors")
    )

    herm = [
        e
        for e in connectome.gene_expressions
        if strip(e.dataset) == "morillo_2025_hermaphrodite_receptors"
    ]
    male = [
        e for e in connectome.gene_expressions if strip(e.dataset) == "morillo_2025_male_receptors"
    ]
    assert len(herm) == 1664 and len(male) == 1603
    assert all(str(e.confidence) == "reported" for e in herm + male)

    # 45 NPP-GPCR receptors, category neuropeptide_receptor, keyed by WBGene so they share with the
    # Ripoll receptor gene nodes (no duplication). 6 (npr-37/39/40/41/42/43) are below mmc12's 0.5%
    # inclusion threshold and absent.
    recv = [g for g in connectome.genes if str(g.category) == "neuropeptide_receptor"]
    assert len({g.id for g in recv}) == len(recv)  # receptor gene nodes deduped by WBGene
    assert len({strip(e.gene) for e in herm + male}) == 45

    herm_cells = {strip(e.cell) for e in herm}
    male_cells = {strip(e.cell) for e in male}
    # cross-sex columns (mmc12 lists all clusters in every sheet) are resolved per sex, not errored:
    # hermaphrodite-specific VCn/HSN only in herm; male-specific CEM only in male.
    male_only = {c.name for c in connectome.cells if {str(s) for s in c.sexes} == {"male"}}
    assert male_cells & male_only  # male-specific receptor coverage
    assert not (herm_cells & male_only)  # no male-specific cells leak into the herm dataset
    assert "HSNL" not in male_cells  # herm-specific neuron absent from the male dataset


def test_neuropeptide_pairs_map_cell_level_matches_weight(built) -> None:
    """Per-edge pair attribution reproduces the published edge weights exactly (validation), and
    the class-level viz map surfaces the mediating pairs by gene symbol."""
    from celegans_connectome_kg.export.neuron_graph_json import neuropeptide_pairs_map

    connectome, _ = built
    edge_pairs = REPO / "data" / "ripoll-2023-neuropeptide" / "mechanistic" / "edge_pairs.csv"

    # Cell-level check: for every edge, #mediating pairs == long-range weight (0 mismatches).
    import csv
    from collections import defaultdict

    per_edge = defaultdict(int)
    with open(edge_pairs, newline="") as f:
        for r in csv.DictReader(f):
            per_edge[(r["source"], r["target"])] += 1
    assert sum(per_edge.values()) == 145834

    lr = {}
    for c in connectome.connections:
        if str(c.connection_type) != "neuropeptidergic":
            continue
        if str(c.dataset).rsplit("/", 1)[-1] != "ripoll_2023_neuropeptide_lr":
            continue
        pre = str(c.pre).rsplit("/", 1)[-1]
        post = str(c.post).rsplit("/", 1)[-1]
        lr[(pre, post)] = int(c.weight)
    assert per_edge == lr  # exact decomposition of the published long-range weights

    # Class-level viz map: 92 pairs, and AVA->ADA is mediated by flp-18 -> npr-4.
    m = neuropeptide_pairs_map(connectome, edge_pairs)
    assert len(m["pairs"]) == 92
    ava_ada = m["conn"]["AVA"]["ADA"]
    sym = {(m["pairs"][i][0], m["pairs"][i][1]) for i in ava_ada}
    assert ("flp-18", "npr-4") in sym
