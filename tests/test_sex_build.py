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
        gene_expr_xlsx_path=GENE_EXPR_XLSX,
        gene_map_path=GENE_MAP,
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
    # weight = EM sections; the novel I2 -> pharyngeal-muscle edges are present
    edge = {(strip(c.pre), strip(c.post)): c.weight for c in bh}
    assert edge[("I2L", "pm3VL")] == 133.0
    assert any(strip(c.post).startswith("pm") for c in bh)


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
    genes = connectome.genes
    assert len(genes) == 46 == stats.genes
    assert all(g.id.startswith("WB:WBGene") for g in genes)
    assert {str(g.category) for g in genes} == {
        "metabotropic_receptor",
        "ionotropic_receptor",
        "innexin",
        "neuropeptide",
    }
    assert len(connectome.gene_expressions) == stats.gene_expressions == 309
    assert any(d.id.endswith("cook_2020_pharynx_expression") for d in connectome.datasets)
    assert {str(e.confidence) for e in connectome.gene_expressions} <= {"reported", "putative"}


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
