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
    )
    return connectome, stats


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
