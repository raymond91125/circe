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


def test_build_stats_partition(built) -> None:
    _, stats = built
    assert stats.datasets_by_sex["male"] == 1  # cook_2019_male
    assert stats.datasets_by_sex["hermaphrodite"] >= 15
    # most cells are shared; a substantial male-only and herm-only tail exist
    assert stats.cells_by_sex["hermaphrodite+male"] > 400
    assert stats.cells_by_sex.get("male", 0) > 100
