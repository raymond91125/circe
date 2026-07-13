"""Phase 3 build tests: cell-type classification and full-snapshot assembly."""

from pathlib import Path

import pytest

from celegans_connectome_kg.build.assemble import assemble, classify_cell_type

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
NEURON_GRAPH = DATA_DIR / "neuron-graph"
WBBT_JSON = DATA_DIR / "wbbt" / "wbbt.json"


@pytest.mark.parametrize(
    "code,expected",
    [
        ("b", "muscle"),  # body wall + pharyngeal + defecation muscles
        ("s", "neuron"),
        ("i", "neuron"),
        ("m", "neuron"),
        ("sm", "neuron"),
        ("imn", "neuron"),
        ("u", "other"),  # marginal cells (mc*)
        ("", "other"),
        (None, "other"),
    ],
)
def test_classify_cell_type(code: str | None, expected: str) -> None:
    assert classify_cell_type(code) == expected


@pytest.fixture(scope="module")
def built() -> tuple[object, object]:
    return assemble(NEURON_GRAPH, WBBT_JSON)


def test_build_stats(built: tuple[object, object]) -> None:
    _, stats = built
    assert stats.cells == 447
    assert stats.datasets == 14  # 15 neuron-graph datasets minus de-duplicated randi_funconn_wildcp
    # Source records aggregate by (dataset,pre,post,type); duplicate listings are summed (mirrors
    # neuron-graph populate). randi_funconn_wildcp is dropped as a byte-identical duplicate of
    # randi_funconn_wildty (assemble._REDUNDANT_NG_DATASETS), leaving 23895 connections.
    assert stats.connections == 23895
    assert stats.cells_with_anatomy == 333  # equals the matched count from Phase 2
    assert stats.connections_by_type == {
        "chemical": 16945,
        "gap_junction": 5467,
        "functional": 1483,  # 2614 minus randi_funconn_wildcp's 1131 (functional-only duplicate)
    }
    # Distinct class-level / fragment endpoint names absent from neurons.json.
    assert stats.unknown_connection_cells == 14


def test_cells_have_curie_ids_and_typed(built: tuple[object, object]) -> None:
    connectome, _ = built
    adal = next(c for c in connectome.cells if c.id == "cckg:cell/ADAL")
    assert str(adal.cell_type) == "neuron"
    assert str(adal.anatomy) == "WBbt:0004013"
    bwm = next(c for c in connectome.cells if c.id == "cckg:cell/BWM-DL01")
    assert str(bwm.cell_type) == "muscle"
    assert bwm.anatomy is None  # unmatched -> no anatomy


def test_connections_reference_cells_and_carry_evidence(built: tuple[object, object]) -> None:
    connectome, _ = built
    conn = connectome.connections[0]
    assert conn.id.startswith("cckg:conn/")
    assert conn.pre.startswith("cckg:cell/") and conn.post.startswith("cckg:cell/")
    assert str(conn.dataset).startswith("cckg:dataset/")
    assert len(conn.evidence) == 1
    assert str(conn.evidence[0].evidence_type) == "aggregated_weight"


def test_connection_ids_are_unique(built: tuple[object, object]) -> None:
    connectome, _ = built
    ids = [c.id for c in connectome.connections]
    assert len(ids) == len(set(ids))
