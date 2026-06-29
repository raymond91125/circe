"""Tests for stub cells minted for class-level connection endpoints."""

from pathlib import Path

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.export.neuron_graph_json import cells_projection
from celegans_connectome_kg.match.curation import load_endpoint_cells

REPO = Path(__file__).resolve().parents[1]
NEURON_GRAPH = REPO / "data" / "neuron-graph"
WBBT_JSON = REPO / "data" / "wbbt" / "wbbt.json"
CURATION = REPO / "data" / "curation" / "anatomy_curation.csv"
ENDPOINT_CELLS = REPO / "data" / "curation" / "connection_endpoint_cells.csv"


def test_load_endpoint_cells() -> None:
    endpoints = load_endpoint_cells(ENDPOINT_CELLS)
    assert len(endpoints) == 11
    by_name = {e.name: e for e in endpoints}
    assert by_name["pm5"].wbbt_id == "WBbt:0003737" and by_name["pm5"].cell_type == "muscle"
    assert by_name["G1"].wbbt_id == "WBbt:0003712" and by_name["G1"].cell_type == "other"
    assert by_name["VAn"].cell_type == "neuron"
    # The 3 reconstruction fragments are intentionally not mapped.
    assert {"Fragment", "NR_fragment", "vncfrag"}.isdisjoint(by_name)


def test_build_mints_stub_cells_and_reduces_unknowns() -> None:
    connectome, stats = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    assert stats.cells == 458  # 447 neuron-graph cells + 11 stubs
    assert stats.cells_with_anatomy == 458  # every cell grounded
    assert stats.unknown_connection_cells == 3  # only the unmapped fragments remain

    by_id = {c.id: c for c in connectome.cells}
    pm5 = by_id["cckg:cell/pm5"]
    assert str(pm5.anatomy) == "WBbt:0003737" and str(pm5.cell_type) == "muscle"
    assert pm5.nemanode_type is None  # stub, not a neuron-graph cell


def test_stub_cells_excluded_from_viz_projection() -> None:
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    names = {c["name"] for c in cells_projection(connectome)}
    assert len(names) == 447  # stubs excluded -> still matches neuron-graph /api/cells
    assert "pm5" not in names and "G1" not in names
