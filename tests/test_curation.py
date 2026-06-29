"""Tests for manual anatomy curation and its precedence over lexical matching."""

from pathlib import Path

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.ingest.neuron_graph import read_cells
from celegans_connectome_kg.match.curation import load_curation
from celegans_connectome_kg.match.matcher import match_cells, summarize
from celegans_connectome_kg.match.wbbt import WBBTIndex

REPO = Path(__file__).resolve().parents[1]
NEURON_GRAPH = REPO / "data" / "neuron-graph"
WBBT_JSON = REPO / "data" / "wbbt" / "wbbt.json"
CURATION = REPO / "data" / "curation" / "anatomy_curation.csv"


def test_curation_covers_the_whole_worklist() -> None:
    curation = load_curation(CURATION)
    assert len(curation) == 114  # 8 ambiguous + 106 unmatched
    assert all(v.startswith("WBbt:") for v in curation.values())


def test_curation_resolves_all_buckets() -> None:
    cells = read_cells(NEURON_GRAPH / "neurons.json")
    index = WBBTIndex.from_obograph(WBBT_JSON)
    counts = summarize(match_cells(cells, index, load_curation(CURATION)))
    assert counts["matched"] == 333
    assert counts["curated"] == 114
    assert counts.get("ambiguous", 0) == 0
    assert counts.get("unmatched", 0) == 0


def test_curation_corrects_false_lexical_match() -> None:
    """M1 lexically collides with the pm1 muscle synonym; curation points to the M1 neuron."""
    cells = read_cells(NEURON_GRAPH / "neurons.json")
    index = WBBTIndex.from_obograph(WBBT_JSON)

    uncurated = {m.cell_name: m for m in match_cells(cells, index)}
    assert uncurated["M1"].status == "ambiguous"  # only the wrong pm1 related-synonym hit

    curated = {m.cell_name: m for m in match_cells(cells, index, load_curation(CURATION))}
    assert curated["M1"].status == "curated"
    assert curated["M1"].wbbt_id == "WBbt:0004488"  # M1 neuron, not pm1 muscle


def test_build_with_curation_grounds_every_cell() -> None:
    _, stats = assemble(NEURON_GRAPH, WBBT_JSON, CURATION)
    assert stats.cells == 447
    assert stats.cells_with_anatomy == 447  # 333 matched + 114 curated
