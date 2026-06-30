"""Tests for the cell-class → WBbt map that drives the viz cell-info WormBase link."""

from pathlib import Path

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.export.neuron_graph_json import anatomy_terms_map
from celegans_connectome_kg.match.curation import load_class_curation
from celegans_connectome_kg.match.wbbt import WBBTIndex

REPO = Path(__file__).resolve().parents[1]
NEURON_GRAPH = REPO / "data" / "neuron-graph"
WBBT_JSON = REPO / "data" / "wbbt" / "wbbt.json"
CURATION = REPO / "data" / "curation" / "anatomy_curation.csv"
ENDPOINT_CELLS = REPO / "data" / "curation" / "connection_endpoint_cells.csv"
CLASS_CURATION = REPO / "data" / "curation" / "class_anatomy_curation.csv"


def test_load_class_curation() -> None:
    cc = load_class_curation(CLASS_CURATION)
    assert len(cc) == 28
    assert cc["VAn"] == "WBbt:0005339"  # VA neuron class, not a VA muscle/placeholder
    assert cc["BWM01"] == "WBbt:0006804"  # positional BWM -> generic body wall muscle cell


def test_anatomy_terms_cover_classes_cells_and_stubs() -> None:
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    index = WBBTIndex.from_obograph(WBBT_JSON)
    terms = anatomy_terms_map(connectome, index, load_class_curation(CLASS_CURATION))

    # every cell class is covered (keys upper-cased to match cellClass())
    classes = {c.cell_class.upper() for c in connectome.cells if c.cell_class}
    assert classes <= set(terms)
    # every grounded cell name is covered too (handles individual display + non-neuron nodes)
    cell_names = {c.name.upper() for c in connectome.cells if c.anatomy}
    assert cell_names <= set(terms)
    assert all(k == k.upper() for k in terms)
    assert all(v.startswith("WBbt:") for v in terms.values())


def test_anatomy_terms_resolution_sources() -> None:
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    index = WBBTIndex.from_obograph(WBBT_JSON)
    terms = anatomy_terms_map(connectome, index, load_class_curation(CLASS_CURATION))

    assert terms["AVA"] == "WBbt:0005842"  # strong class-name match (the reported example)
    assert terms["M1"] == "WBbt:0004488"  # corrected cell curation (M1 neuron, not pm1)
    assert terms["MC"] == "WBbt:0003638"  # class curation (MC neuron)
    # previously-broken non-neuron cell / endpoint-stub names now resolve
    assert terms["PM2D"] == "WBbt:0003632"  # individual pharyngeal muscle cell
    assert terms["INTMUL"] == "WBbt:0003833"
    assert terms["PM5"] == "WBbt:0003737"  # endpoint stub (pharyngeal muscle class)
    assert terms["BODYWALLMUSCLES"] == "WBbt:0006804"


def test_class_curation_takes_precedence_over_weak_match() -> None:
    """M2/M3 lexically collide with pm2/pm3 muscles (weak); curation forces the neuron term."""
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    index = WBBTIndex.from_obograph(WBBT_JSON)
    terms = anatomy_terms_map(connectome, index, load_class_curation(CLASS_CURATION))
    assert terms["M2"] == "WBbt:0003634"  # M2 neuron
    assert terms["M3"] == "WBbt:0003754"  # M3 neuron
