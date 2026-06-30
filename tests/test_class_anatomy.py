"""Tests for the cell-class → WBbt map that drives the viz cell-info WormBase link."""

from pathlib import Path

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.export.neuron_graph_json import class_anatomy_map
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


def test_class_anatomy_map_covers_every_class() -> None:
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    index = WBBTIndex.from_obograph(WBBT_JSON)
    terms = class_anatomy_map(connectome, index, load_class_curation(CLASS_CURATION))

    # every neuron-graph cell class is mapped
    classes = {c.cell_class for c in connectome.cells if c.cell_class}
    assert classes <= set(terms)
    assert len(terms) == 147
    assert all(v.startswith("WBbt:") for v in terms.values())


def test_class_map_resolution_sources() -> None:
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    index = WBBTIndex.from_obograph(WBBT_JSON)
    terms = class_anatomy_map(connectome, index, load_class_curation(CLASS_CURATION))

    assert terms["AVA"] == "WBbt:0005842"  # strong class-name match (the reported example)
    assert terms["DVB"].startswith("WBbt:")  # singleton-class reuse of the cell's anatomy
    assert terms["M1"] == "WBbt:0004488"  # singleton reuse of corrected cell curation (M1 neuron)
    assert terms["MC"] == "WBbt:0003638"  # class curation (MC neuron)


def test_class_curation_takes_precedence_over_weak_match() -> None:
    """M2/M3 lexically collide with pm2/pm3 muscles (weak); curation forces the neuron term."""
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    index = WBBTIndex.from_obograph(WBBT_JSON)
    terms = class_anatomy_map(connectome, index, load_class_curation(CLASS_CURATION))
    assert terms["M2"] == "WBbt:0003634"  # M2 neuron
    assert terms["M3"] == "WBbt:0003754"  # M3 neuron
