"""Tests for the node-name → WormAtlas URL map used by the viz cell-info link."""

from pathlib import Path

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.export.neuron_graph_json import wormatlas_links_map

REPO = Path(__file__).resolve().parents[1]
NEURON_GRAPH = REPO / "data" / "neuron-graph"
WBBT_JSON = REPO / "data" / "wbbt" / "wbbt.json"
CURATION = REPO / "data" / "curation" / "anatomy_curation.csv"
ENDPOINT_CELLS = REPO / "data" / "curation" / "connection_endpoint_cells.csv"

BWM_URL = "https://wormatlas.org/hermaphrodite/musclesomatic/MusSomaticframeset.html"


def _map() -> dict[str, str]:
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON, CURATION, ENDPOINT_CELLS)
    return wormatlas_links_map(connectome)


def test_neuron_links_use_class_page() -> None:
    m = _map()
    assert m["AVA"] == "http://www.wormatlas.org/neurons/Individual%20Neurons/AVAframeset.html"
    # an individual neuron cell points at its class page
    assert m["AVAL"].endswith("/AVAframeset.html")


def test_body_wall_muscles_use_somatic_muscle_page() -> None:
    m = _map()
    for key in ("BODYWALLMUSCLES", "BWM-DL01", "BWM01", "LEGACYBODYWALLMUSCLES"):
        assert m[key] == BWM_URL


def test_other_non_neurons_are_omitted() -> None:
    """Non-neuron, non-body-wall entities have no mapped WormAtlas page (link hidden)."""
    m = _map()
    # pharyngeal muscle / gland / hypodermis have no curated WormAtlas page yet
    for key in ("PM5", "G1", "HYP"):
        assert key not in m


def test_keys_are_upper_cased() -> None:
    m = _map()
    assert all(k == k.upper() for k in m)
