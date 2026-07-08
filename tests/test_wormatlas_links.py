"""Tests for the node-name → WormAtlas URL map used by the viz cell-info link."""

from pathlib import Path

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.build.datamodel import datamodel
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


def test_male_class_curation_supplies_neuron_class_page() -> None:
    """Male-specific neuron cells (no cell_class) use the curated WormAtlas class slug."""
    dm = datamodel()
    connectome = dm.Connectome(
        cells=[
            dm.Cell(id="cckg:cell/R1AL", name="R1AL", cell_type="neuron", sexes=["male"]),
            dm.Cell(id="cckg:cell/CEMDL", name="CEMDL", cell_type="neuron", sexes=["male"]),
            dm.Cell(id="cckg:cell/HOA", name="HOA", cell_type="neuron", sexes=["male"]),
        ],
        datasets=[],
        connections=[],
    )
    curation = {"R1AL": "R1A", "CEMDL": "CEM"}  # HOA intentionally uncurated
    m = wormatlas_links_map(connectome, curation)
    assert m["R1AL"].endswith("/R1Aframeset.html")
    assert m["CEMDL"].endswith("/CEMframeset.html")
    # uncurated male neuron falls back to its own name (no crash, still a neuron page)
    assert m["HOA"].endswith("/HOAframeset.html")


def test_no_male_curation_falls_back_to_name() -> None:
    """Without curation, a classless male neuron uses its own name (prior behavior)."""
    dm = datamodel()
    connectome = dm.Connectome(
        cells=[dm.Cell(id="cckg:cell/R1AL", name="R1AL", cell_type="neuron", sexes=["male"])],
        datasets=[],
        connections=[],
    )
    m = wormatlas_links_map(connectome)
    assert m["R1AL"].endswith("/R1ALframeset.html")
