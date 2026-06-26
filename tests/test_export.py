"""Phase 3 export tests: JSON round-trip and RDF serialization.

RDF is dumped only for a tiny hand-built Connectome (the full 25k-connection dump is slow);
the JSON round-trip is exercised on the real assembled graph.
"""

from pathlib import Path

import pytest
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import XSD

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.build.datamodel import datamodel
from celegans_connectome_kg.export.rdf import load_json, to_turtle, write_json

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
NEURON_GRAPH = DATA_DIR / "neuron-graph"
WBBT_JSON = DATA_DIR / "wbbt" / "wbbt.json"
BASE = "https://wormbase.org/resources/connectome/"


def _tiny_connectome() -> object:
    dm = datamodel()
    return dm.Connectome(
        cells=[
            dm.Cell(id="cckg:cell/ADAL", name="ADAL", cell_type="neuron", anatomy="WBbt:0004013")
        ],
        datasets=[dm.Dataset(id="cckg:dataset/white_1986_whole", name="White 1986 whole")],
        connections=[
            dm.Connection(
                id="cckg:conn/white_1986_whole.ADAL.ADAL.0",
                pre="cckg:cell/ADAL",
                post="cckg:cell/ADAL",
                connection_type="chemical",
                weight=16.0,
                dataset="cckg:dataset/white_1986_whole",
                evidence=[dm.Evidence(evidence_type="aggregated_weight", value="16")],
            )
        ],
    )


def test_turtle_parses_and_has_expected_triples() -> None:
    g = Graph()
    g.parse(data=to_turtle(_tiny_connectome()), format="turtle")

    cell = URIRef(f"{BASE}cell/ADAL")
    conn = URIRef(f"{BASE}conn/white_1986_whole.ADAL.ADAL.0")
    # Cell typed and grounded in WBBT.
    assert (
        cell,
        URIRef(f"{BASE}anatomy"),
        URIRef("http://purl.obolibrary.org/obo/WBbt_0004013"),
    ) in g
    # Connection weight serialized as xsd:float.
    assert (conn, URIRef(f"{BASE}weight"), Literal("16.0", datatype=XSD.float)) in g
    # pre/post reference the cell URI.
    assert (conn, URIRef(f"{BASE}pre"), cell) in g


@pytest.fixture(scope="module")
def built() -> object:
    connectome, _ = assemble(NEURON_GRAPH, WBBT_JSON)
    return connectome


def test_json_roundtrip_preserves_counts(built: object, tmp_path: Path) -> None:
    path = tmp_path / "connectome.json"
    write_json(built, path)
    reloaded = load_json(path)
    assert len(reloaded.cells) == len(built.cells) == 447
    assert len(reloaded.connections) == len(built.connections) == 25036
    assert len(reloaded.datasets) == len(built.datasets) == 15
    # A spot-checked cell survives the round-trip with its anatomy.
    adal = next(c for c in reloaded.cells if c.id == "cckg:cell/ADAL")
    assert str(adal.anatomy) == "WBbt:0004013"
