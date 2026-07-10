"""Phase 3 part 2: neuron-graph JSON projection — shape + round-trip fidelity.

The cells projection must round-trip losslessly against the raw neurons.json. The connections
projection must reproduce, from the assembled KG, the same grouping neuron-graph derives
directly from the source records — verifying build (id minting) and projection (id stripping
+ regrouping + gap-junction merge) preserve the data.
"""

import json
from collections import defaultdict
from pathlib import Path

import pytest

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.export.neuron_graph_json import (
    cells_projection,
    connections_projection,
)
from celegans_connectome_kg.ingest.neuron_graph import load_neuron_graph

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
NEURON_GRAPH = DATA_DIR / "neuron-graph"
WBBT_JSON = DATA_DIR / "wbbt" / "wbbt.json"

_API_TYPE = {"chemical": "chemical", "gap_junction": "electrical", "functional": "functional"}


@pytest.fixture(scope="module")
def connectome() -> object:
    c, _ = assemble(NEURON_GRAPH, WBBT_JSON)
    return c


def test_cells_projection_roundtrips_against_source(connectome: object) -> None:
    raw = json.loads((NEURON_GRAPH / "neurons.json").read_text())
    expected = {
        r["name"]: {
            "name": r["name"],
            "class": r.get("classes"),
            "type": r.get("typ"),
            "neurotransmitter": r.get("nt"),
            "embryonic": bool(r.get("emb")),
            "inhead": bool(r.get("inhead")),
            "intail": bool(r.get("intail")),
        }
        for r in raw
    }
    projected = cells_projection(connectome)
    assert len(projected) == len(expected) == 447
    assert {c["name"]: c for c in projected} == expected


def test_connection_objects_have_exact_shape(connectome: object) -> None:
    conns = connections_projection(connectome)
    for c in conns:
        assert set(c) == {"pre", "post", "type", "annotations", "synapses"}
        assert c["type"] in {"chemical", "electrical", "functional"}
        assert c["annotations"] == []
        assert c["synapses"] and all(isinstance(v, int) for v in c["synapses"].values())
        assert all(isinstance(k, str) for k in c["synapses"])


def test_gap_junctions_fully_merged(connectome: object) -> None:
    conns = connections_projection(connectome)
    electrical = [c for c in conns if c["type"] == "electrical"]
    keys = ["$".join(sorted([c["pre"], c["post"]])) for c in electrical]
    assert len(keys) == len(set(keys))  # no unmerged flipped pairs remain


def test_projection_matches_source_grouping(connectome: object) -> None:
    """Independently regroup the ingest records and compare to the KG-derived projection."""
    from celegans_connectome_kg.build.assemble import _NEURON_GRAPH_ALIAS as GA

    data = load_neuron_graph(NEURON_GRAPH)
    grouped: dict[tuple[str, str, str], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for r in data.connections:
        # Sum duplicate listings, as neuron-graph populate does; apply the same G1/G2->g1/g2
        # gland-misnomer rename the build applies to neuron-graph endpoints.
        pre, post = GA.get(r.pre, r.pre), GA.get(r.post, r.post)
        grouped[(pre, post, _API_TYPE[r.connection_type])][r.dataset_id] += int(r.weight)

    # Merge electrical the same way neuron-graph does.
    def regroup(items: dict) -> dict:
        electrical: dict[str, dict] = {}
        others: dict[tuple, dict] = {}
        for (pre, post, t), syn in items.items():
            if t == "electrical":
                key = "$".join(sorted([pre, post]))
                tgt = electrical.setdefault(key, {})
                for ds, n in syn.items():
                    tgt[ds] = tgt.get(ds, 0) + n
            else:
                others[(pre, post, t)] = syn
        return {"electrical": electrical, "others": others}

    expected = regroup(grouped)
    projected = connections_projection(connectome)

    proj_elec = {
        "$".join(sorted([c["pre"], c["post"]])): c["synapses"]
        for c in projected
        if c["type"] == "electrical"
    }
    proj_others = {
        (c["pre"], c["post"], c["type"]): c["synapses"]
        for c in projected
        if c["type"] != "electrical"
    }
    assert proj_elec == expected["electrical"]
    assert proj_others == expected["others"]


def test_total_synapse_weight_conserved(connectome: object) -> None:
    data = load_neuron_graph(NEURON_GRAPH)
    source_total = sum(int(c.weight) for c in data.connections)
    proj_total = sum(v for c in connections_projection(connectome) for v in c["synapses"].values())
    assert proj_total == source_total
