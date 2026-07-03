"""Tests for the neurotransmitter curation overlay (corrects neuron-graph `nt`)."""

from pathlib import Path

from celegans_connectome_kg.build.assemble import assemble
from celegans_connectome_kg.match.curation import load_nt_curation

REPO = Path(__file__).resolve().parents[1]
NEURON_GRAPH = REPO / "data" / "neuron-graph"
WBBT_JSON = REPO / "data" / "wbbt" / "wbbt.json"
NT_CURATION = REPO / "data" / "curation" / "neurotransmitter_curation.csv"


def test_load_nt_curation() -> None:
    cur = load_nt_curation(NT_CURATION)
    assert cur["HSNL"] == "as" and cur["HSNR"] == "as"  # glutamate+serotonin -> ACh+serotonin
    # atlas gap-fills (neuron-graph nt=u -> confident Wang 2024 assignment)
    assert cur["ALA"] == "g" and cur["AWAL"] == "a" and cur["I4"] == "l"
    # uptake-only cells deliberately NOT assigned
    assert "AVFL" not in cur and "ASIL" not in cur


def test_overlay_corrects_hsn_and_leaves_others() -> None:
    # without the overlay: neuron-graph's original nt (glutamate+serotonin)
    base, _ = assemble(NEURON_GRAPH, WBBT_JSON)
    base_by_name = {c.name: c for c in base.cells}
    assert base_by_name["HSNL"].neurotransmitter == "ls"

    # with the overlay: HSN corrected to ACh+serotonin
    curated, _ = assemble(NEURON_GRAPH, WBBT_JSON, nt_curation_path=NT_CURATION)
    by_name = {c.name: c for c in curated.cells}
    assert by_name["HSNL"].neurotransmitter == "as"
    assert by_name["HSNR"].neurotransmitter == "as"
    # gap-fills: unknown -> confident assignment
    assert base_by_name["ALA"].neurotransmitter == "u"
    assert by_name["ALA"].neurotransmitter == "g"
    assert by_name["AWAL"].neurotransmitter == "a"
    assert by_name["I4"].neurotransmitter == "l"
    # uptake-only cells left unknown
    assert by_name["AVFL"].neurotransmitter == "u"
    assert by_name["ASIL"].neurotransmitter == "u"
    # a non-curated neuron is untouched
    assert by_name["AVAL"].neurotransmitter == base_by_name["AVAL"].neurotransmitter
