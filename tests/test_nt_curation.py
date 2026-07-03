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
    # applied gap-fills (atlas clean, non-hedged calls)
    assert cur["ALA"] == "g" and cur["RIPL"] == "a" and cur["RIPR"] == "a"
    # atlas-hedged ("potential") and uptake-only calls are NOT applied
    for held in ("AVJL", "AWAL", "I4", "AVFL", "ASIL"):
        assert held not in cur


def test_overlay_applies_clean_calls_and_holds_hedged() -> None:
    # without the overlay: neuron-graph's original nt
    base, _ = assemble(NEURON_GRAPH, WBBT_JSON)
    base_by_name = {c.name: c for c in base.cells}
    assert base_by_name["HSNL"].neurotransmitter == "ls"

    # with the overlay
    curated, _ = assemble(NEURON_GRAPH, WBBT_JSON, nt_curation_path=NT_CURATION)
    by_name = {c.name: c for c in curated.cells}
    # applied: HSN correction + ALA/RIP gap-fills
    assert by_name["HSNL"].neurotransmitter == "as"
    assert by_name["HSNR"].neurotransmitter == "as"
    assert base_by_name["ALA"].neurotransmitter == "u" and by_name["ALA"].neurotransmitter == "g"
    assert by_name["RIPL"].neurotransmitter == "a"
    # held (atlas-hedged "potential" / uptake-only): left unknown
    for held in ("AVJL", "AWAL", "I4", "AVFL", "ASIL"):
        assert by_name[held].neurotransmitter == "u"
    # a non-curated neuron is untouched
    assert by_name["AVAL"].neurotransmitter == base_by_name["AVAL"].neurotransmitter
