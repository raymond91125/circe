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
    # a non-curated neuron is untouched
    assert by_name["AVAL"].neurotransmitter == base_by_name["AVAL"].neurotransmitter
