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
    # betaine (b) added for RIM: glutamate+tyramine+betaine (Hardege et al. 2022)
    assert cur["RIML"] == "blt" and cur["RIMR"] == "blt"
    # betaine (uptake) neurons from Wang 2024 Table S2 (Neurotransmitter(s) block U/V/W):
    # cat-1/VMAT + snf-3 co-expression -> betaine uptake+release (cf. RIH 5-HT uptake)
    assert cur["ASIL"] == "b" and cur["ASIR"] == "b"  # primary betaine
    assert cur["AUAL"] == "bl" and cur["AUAR"] == "bl"  # glutamate + betaine
    assert cur["CANL"] == "b" and cur["CANR"] == "b"  # orphan -> betaine
    assert cur["RIR"] == "ab"  # ACh + betaine
    # firmer (un-hedged) U/V/W co-transmitters from Wang Table S2
    assert cur["AFDL"] == "al" and cur["DVA"] == "al" and cur["M5"] == "al"
    assert cur["PDEL"] == "dl" and cur["RICL"] == "lo" and cur["MI"] == "ls"
    assert cur["AIBL"] == "bl" and cur["PHCL"] == "bl" and cur["RIS"] == "bg"
    assert cur["PVNL"] == "abl"
    # atlas-hedged ("*") calls, uptake-without-VGAT (AVF), and male-only (PVW) are NOT applied
    for held in ("AVJL", "AWAL", "I4", "AVFL", "PVW", "ASGL", "URXL", "DVB", "SMDDL"):
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
    assert by_name["RIML"].neurotransmitter == "blt"  # + betaine (Hardege 2022)
    # betaine (uptake) neurons (Wang 2024 Table S2): cat-1/VMAT + snf-3
    assert base_by_name["ASIL"].neurotransmitter == "u" and by_name["ASIL"].neurotransmitter == "b"
    assert base_by_name["AUAL"].neurotransmitter == "l" and by_name["AUAL"].neurotransmitter == "bl"
    assert base_by_name["CANL"].neurotransmitter == "n" and by_name["CANL"].neurotransmitter == "b"
    assert base_by_name["RIR"].neurotransmitter == "a" and by_name["RIR"].neurotransmitter == "ab"
    # firmer U/V/W co-transmitters (Wang Table S2), non-betaine and betaine-uptake
    assert base_by_name["AFDL"].neurotransmitter == "l" and by_name["AFDL"].neurotransmitter == "al"
    assert base_by_name["PDEL"].neurotransmitter == "d" and by_name["PDEL"].neurotransmitter == "dl"
    assert base_by_name["AIBL"].neurotransmitter == "l" and by_name["AIBL"].neurotransmitter == "bl"
    assert base_by_name["PVNL"].neurotransmitter == "a" and by_name["PVNL"].neurotransmitter == "abl"
    # held: hedged ("*") and uptake-without-VGAT (AVF) -> unchanged
    for held in ("AVJL", "AWAL", "I4", "AVFL"):
        assert by_name[held].neurotransmitter == "u"
    # a non-curated neuron is untouched
    assert by_name["AVAL"].neurotransmitter == base_by_name["AVAL"].neurotransmitter
