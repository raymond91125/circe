"""Unit tests for the Wang 2024 (eLife 95402) per-sex neurotransmitter ingest."""

from pathlib import Path

from celegans_connectome_kg.ingest.neurotransmitter import SOURCE, read_neurotransmitters

CSV = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "wang-neurotransmitter-atlas"
    / "sex_neurotransmitters.csv"
)


def test_read_neurotransmitters() -> None:
    recs = read_neurotransmitters(CSV)
    assert "95402" in SOURCE
    by = {(r.cell, r.sex): r for r in recs}
    # male-specific neurons (codes follow the neuron-graph scheme)
    assert by[("CEMDL", "male")].neurotransmitter == "a"
    assert by[("R7AL", "male")].neurotransmitter == "d"
    # dimorphic sex-shared neuron recorded for both sexes
    assert by[("AIML", "hermaphrodite")].neurotransmitter == "ls"
    assert by[("AIML", "male")].neurotransmitter == "as"
    # every row has a confidence and a note; sexes are constrained
    assert all(r.confidence in {"reported", "putative"} for r in recs)
    assert {r.sex for r in recs} == {"male", "hermaphrodite"}
    # no duplicate (cell, sex) after the parser's de-duplication
    keys = [(r.cell, r.sex) for r in recs]
    assert len(keys) == len(set(keys))
