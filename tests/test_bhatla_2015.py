"""Unit tests for the Bhatla et al. 2015 I2 EM synapse ingest."""

from pathlib import Path

from celegans_connectome_kg.ingest.bhatla_2015 import DATASET_ID, read_bhatla_i2

CSV = Path(__file__).resolve().parents[1] / "data" / "bhatla-2015-i2" / "i2_synapses.csv"


def test_read_bhatla_i2() -> None:
    data = read_bhatla_i2(CSV)
    assert data.dataset_id == DATASET_ID == "bhatla_2015_i2"
    assert data.sex == "hermaphrodite"
    conns = data.connections
    assert len(conns) == 26
    # all presynaptic I2, all chemical, weight = EM sections (positive)
    assert {c.pre for c in conns} == {"I2L", "I2R"}
    assert all(c.connection_type == "chemical" and c.weight > 0 for c in conns)
    # spot-check the strongest edge and the basal-lamina ("bm") synapse
    by = {(c.pre, c.post): c.weight for c in conns}
    assert by[("I2L", "pm3VL")] == 133.0
    assert by[("I2R", "bm")] == 2.0  # bm = basal lamina (WBbt:0005756), not a muscle
    # I2->pharyngeal-muscle synapses are present (heavily weighted here vs. minimal in Cook 2020)
    assert any(c.post.startswith("pm") for c in conns)
