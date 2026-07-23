"""Unit tests for the Yim et al. 2024 dauer nerve-ring connectome ingest."""

from pathlib import Path

from celegans_connectome_kg.ingest.yim_2024_dauer import DATASET_ID, read_dauer

CSV = Path(__file__).resolve().parents[1] / "data" / "yim-2024-dauer" / "dauer_connections.csv"


def test_read_dauer() -> None:
    data = read_dauer(CSV)
    assert data.dataset_id == DATASET_ID == "yim_2024_dauer"
    assert data.sex == "hermaphrodite"
    conns = data.connections
    # 2,200 directed edges aggregated from 6,371 synapses; chemical only; weight = synapse count
    assert len(conns) == 2200
    assert all(c.connection_type == "chemical" and c.weight > 0 for c in conns)
    assert int(sum(c.weight for c in conns)) == 6371
    # spot-check a specific edge and the excretory-duct target (a non-neuronal partner)
    by = {(c.pre, c.post): c.weight for c in conns}
    assert by[("ADAL", "AIAL")] == 1.0
    assert ("AVM", "exc_duct") in by
