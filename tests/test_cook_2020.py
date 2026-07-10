"""Tests for the Cook 2020 pharyngeal connectome ingest (SI 3 aggregation)."""

from celegans_connectome_kg.ingest.cook_2020 import DATASET_ID, read_cook_2020


def test_dataset_metadata() -> None:
    data = read_cook_2020()
    assert len(data.datasets) == 1
    d = data.datasets[0]
    assert d.id == DATASET_ID == "cook_2020_pharynx"
    assert d.sex == "hermaphrodite"  # reconstructed from N2 hermaphrodite specimens
    assert "2020" in d.name and "pharynx" in d.name.lower()


def test_edges_counts_and_types() -> None:
    data = read_cook_2020()
    conns = data.connections
    chem = [c for c in conns if c.connection_type == "chemical"]
    gap = [c for c in conns if c.connection_type == "gap_junction"]
    # SI 3 aggregation: 332 edges (276 chemical, 56 electrical->gap_junction)
    assert len(conns) == 332
    assert len(chem) == 276
    assert len(gap) == 56
    assert all(c.dataset_id == DATASET_ID for c in conns)


def test_weight_is_summed_sections() -> None:
    """Weights are total EM sections (Cook-2019 metric): integer, e.g. NSMR->bm = 72."""
    conns = read_cook_2020().connections
    by = {(c.pre, c.post, c.connection_type): c.weight for c in conns}
    assert by[("NSMR", "bm", "chemical")] == 72
    assert by[("M5", "pm7D", "chemical")] == 8
    assert all(float(c.weight).is_integer() for c in conns)


def test_gap_junctions_canonical_and_no_reciprocals() -> None:
    conns = read_cook_2020().connections
    gap = [(c.pre, c.post) for c in conns if c.connection_type == "gap_junction"]
    assert all(a <= b for a, b in gap)  # canonical (sorted) orientation
    assert not any((b, a) in set(gap) for a, b in gap if a != b)  # undirected, one per pair


def test_noncell_endpoints_excluded() -> None:
    """Raw object ids / unk are dropped by the aggregator (g1vl/g1vr are kept, compiled to g1
    by cook_name_aliases at assemble time — so they still appear in the raw edge list here)."""
    conns = read_cook_2020().connections
    names = {c.pre for c in conns} | {c.post for c in conns}
    assert not any(n.lower().startswith(("obj", "unk")) for n in names)
