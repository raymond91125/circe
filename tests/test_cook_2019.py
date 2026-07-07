"""Tests for the Cook et al. 2019 SI 5 adjacency-matrix ingest (sex-aware) [M3]."""

from collections import Counter

from celegans_connectome_kg.ingest.cook_2019 import DEFAULT_COOK_XLSX, read_cook


def _load():
    assert DEFAULT_COOK_XLSX.exists(), f"pinned SI 5 missing at {DEFAULT_COOK_XLSX}"
    return read_cook()


def test_datasets_are_sex_tagged() -> None:
    data = _load()
    by_id = {d.id: d for d in data.datasets}
    assert by_id["cook_2019_hermaphrodite"].sex == "hermaphrodite"
    assert by_id["cook_2019_male"].sex == "male"


def test_connection_counts_and_types() -> None:
    data = _load()
    counts = Counter((c.dataset_id, c.connection_type) for c in data.connections)
    # exact counts against the pinned (checksummed) SI 5 workbook
    assert counts[("cook_2019_hermaphrodite", "chemical")] == 4879
    assert counts[("cook_2019_hermaphrodite", "gap_junction")] == 1450
    assert counts[("cook_2019_male", "chemical")] == 5306
    assert counts[("cook_2019_male", "gap_junction")] == 1758
    # only the two EM connection types, all weights positive
    assert {c.connection_type for c in data.connections} == {"chemical", "gap_junction"}
    assert all(c.weight > 0 for c in data.connections)
    # Cook adjacency has no per-contact array (that detail isn't in SI 5)
    assert all(c.syn == () for c in data.connections)


def test_gap_junctions_deduped_to_one_per_unordered_pair() -> None:
    data = _load()
    for dataset_id in ("cook_2019_hermaphrodite", "cook_2019_male"):
        gaps = [
            c
            for c in data.connections
            if c.dataset_id == dataset_id and c.connection_type == "gap_junction"
        ]
        keys = [frozenset((c.pre, c.post)) for c in gaps]
        assert len(keys) == len(set(keys)), "symmetric gap matrix should yield one edge per pair"
        # canonical orientation: pre <= post
        assert all(c.pre <= c.post for c in gaps)


def test_sex_specific_cells() -> None:
    data = _load()
    male, herm = data.cells_by_sex["male"], data.cells_by_sex["hermaphrodite"]
    # male-specific
    for cell in ("CEMDL", "CA01", "CP01"):
        assert cell in male and cell not in herm
    # hermaphrodite-specific
    for cell in ("HSNL", "VC01"):
        assert cell in herm and cell not in male
    # shared core neuron present in both
    assert "AVAL" in male and "AVAL" in herm


def test_known_male_chemical_edge() -> None:
    data = _load()
    male_chem = {
        (c.pre, c.post): c.weight
        for c in data.connections
        if c.dataset_id == "cook_2019_male" and c.connection_type == "chemical"
    }
    assert male_chem[("AVAL", "AVAR")] == 12.0
