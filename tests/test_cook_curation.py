"""Tests for the Cook 2019 name reconciliation + anatomy grounding [M4]."""

from celegans_connectome_kg.ingest.cook_2019 import read_cook
from celegans_connectome_kg.ingest.neuron_graph import read_cells
from celegans_connectome_kg.match.curation import (
    DEFAULT_COOK_ALIASES_PATH,
    DEFAULT_COOK_ANATOMY_CURATION_PATH,
    load_cook_aliases,
    load_curation,
)

NG = {c.name for c in read_cells("data/neuron-graph/neurons.json")}


def test_aliases_reconcile_naming() -> None:
    aliases = load_cook_aliases(DEFAULT_COOK_ALIASES_PATH)
    # body-wall muscle renaming -> the neuron-graph BWM cells (shared, already grounded)
    assert aliases["dBWML1"] == "BWM-DL01"
    assert aliases["vBWMR9"] == "BWM-VR09"
    # zero-padded motor-neuron series -> canonical (shared)
    assert aliases["DA01"] == "DA1" and aliases["AS09"] == "AS9"
    # male-specific series normalized to canonical (still absent from the herm registry)
    assert aliases["CA01"] == "CA1" and "CA1" not in NG
    # every alias key is a real Cook cell, and the map is non-trivial
    cook_cells = set().union(*read_cook().cells_by_sex.values())
    assert set(aliases) <= cook_cells
    assert len(aliases) >= 150


def test_aliases_shared_targets_exist_in_registry() -> None:
    aliases = load_cook_aliases(DEFAULT_COOK_ALIASES_PATH)
    # BWM aliases must all resolve to an existing neuron-graph cell
    bwm = {c: t for c, t in aliases.items() if t.startswith("BWM-")}
    assert len(bwm) == 95
    assert all(t in NG for t in bwm.values())


def test_cook_anatomy_curation_grounds_male_cells() -> None:
    curated = load_curation(DEFAULT_COOK_ANATOMY_CURATION_PATH)
    # every value is a WBbt CURIE
    assert curated and all(v.startswith("WBbt:") for v in curated.values())
    # known male-specific / new cells are grounded
    for cell in ("CEMDL", "CA9", "CP9", "HOA", "PCAL", "proctodeum"):
        assert cell in curated
    # ray sheath cells map to their WBbt ray structural term (Cook sh == WBbt st)
    assert curated["R1shL"] == "WBbt:0004044"  # R1stL
    assert curated["R8shR"] == "WBbt:0003972"  # R8stR
    # M4b external-term-lookup groundings
    assert curated["ailL"] == "WBbt:0003790"  # anterior inner longitudinal muscle L (male)
    assert curated["vsrR"] == "WBbt:0004908"  # ventral spicule retractor R
    assert curated["dglL1"] == "WBbt:0008381"  # male diagonal muscle left (Cook 'dgl')
    assert curated["um1AL"] == "WBbt:0006915" and curated["vm1AL"] == "WBbt:0006917"
    assert curated["exc_cell"] == "WBbt:0005812" and curated["int"] == "WBbt:0005772"
    assert curated["bm"] == "WBbt:0005756"  # basement membrane (NSM secretion target)
    # e2 pharyngeal cells: Cook 2019 mislabels; connectivity-confirmed canonical identity
    # (Cook 2020 pharynx connectome / WBBT): e2D=e2DR, e2VL=e2DL, e2VR=e2V
    assert curated["e2D"] == "WBbt:0004552"  # e2DR
    assert curated["e2VL"] == "WBbt:0004554"  # e2DL
    assert curated["e2VR"] == "WBbt:0004550"  # e2V
    # curated Cook cells are genuinely outside the hermaphrodite registry
    assert all(cell not in NG for cell in curated)


def test_worklist_fully_curated() -> None:
    import csv
    from pathlib import Path

    wl = Path("data/curation/cook_curation_worklist.csv")
    rows = list(csv.DictReader(wl.open()))
    assert rows == [], f"M4 worklist should be empty; {len(rows)} residual cells remain"


def test_naming_variant_aliases() -> None:
    aliases = load_cook_aliases(DEFAULT_COOK_ALIASES_PATH)
    # Cook naming variants of existing neuron-graph cells
    assert aliases["hmc"] == "HMC" and aliases["exc_gl"] == "excgl"
    assert aliases["mu_intL"] == "intmuL" and aliases["mu_intR"] == "intmuR"
    assert all(t in NG for t in ("HMC", "excgl", "intmuL", "intmuR"))
