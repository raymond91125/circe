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
    # ray structural cells (Cook 'Rnsh' == Rnst; aliased sh->st, then grounded under the st name)
    assert curated["R1stL"] == "WBbt:0004044"
    assert curated["R8stR"] == "WBbt:0003972"
    assert "R1shL" not in curated and "R8shR" not in curated  # sh names no longer curated cells
    # M4b external-term-lookup groundings
    assert curated["ailL"] == "WBbt:0003790"  # anterior inner longitudinal muscle L (male)
    assert curated["vsrR"] == "WBbt:0004908"  # ventral spicule retractor R
    assert curated["dglL1"] == "WBbt:0008381"  # male diagonal muscle left (Cook 'dgl')
    assert curated["um1AL"] == "WBbt:0006915" and curated["vm1AL"] == "WBbt:0006917"
    assert curated["exc_cell"] == "WBbt:0005812" and curated["int"] == "WBbt:0005772"
    assert curated["bm"] == "WBbt:0005756"  # basement membrane (NSM secretion target)
    # e2 pharyngeal cells: Cook 2019 mislabels, relabeled to canonical (aliased e2D->e2DR etc.);
    # connectivity-confirmed via Cook 2020 pharynx connectome / WBBT.
    assert curated["e2DR"] == "WBbt:0004552"
    assert curated["e2DL"] == "WBbt:0004554"
    assert curated["e2V"] == "WBbt:0004550"
    assert not {"e2D", "e2VL", "e2VR"} & set(curated)  # old Cook labels no longer curated cells
    # curated Cook cells are genuinely outside the hermaphrodite registry
    assert all(cell not in NG for cell in curated)


def test_worklist_fully_curated() -> None:
    import csv
    from pathlib import Path

    wl = Path("data/curation/cook_curation_worklist.csv")
    rows = list(csv.DictReader(wl.open()))
    assert rows == [], f"M4 worklist should be empty; {len(rows)} residual cells remain"


def test_male_nonneuron_wormatlas_pages() -> None:
    """Male-specific non-neurons link to the general WormAtlas male-anatomy section pages."""
    from celegans_connectome_kg.match.curation import load_wormatlas_urls

    urls = load_wormatlas_urls("data/curation/cook_wormatlas_class.csv")
    # ray structural/support cells -> male rays page
    assert "male/rays/Rayframeset" in urls["R3stL"] and "male/rays/Rayframeset" in urls["R9stL"]
    # M-derived sex muscles -> male-specific muscle page (spicule, gubernacular, diagonal, ...)
    for muscle in ("dglL1", "vsrR", "dspL", "gecL", "ailL"):
        assert "male/musclemale/Musmaleframeset" in urls[muscle]
    # gonad / proctodeum -> male reproductive-system page
    assert "male/reproductive" in urls["gonad"] and "male/reproductive" in urls["proctodeum"]
    # every curated WormAtlas URL is https
    assert all(u.startswith("https://") for u in urls.values())


def test_shared_nonneuron_wormatlas_pages() -> None:
    """Shared (male+herm) non-neurons link to the hermaphrodite/general handbook pages."""
    from celegans_connectome_kg.match.curation import load_wormatlas_urls

    urls = load_wormatlas_urls("data/curation/cook_wormatlas_class.csv")
    # pharyngeal cells (muscle/marginal/epithelial/gland) -> pharynx handbook
    for cell in ("pm3D", "mc2DL", "e2DL", "g1AL"):
        assert "hermaphrodite/pharynx/" in urls[cell]
    assert "muscleGLR" in urls["GLRDL"]  # GLR cells
    assert "neuronalsupport" in urls["CEPshDL"]  # CEP sheath (glia)
    assert "Individual%20Neurons/CANframeset" in urls["CANL"]  # CAN is a neuron
    assert "musclenonstriated" in urls["mu_anal"]  # enteric muscle
    assert "excretory" in urls["exc_cell"] and "intestine" in urls["int"]
    assert "muscleheadcell" in urls["HMC"] and "hypodermis" in urls["hyp"]
    assert "bm" not in urls  # basal lamina = ECM, no handbook cell page


def test_naming_variant_aliases() -> None:
    aliases = load_cook_aliases(DEFAULT_COOK_ALIASES_PATH)
    # Cook naming variants of existing neuron-graph cells
    assert aliases["hmc"] == "HMC" and aliases["exc_gl"] == "excgl"
    assert aliases["mu_intL"] == "intmuL" and aliases["mu_intR"] == "intmuR"
    assert all(t in NG for t in ("HMC", "excgl", "intmuL", "intmuR"))
    # Cook 'Rnsh' (SI-only) == ray structural cell Rnst; aliased to merge connectivity
    assert aliases["R1shL"] == "R1stL" and aliases["R8shR"] == "R8stR"
    assert "R1shL" not in NG
    # Cook e2 relabeling to canonical pharyngeal-muscle names
    assert aliases["e2D"] == "e2DR" and aliases["e2VL"] == "e2DL" and aliases["e2VR"] == "e2V"
