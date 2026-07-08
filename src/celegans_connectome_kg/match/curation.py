"""Manual anatomy curation: human-resolved cell → WBbt mappings.

The lexical matcher leaves a work-list tail (ambiguous + unmatched). A curator resolves those
against WBBT and records the confirmed mappings in ``data/curation/anatomy_curation.csv``.
This loader feeds them back into the match stage, where they take precedence over (and
correct) the lexical result — e.g. M1/M4/M5 are pharyngeal neurons that lexically collide
with the pm1/pm4/pm5 muscle synonyms.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

DEFAULT_CURATION_PATH = Path("data/curation/anatomy_curation.csv")
DEFAULT_ENDPOINT_CELLS_PATH = Path("data/curation/connection_endpoint_cells.csv")
DEFAULT_CLASS_CURATION_PATH = Path("data/curation/class_anatomy_curation.csv")
DEFAULT_NT_CURATION_PATH = Path("data/curation/neurotransmitter_curation.csv")
# Cook et al. 2019 (sex extension): naming reconciliation + male-specific anatomy grounding.
DEFAULT_COOK_ALIASES_PATH = Path("data/curation/cook_name_aliases.csv")
DEFAULT_COOK_ANATOMY_CURATION_PATH = Path("data/curation/cook_anatomy_curation.csv")


@dataclass(frozen=True)
class EndpointCell:
    """A curated stub cell for a class-level / aggregate connection endpoint.

    These names appear as connection pre/post but are not in neuron-graph's cell list; we
    mint a Cell so the references resolve and carry WBBT grounding. They are KG-only (not in
    neuron-graph's /api/cells), so they are excluded from the viz cells projection.
    """

    name: str
    wbbt_id: str
    cell_type: str


def load_curation(path: Path) -> dict[str, str]:
    """Load cell_name → WBbt CURIE from a curation CSV (rows without a wbbt_id are skipped)."""
    curated: dict[str, str] = {}
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            wbbt_id = (row.get("wbbt_id") or "").strip()
            if wbbt_id:
                curated[row["cell_name"]] = wbbt_id
    return curated


def load_cook_aliases(path: Path) -> dict[str, str]:
    """Load Cook 2019 cell name → canonical KG name (naming reconciliation, sex extension).

    Reconciles Cook's naming to the KG's canonical names: zero-padded series (``DA01`` → ``DA1``)
    and body-wall muscles (``dBWML1`` → ``BWM-DL01``). A target already in the neuron-graph
    registry is a shared cell; a target absent from it is a male-specific cell renamed to
    canonical form (grounded via ``cook_anatomy_curation.csv``).
    """
    aliases: dict[str, str] = {}
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            canonical = (row.get("canonical_name") or "").strip()
            if canonical:
                aliases[row["cook_name"]] = canonical
    return aliases


def load_class_curation(path: Path) -> dict[str, str]:
    """Load cell-class → WBbt CURIE for classes the lexical class match can't resolve."""
    curated: dict[str, str] = {}
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            wbbt_id = (row.get("wbbt_id") or "").strip()
            if wbbt_id:
                curated[row["class_name"]] = wbbt_id
    return curated


def load_wormatlas_class(path: Path) -> dict[str, str]:
    """Load cell_name → WormAtlas neuron-class page slug for male-specific cells.

    Male-specific cells (Cook) are minted without a ``cell_class``, so the WormAtlas
    Individual-Neurons URL (keyed by neuron *class*, e.g. ``CEMDL`` → ``CEM``, ``R1AL`` → ``R1A``)
    can't be derived from the KG. This hand-curated map supplies the class slug; rows without a
    ``wormatlas_class`` are skipped (no link emitted).
    """
    curated: dict[str, str] = {}
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            slug = (row.get("wormatlas_class") or "").strip()
            if slug:
                curated[row["cell_name"]] = slug
    return curated


def load_nt_curation(path: Path) -> dict[str, str]:
    """Load cell_name → corrected neurotransmitter code, overriding neuron-graph's `nt`.

    Corrections are evidence-based (WormAtlas / Wang et al. atlas); see the reconciliation in
    ``analysis/neurotransmitter_reconciliation.md``.
    """
    curated: dict[str, str] = {}
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            nt = (row.get("neurotransmitter") or "").strip()
            if nt:
                curated[row["cell_name"]] = nt
    return curated


def load_endpoint_cells(path: Path) -> list[EndpointCell]:
    """Load curated stub cells for class-level/aggregate connection endpoints."""
    out: list[EndpointCell] = []
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            wbbt_id = (row.get("wbbt_id") or "").strip()
            if wbbt_id:
                out.append(
                    EndpointCell(
                        name=row["cell_name"],
                        wbbt_id=wbbt_id,
                        cell_type=(row.get("cell_type") or "other").strip(),
                    )
                )
    return out
