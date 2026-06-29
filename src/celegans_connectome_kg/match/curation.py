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
