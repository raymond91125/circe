"""Manual anatomy curation: human-resolved cell → WBbt mappings.

The lexical matcher leaves a work-list tail (ambiguous + unmatched). A curator resolves those
against WBBT and records the confirmed mappings in ``data/curation/anatomy_curation.csv``.
This loader feeds them back into the match stage, where they take precedence over (and
correct) the lexical result — e.g. M1/M4/M5 are pharyngeal neurons that lexically collide
with the pm1/pm4/pm5 muscle synonyms.
"""

from __future__ import annotations

import csv
from pathlib import Path

DEFAULT_CURATION_PATH = Path("data/curation/anatomy_curation.csv")


def load_curation(path: Path) -> dict[str, str]:
    """Load cell_name → WBbt CURIE from a curation CSV (rows without a wbbt_id are skipped)."""
    curated: dict[str, str] = {}
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            wbbt_id = (row.get("wbbt_id") or "").strip()
            if wbbt_id:
                curated[row["cell_name"]] = wbbt_id
    return curated
