"""Ingest sex-qualified neurotransmitter assignments from the Wang et al. 2024 atlas (eLife 95402).

The build input is a vendored, pre-parsed CSV (``data/wang-neurotransmitter-atlas/
sex_neurotransmitters.csv``) produced by ``parse_male_atlas.py`` from Supplementary Files 3
(male-specific neurons) and 4 (sex-shared dimorphic neurons). Each row is one assertion:
``cell, sex, neurotransmitter, confidence, note``. The neurotransmitter code follows the
neuron-graph scheme (a=ACh, l=Glu, g=GABA, d=dopamine, s=serotonin, u=unknown). This module only
parses; the build maps cell names to KG cells and mints reified NeurotransmitterAssignment records.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

#: DOI of the source atlas, recorded as each assignment's provenance.
SOURCE = "https://doi.org/10.7554/eLife.95402"


@dataclass(frozen=True)
class NeurotransmitterRecord:
    cell: str
    sex: str  # "male" / "hermaphrodite"
    neurotransmitter: str  # neuron-graph-style code, e.g. "a", "ag", "u"
    confidence: str  # "reported" / "putative"
    note: str


def read_neurotransmitters(csv_path: Path) -> list[NeurotransmitterRecord]:
    """Read the vendored per-sex neurotransmitter assignments."""
    out: list[NeurotransmitterRecord] = []
    with Path(csv_path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            out.append(
                NeurotransmitterRecord(
                    cell=row["cell"].strip(),
                    sex=row["sex"].strip(),
                    neurotransmitter=row["neurotransmitter"].strip(),
                    confidence=(row.get("confidence") or "reported").strip(),
                    note=(row.get("note") or "").strip(),
                )
            )
    return out
