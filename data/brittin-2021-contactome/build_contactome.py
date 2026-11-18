"""Emit the L4 + adult nerve-ring contactome edge lists from Brittin et al. 2021.

The "contactome" (term coined here) is the set of physical membrane contacts (appositions) between
neurons -- the substrate that bounds the synaptic connectome, since a synapse can only form where
membranes touch (so synapses are a subset of contacts). Brittin et al. volumetrically reconstructed
membrane contacts in the L4 (JSH) and adult (N2U) nerve rings -- the SAME two EM series White et al.
1986 scored synapses in -- so the synapse-subset-of-contacts property is a clean within-animal check
(holds for ~99% of L4 and ~98% of adult White synapse pairs). This reads the raw per-section contact
scoring and emits the two committed edge lists the build ingests as ``adjacency`` connections.

Source (not vendored -- large third-party file, ~8.9 MB):
  Brittin CA, Cook SJ, Hall DH, Emmons SW, Cohen N. "A multiscale brain map derived from whole-brain
  volumetric reconstructions." Nature 591:105 (2021). PMID 33627874; doi:10.1038/s41586-021-03284-x.
  Supplementary Information 1 (``NIHMS2028296-supplement-Supp_Info_1.csv``): a header-less plain-text
  CSV, one row per membrane-contact instance scored in one EM section, with columns
  ``[dataset, cell1, cell2, cell1_index, cell2_index, EM_section_number, membrane_contact_in_pixels]``
  where ``dataset`` is ``l4`` (JSH) or ``adult`` (N2U). A single cell pair recurs across EM sections
  and across neurite-segment index pairs, so the total contact is the sum of the pixel column over all
  those rows; pixels convert to membrane-contact area in um^2 by scaling by 4.5e-4 (both series).

    python -I build_contactome.py <...Supp_Info_1.csv>    # or set BRITTIN_CONTACTOME_CSV
"""

from __future__ import annotations

import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
PX_TO_UM2 = 4.5e-4  # Brittin et al. 2021 Supp Info 1: pixels -> um^2 (same scale for L4 and adult)

#: Brittin zero-pads single-digit ventral-cord motor-neuron indices (VB01, VC01, ...); CIRCE and
#: White/Cook use the unpadded form (VB1, VC1). Strip the pad so every label grounds to a CIRCE cell.
_UNPAD = re.compile(r"(?<=[A-Z])0(\d)$")

#: dataset tag in the CSV -> (committed filename, CIRCE dataset id) for the message only
_DATASETS = {"l4": "l4_contactome.csv", "adult": "adult_contactome.csv"}


def norm(name: str) -> str:
    return _UNPAD.sub(r"\1", name.strip())


def src_path() -> Path:
    if len(sys.argv) > 1:
        return Path(sys.argv[1]).expanduser()
    env = os.environ.get("BRITTIN_CONTACTOME_CSV")
    if env:
        return Path(env).expanduser()
    sys.exit(
        "usage: python -I build_contactome.py <...Supp_Info_1.csv>  (or set BRITTIN_CONTACTOME_CSV)"
    )


def main() -> None:
    # (dataset, frozenset{cell_a, cell_b}) -> summed contact pixels
    px: dict[tuple[str, frozenset], int] = defaultdict(int)
    with open(src_path(), newline="") as f:
        for r in csv.reader(f):
            if len(r) != 7:
                continue
            ds, c1, c2 = r[0].strip(), norm(r[1]), norm(r[2])
            if ds not in _DATASETS or c1 == c2:
                continue
            px[(ds, frozenset((c1, c2)))] += int(r[6])

    for ds, fname in _DATASETS.items():
        edges = sorted(
            (tuple(sorted(pair)) + (round(p * PX_TO_UM2, 6),))
            for (d, pair), p in px.items()
            if d == ds
        )
        with open(HERE / fname, "w", newline="") as f:
            wri = csv.writer(f)
            wri.writerow(["cell_a", "cell_b", "contact_area"])
            for a, b, area in edges:
                wri.writerow([a, b, area])
        cells = {c for a, b, _ in edges for c in (a, b)}
        print(f"wrote {fname}: {len(edges)} undirected contacts over {len(cells)} cells")


if __name__ == "__main__":
    main()
