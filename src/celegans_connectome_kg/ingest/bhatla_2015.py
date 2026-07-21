"""Ingest the Bhatla et al. 2015 I2 EM synapse dataset.

Bhatla N, Droste R, Sando SR, Huang A, Horvitz HR (2015). "Distinct neural circuits control
rhythm inhibition and spitting by the myogenic pharynx of C. elegans." Curr Biol 25(16):2075-2089
(PMID 26212880). Supplemental Data 1, Figure S6 panel E: per-partner chemical synapses of the
pharyngeal I2 neurons (I2L, I2R), reconstructed by electron microscopy.

The source table is image-only in the PDF, so it is transcribed to a vendored CSV
(``i2_synapses.csv``: pre, post, sections, synapses) with partner names normalized to CIRCE cell
names. Weight = the number of EM sections over which the synapses onto each recipient are
distributed. Rows with zero synapses in this work (reported by Albertson & Thomson 1976 but not
confirmed here) are not vendored.

This dataset adds no new I2 cell partners: every partner (pm1-pm5, e3, I1, I4, I6, M1, M3, NSM,
and the basal lamina) already appears in the White and/or Cook connectomes. Its distinctive
contribution is the much larger synaptic weight on the I2->pharyngeal-muscle projection, which
Cook 2020 records only minimally. The postsynaptic ``bm`` partner is the basal lamina
(WBbt:0005756), not a muscle.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from celegans_connectome_kg.ingest.neuron_graph import ConnectionRecord

DATASET_ID = "bhatla_2015_i2"


@dataclass(frozen=True)
class BhatlaI2Data:
    connections: list[ConnectionRecord]
    dataset_id: str
    dataset_name: str
    dataset_description: str
    sex: str


def read_bhatla_i2(csv_path: Path) -> BhatlaI2Data:
    """Read the vendored I2 synapse table into ConnectionRecords (weight = EM sections)."""
    connections = []
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            connections.append(
                ConnectionRecord(
                    dataset_id=DATASET_ID,
                    pre=row["pre"],
                    post=row["post"],
                    connection_type="chemical",
                    weight=float(row["sections"]),
                    syn=(),
                    ids=None,
                    pre_tid=None,
                    post_tid=None,
                )
            )
    return BhatlaI2Data(
        connections=connections,
        dataset_id=DATASET_ID,
        dataset_name="Bhatla et al. 2015 (I2 EM synapses)",
        dataset_description=(
            "Electron-microscopy reconstruction of all chemical synapses of the pharyngeal I2 "
            "neurons (I2L, I2R), from Bhatla et al. 2015 (Curr Biol 25:2075-2089, Fig S6E). "
            "Weight = EM serial sections spanned by the synapses onto each partner. Every partner "
            "already appears in the White/Cook connectomes; the distinctive contribution is the "
            "much greater weight on the I2->pharyngeal-muscle projection, recorded only minimally "
            "by Cook 2020 (the 'bm' partner is the basal lamina, not a muscle)."
        ),
        sex="hermaphrodite",
    )
