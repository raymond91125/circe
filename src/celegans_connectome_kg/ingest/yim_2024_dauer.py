"""Ingest the Yim, Choe, Bae et al. 2024 dauer nerve-ring connectome.

Yim H, Choe DT, Bae JA, Choi M-K, Kang H-M, Nguyen KCQ, et al. (2024). "Comparative connectomics
of dauer reveals developmental plasticity." Nat Commun 15:1546 (PMID 38413604). A complete
electron-microscopy reconstruction of the C. elegans dauer nerve ring, with deep-learning synapse
detection. Dauer is the stress-induced alternative third larval stage; the animal is a
hermaphrodite. Scope is the nerve ring only.

The vendored ``dauer_connections.csv`` (pre, post, synapses, size_sum) is aggregated from the
per-synapse Supplementary Data table. Weight = the number of synapses per ordered pair, matching
the synapse-count convention of the White/Witvliet datasets already in the KG; the paper's own
active-zone size metric is retained in ``size_sum`` for reference but not loaded. All connections
are chemical: the study did not reconstruct gap junctions.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from celegans_connectome_kg.ingest.neuron_graph import ConnectionRecord

DATASET_ID = "yim_2024_dauer"


@dataclass(frozen=True)
class DauerData:
    connections: list[ConnectionRecord]
    dataset_id: str
    dataset_name: str
    dataset_description: str
    sex: str


def read_dauer(csv_path: Path) -> DauerData:
    """Read the vendored dauer connectome into ConnectionRecords (weight = synapse count)."""
    connections = []
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            connections.append(
                ConnectionRecord(
                    dataset_id=DATASET_ID,
                    pre=row["pre"],
                    post=row["post"],
                    connection_type="chemical",
                    weight=float(row["synapses"]),
                    syn=(),
                    ids=None,
                    pre_tid=None,
                    post_tid=None,
                )
            )
    return DauerData(
        connections=connections,
        dataset_id=DATASET_ID,
        dataset_name="Yim et al. 2024 (dauer nerve ring)",
        dataset_description=(
            "Electron-microscopy reconstruction of the C. elegans dauer nerve ring, from Yim, "
            "Choe, Bae et al. 2024 (Nat Commun 15:1546). Chemical synapses only (gap junctions "
            "were not reconstructed). Weight = number of synapses per connection. Dauer is the "
            "stress-induced alternative third larval stage (hermaphrodite)."
        ),
        sex="hermaphrodite",
    )
