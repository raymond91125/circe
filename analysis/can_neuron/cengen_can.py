#!/usr/bin/env python3
"""CAN expression of neurotransmission-machinery genes, from CeNGEN scRNA-seq.

CIRCE's own gene layer carries only the 155 connectome-signaling genes (neuropeptides,
their receptors, innexins, monoamine/ionotropic/metabotropic receptors) — NOT the
synaptic-vesicle / active-zone / dense-core release machinery. So the "does CAN express
synaptotagmin?" question is answered from CeNGEN (Taylor et al. 2021), not from CIRCE.

This pulls CeNGEN's binarized threshold matrices + TPM table from the cengenDataSC R
package and extracts the CAN column (CeNGEN merges CANL/CANR into class "CAN";
hermaphrodite L4/adult, 128 neuron types). Needs pyreadr + pandas and network access:

    uv run --with pyreadr --with pandas python analysis/can_neuron/cengen_can.py

Data source: https://github.com/AlexWeinreb/cengenDataSC (data/cengen_sc_{1..4}.rda,
cengen_TPM.rda). Thresholds per Taylor et al. 2021 (bioRxiv 737577): "2" = recommended
compromise, "4" = most stringent (the level CIRCE's neuropeptide layer uses).

Expected output (recorded 2026-11-11, cengenDataSC master) is reproduced in RECORDED
below and in article.md; WBGene ids resolved via mygene.info (species 6239).
"""

from __future__ import annotations
import sys
import tempfile
import urllib.request
from pathlib import Path

BASE = "https://raw.githubusercontent.com/AlexWeinreb/cengenDataSC/master/data"

# gene -> WBGene id (mygene.info, C. elegans / taxid 6239)
PANEL = {
    "fast_synaptic": {
        "snt-1 (synaptotagmin-1, Ca2+ sensor)": "WBGene00004921",
        "unc-10 (RIM, active-zone)": "WBGene00006750",
        "unc-13 (priming)": "WBGene00006752",
        "rab-3 (SV GTPase)": "WBGene00004267",
        "snb-1 (synaptobrevin, v-SNARE)": "WBGene00004897",
    },
    "dense_core_peptidergic": {
        "unc-31 (CAPS, DCV exocytosis)": "WBGene00006767",
        "egl-3 (proprotein convertase)": "WBGene00001172",
        "egl-21 (carboxypeptidase)": "WBGene00001189",
        "ida-1 (DCV membrane)": "WBGene00002048",
        "snt-2": "WBGene00004922",
        "snt-4": "WBGene00004924",
    },
    "controls": {
        "eef-1A.1 (housekeeping)": "WBGene00001167",  # high in CAN -> data not dropout
    },
}

# Recorded CAN values (thr2 ON/OFF, thr4 ON/OFF, CAN TPM, percentile among 130 cell types,
# prevalence at thr2 across 128 neuron types). See article.md for interpretation.
RECORDED = """\
FAST SYNAPTIC            thr2 thr4   TPM  pctile  prev(thr2)
snt-1  (Ca2+ sensor)      OFF  OFF    0.5    1%    126/128   <- CAN = 1 of only 2 neuron types lacking it
unc-10 (RIM)             OFF  OFF    1.5    2%    118/128
unc-13 (priming)         ON   ON    52.2    9%    124/128
rab-3  (SV)              ON   ON    52.1    3%    128/128   (reporter: ABSENT in CAN, Hobert 2015)
snb-1  (v-SNARE)         ON   ON   150.3    5%    128/128
DENSE-CORE / PEPTIDE
unc-31 (CAPS)            OFF  OFF    0.7   12%     91/128
egl-3  (convertase)      OFF  OFF    4.4    1%    127/128
egl-21 (carboxypept.)    ON   ON    45.0    1%    128/128
ida-1  (DCV)             ON   ON    25.7    1%    128/128
snt-2                    OFF  OFF    0.0   24%     45/128
snt-4                    ON   ON   319.2    2%    128/128
CONTROL
eef-1A.1 (housekeeping)  ON    -   9254.0    -        -     (high -> CAN transcriptome well-sampled)
"""


def main() -> None:
    print(RECORDED)
    if "--fetch" not in sys.argv:
        print("[recorded values above] pass --fetch to re-download CeNGEN and recompute.")
        return
    import pandas as pd  # noqa: F401
    import pyreadr

    with tempfile.TemporaryDirectory() as td:

        def load(fn):
            p = Path(td) / fn
            urllib.request.urlretrieve(f"{BASE}/{fn}", p)
            res = pyreadr.read_r(str(p))
            return res[list(res.keys())[0]]

        sc2 = load("cengen_sc_2.rda")
        sc4 = load("cengen_sc_4.rda")
        tpm = load("cengen_TPM.rda")
        print(f"{'gene':<34}{'thr2':>5}{'thr4':>5}{'TPM_CAN':>9}{'pctile':>8}{'prev2':>9}")
        for group in PANEL.values():
            for sym, wb in group.items():
                b2 = int(sc2.loc[wb, "CAN"]) if wb in sc2.index else None
                b4 = int(sc4.loc[wb, "CAN"]) if wb in sc4.index else None
                p2 = int((sc2.loc[wb] > 0).sum()) if wb in sc2.index else None
                if wb in tpm.index:
                    v = tpm.loc[wb]
                    can = float(v["CAN"])
                    pct = round(100 * (v <= can).mean())
                    tv = f"{can:9.1f}{str(pct) + '%':>8}"
                else:
                    tv = f"{'NA':>9}{'NA':>8}"
                print(f"{sym:<34}{str(b2):>5}{str(b4):>5}{tv}{str(p2) + '/128':>9}")


if __name__ == "__main__":
    main()
