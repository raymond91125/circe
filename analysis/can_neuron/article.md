# Is CAN a chemical-transmitting neuron? Three lines of evidence

*Analysis date: 2026-11-11. Data: CIRCE `outputs/connectome.json`; CeNGEN scRNA-seq
(Taylor et al. 2021); Stefanakis, Carrera & Hobert 2015 (Neuron).*

## Question

The **CAN** pair (canal-associated neurons, CANL/CANR) is a bona-fide neuron by lineage,
yet CIRCE types it `cell_type: other` and it falls outside the "300 hermaphrodite neurons"
count (300 + CANL/CANR = the textbook 302). Does CAN actually function as a *chemical-
transmitting* neuron — i.e. does it send chemical synapses — or is it something else?

The discriminating criterion is **presynaptic chemical output**. Gap junctions, receptor
expression, and receiving chemical input are all shared with non-neuronal cells (glia,
muscle, hypodermis) and do not establish neuronal chemical transmission. The defining act
is releasing neurotransmitter onto targets.

## Finding

**Three independent datasets converge: CAN does not do conventional fast chemical synaptic
transmission.**

### 1. Connectome (EM reconstructions) — `can_circe.py` → `results.json`

CAN is **not** synapse-free (a common misconception from the classical White et al. 1986
view), but its only structural partners are:

| Type | Partners (deduped) | Datasets |
|------|--------------------|----------|
| Chemical **input** | ALA → CANL, ALA → CANR | Cook 2019 herm + male; Witvliet 2020 (ALA→CANR) |
| Chemical **output** | **NONE** | — (all datasets) |
| Gap junction | ALA↔CAN, CAN↔PVD, CAN↔excretory cell, CAN↔hypodermis; CANR↔CEPsh (glia) | Cook 2019 herm + male; Witvliet 2020 |

CAN appears **only as a postsynaptic (receiving) partner** — it has **zero chemical-synapse
outputs** in every EM dataset CIRCE holds. Its partners are telling: ALA, PVD, and three
*non-neuronal* cells (excretory cell, hypodermis, CEP sheath glia), consistent with an
osmoregulatory / canal role. (CIRCE also carries 1,466 CAN edges in the *predicted*
neuropeptide layer; those are diffusion-model inferences from expression, not observed
synapses.)

Asymmetry worth noting: in the same volumes where no *output* synapse was found, an *input*
synapse onto CAN **was** resolved — so "outputs just below detection" is weaker than a bare
unreconstructed-region argument (though not excluded; presynaptic specializations can be small).

### 2. Transcriptome (CeNGEN) — `cengen_can.py`

CIRCE's gene layer (155 signaling genes) does not include release machinery, so this comes
from CeNGEN. CAN = merged CANL/CANR class, 128 neuron types, percentile among 130 cell types:

| gene | thr-2 | thr-4 | CAN TPM | pctile | in N/128 |
|------|:---:|:---:|---:|---:|---:|
| **snt-1** (synaptotagmin-1, Ca²⁺ sensor) | **OFF** | OFF | 0.5 | 1st | 126/128 |
| **unc-10** (RIM, active zone) | **OFF** | OFF | 1.5 | 2nd | 118/128 |
| unc-13 (priming) | ON | ON | 52 | 9th | 124/128 |
| rab-3 (SV) | ON | ON | 52 | 3rd | 128/128 |
| snb-1 (v-SNARE) | ON | ON | 150 | 5th | 128/128 |
| unc-31 (CAPS, DCV release) | OFF | OFF | 0.7 | 12th | 91/128 |
| egl-3 (convertase) | OFF | OFF | 4.4 | 1st | 127/128 |
| egl-21 / ida-1 / snt-4 (DCV) | ON | ON | 45 / 26 / 319 | 1st–2nd | 128/128 |
| *eef-1A.1 (housekeeping control)* | ON | — | **9254** | — | — |

Reading:

- CAN **lacks the two genes that define fast, Ca²⁺-triggered, active-zone synaptic release —
  `snt-1` and `unc-10`** — at CeNGEN's *permissive* threshold. For `snt-1`, CAN is literally
  **one of only 2 of 128 neuron types** without it.
- The generic vesicle toolkit (`rab-3`, `snb-1`, `unc-13`) is present but at the **bottom of
  the neuronal range** (3rd–9th percentile).
- The dense-core/peptidergic program is **also weak** (mostly 1st–12th percentile; `unc-31`
  and `egl-3` off). So CAN is not a high-output secretory cell of *either* kind.
- The housekeeping gene is high (9254 TPM) → CAN's transcriptome is **well-sampled**, so the
  low/zero release-gene values are real biology, **not dropout**.

### 3. Reporter transgenics (Stefanakis, Carrera & Hobert 2015, *Neuron*)

The canonical `rab-3` pan-neuronal reporter (`rab-3prom1`) is expressed in **300/302 (99%)**
adult neurons. Quoting: *"The only neurons in which we did not observe rab-3 expression are
the canal-associated neurons (CAN), a neuron pair that was previously noted for its scarcity
of synaptic connections with other neurons."* CAN is **the single exception** across the
nervous system, explicitly linked to its synaptic scarcity. (The 2 missing neurons are
CANL/CANR — the same pair behind the 300-vs-302 count.)

## A reconciled discrepancy (reporter vs scRNA-seq for `rab-3`)

The two molecular methods disagree on the binary call, and the disagreement is informative:

| | Hobert 2015 | CeNGEN 2021 |
|---|---|---|
| Method | `rab-3prom1` **reporter** (4.3 kb promoter + intron 1) | endogenous **transcript** (scRNA-seq) |
| CAN call | **absent** (sole exception of 302) | **present, 3rd percentile** (lowest tier) |

Resolution = the Hobert paper's own thesis: pan-neuronal genes are driven by **multiple,
redundant, distributed cis-regulatory elements**, no single master element, so any one
reporter fragment can miss a neuron. The `rab-3prom1` fragment lacks whatever drives CAN's
faint expression; scRNA-seq still detects the low endogenous transcript. **Both agree CAN is
the `rab-3` outlier**; they differ only because a promoter fragment ≠ the endogenous gene.

## Conclusion & confidence

- **Robust (3 independent lines):** CAN does **not** perform conventional fast chemical
  synaptic transmission — no EM chemical output, no `snt-1`/`unc-10`, and it is the lone
  neuron missing the `rab-3` reporter.
- **Robust:** CAN *is* electrically coupled (gap junctions, innexins `inx-7`/`inx-13` in
  CIRCE) and *receives* chemical input from ALA. It is not synapse-free.
- **Moderate:** CAN is not a strong peptidergic/volume transmitter either — its DCV
  machinery is weakly expressed (`unc-31`/`egl-3` off). An earlier framing of CAN as a
  dedicated "peptide broadcaster" is **only partly supported** and is corrected here.
- **Interpretation (lower confidence):** the molecular profile fits CAN's known essential
  role in **osmoregulation / excretory-canal fluid homeostasis** rather than classical
  neurotransmission, which is consistent with its atypical, near-non-neuronal identity and
  with CIRCE typing it `other`.

## Caveats

- CeNGEN merges CANL/CANR and is hermaphrodite L4/adult; it cannot speak to male-specific or
  developmental states. CIRCE records no sex-specific override for CAN's transmitter.
- This paper's CAN-specific statement concerns `rab-3` (and the reporter panel); it makes
  **no** CAN-specific claim about `snt-1` (which appears only as an smFISH probe in methods).
  The `snt-1`-absent-in-CAN result rests on CeNGEN alone.
- "No output synapse" is bounded by EM detection limits (see asymmetry note above).

## References

- Cook SJ et al. 2019. Whole-animal connectomes of both *C. elegans* sexes. *Nature* 571:63.
- Witvliet D et al. 2021. Connectomes across development. *Nature* 596:257.
- Taylor SR et al. 2021. Molecular topography of an entire nervous system (CeNGEN).
  *Cell* 184:4329 (bioRxiv 737577). Data: github.com/AlexWeinreb/cengenDataSC.
- Stefanakis N, Carrera I, Hobert O. 2015. Regulatory logic of pan-neuronal gene expression
  in *C. elegans*. *Neuron* 87:733. PMC4545498.
- White JG et al. 1986. The structure of the nervous system of *C. elegans*.
  *Phil Trans R Soc Lond B* 314:1.
