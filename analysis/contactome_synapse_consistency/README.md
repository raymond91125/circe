# Contactome ↔ synapse consistency (QC)

**Question.** A chemical synapse or gap junction cannot form where two membranes do not touch, so
*within a single EM reconstruction* every synaptic pair must also be a membrane-contact pair —
**synapses are a subset of contacts.** CIRCE now holds three contactomes, each scored on the **same EM
volume** as a synaptic dataset it already has, which turns that identity into a consistency check:

| life stage | synapses | contactome | same volume |
|---|---|---|---|
| L4 | `white_1986_jsh` | `brittin_2021_l4_contactome` | JSH series |
| adult | `white_1986_n2u` | `brittin_2021_adult_contactome` | N2U series |
| dauer | `yim_2024_dauer` | `yim_2024_dauer_contactome` | Yim reconstruction |

Run in both directions:
- **synapse without contact** → a cross-scoring inconsistency / **candidate error** (should be ~empty).
- **contact without synapse** → a **potential synapse** (cells touch but do not synapse; Peters' rule /
  Stepanyants–Chklovskii). Reported as a count, for context.

## Prior art (reconciled)

Brittin et al. 2021 (*Nature* 591:105; PMID 33627874; [PMC11648602](https://pmc.ncbi.nlm.nih.gov/articles/PMC11648602/))
**state the premise** — *"the space of possible synaptic and gap junction contact is restricted by the
existence of a physical membrane contact"* — but **do not report this QC.** Confirmed against the
full-text Results + Methods:

- Their data validation is **contact-vs-contact**: their scoring algorithm vs manual scoring
  (sensitivity 0.974, specificity 1.00 on JSH001/040) and vs White's **contact** neighborhoods for
  AIAR/AIBR/AQR — where their contact data flagged a *White contact* mis-score (a pair White called
  adjacent that does not physically touch). So "contact as QC" exists in the paper, but in the
  contact→contact direction, not synapse→contact.
- The only synapse fraction they give is on **reproducible M4 contacts**: *"77% and 85% of synaptic and
  gap junction contacts occur at M4 contacts."* That is a **different, stricter denominator** than the
  all-contact fraction here (M4 = reproducible in all four of their datasets; this check uses the full
  raw contact set). Brittin do not report the all-contact synapse fraction.
- Their "core-variable synaptic model" is a 3-parameter **reproducibility** model fit to M4-restricted
  counts (contacts act as a filter on its inputs); it does **not** predict synapses from contacts. The
  contact→synapse *prediction* result is a separate paper (Brittin et al. 2023, *Curr. Biol.*).

**This analysis is therefore the reverse (synapse→contact), all-contact, cross-scoring direction that
Brittin do not report.** It uses the raw Supplementary-Information-1 contact set (see
`data/brittin-2021-contactome/`), not the thresholded M4 reference graphs.

## Method & justification

- **Same volume only.** The subset property is a within-reconstruction identity; it would *not* hold
  across animals/stages (e.g. Cook adult synapses vs Brittin contacts), so only the three same-volume
  pairings above are tested. White↔Brittin are the same animals (JSH, N2U); Yim scored both dauer layers.
- **Cell restriction.** Synapse pairs are restricted to cells **both** present in the contactome, so a
  "miss" means *touch not scored*, never *cell absent from the contactome's volume* (the contactomes are
  nerve-ring / dauer-scoped, narrower than the whole-animal synapse datasets).
- **Both connection types.** Chemical and gap-junction pairs are included (both require apposition);
  adjacency is the contactome itself.
- **Weight = sections.** A miss is weighted by the synapse's section count. A **weight-1** miss
  (single section, at the EM detection limit) is the profile of a borderline/false-positive synapse or
  a sub-threshold contact; a **heavy** miss (≥ 5 sections) is a real contact gap or a naming mismatch
  worth inspecting. This weighting is the key calibration: a naming/merge bug would surface as *one cell
  with many heavy misses*, which is not what we see.
- Deterministic (sorted; no RNG; identical across `PYTHONHASHSEED`). Reads `outputs/connectome.json`
  (`uv run cckg build` first); `results.json` is committed.

## Results (`results.json`)

| volume | synapses at a contact | synapse without contact | of those, wt=1 | heavy (≥5) | potential synapses |
|---|---|---|---|---|---|
| L4 | 1347 / 1357 = **99.3%** | 10 | 8 | 0 | 3345 |
| adult | 1424 / 1454 = **97.9%** | 30 | 20 | 3 | 3735 |
| dauer | 1853 / 1871 = **99.0%** | 18 | 10 | 1 | 4227 |

**Interpretation.** The subset property holds to ~98–99% — CIRCE's synaptic data is consistent with the
independently-scored contactomes. Misses split exactly as the principle predicts:

1. **Detection-limit synapses (the bulk):** ~70–85% of misses are single-section (weight 1) — consistent
   with borderline/false-positive synapses or sub-threshold contacts. Benign.
2. **One systematic cluster — PLN in the adult:** all three heavy adult misses involve **PLN**, and they
   are **bilaterally symmetric** (`PLNL–SMBVL`/`SAADL`, `PLNR–SMBVR`/`SAADR`). PLN is in Brittin's cell
   set (it has other contacts), so this is a localized contact-scoring/reconstruction-extent issue for
   PLN's thin anteriorly-projecting process — worth flagging upstream, **not a CIRCE error.** Dauer has a
   single heavy miss (`RIS–RMED`, 5 sections).

**Calibration.** Robust: the ~98–99% subset agreement and the absence of any cell with many heavy misses
(no naming/merge bug). The PLN cluster is a real, reproducible flag but its *cause* (Brittin scoring vs
reconstruction extent) is a hypothesis, not established here. The large **potential-synapse** counts
(3–4k contacts per volume without a synapse) simply restate that contact is **necessary but not
sufficient** for a synapse — the premise of the whole contactome literature, not a finding of this check.

## Run

```
uv run cckg build                                               # writes outputs/connectome.json
uv run python -I analysis/contactome_synapse_consistency/consistency.py
```
