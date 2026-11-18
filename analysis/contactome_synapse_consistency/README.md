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

### Scope caveat — read this before the error interpretation

The contactomes cover only the **nerve ring**; the White synapse datasets from the same series extend
**beyond** it (sublateral / ventral cord, etc.). So a White↔Brittin "synapse without contact" has a
cause besides scoring error: the **synapse was made outside the nerve-ring volume Brittin segmented for
contacts.** Restricting to cells both present in the contactome is *necessary but not sufficient* — a
nerve-ring cell can still synapse elsewhere. Consequences:

- **Dauer is the clean control, the only one without this confound.** Yim scored synapses *and* contacts
  from the **same reconstruction** (identical volume *and* scope), so its **99.0%** is the one clean
  subset estimate. The two White↔Brittin fractions are **lower bounds** (misses inflated by the scope
  mismatch).
- **The cause of an individual White↔Brittin miss is NOT resolvable here.** CIRCE carries no per-synapse
  coordinates, so we cannot tell an out-of-nerve-ring synapse from a Brittin contact-scoring gap from a
  White false-positive. The analysis records, as **context only**, whether a miss pair is scored as a
  contact in *another* contactome — but that is a **different animal and stage**, neurites remodel, and
  "apposed elsewhere" is equally consistent with a scoring gap as with an out-of-volume synapse, so it
  does **not** establish cause. It only shows the cells are not categorically non-adjacent. (~half of
  misses are such pairs: adult 14/30, L4 5/10, dauer 8/18.) Do not read this field as an out-of-volume
  attribution — an earlier draft did, wrongly.

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
  a sub-threshold contact; a **heavy** miss (≥ 5 sections) is a substantive discrepancy worth inspecting.
  This weighting is the key calibration: a naming/merge bug would surface as *one cell with many heavy
  misses*, which is not what we see.
- Deterministic (sorted; no RNG; identical across `PYTHONHASHSEED`). Reads `outputs/connectome.json`
  (`uv run cckg build` first); `results.json` is committed.

## Results (`results.json`)

| volume | synapses at a contact | synapse w/o contact | of those, wt=1 | heavy (≥5) | apposed in another stage† | potential synapses |
|---|---|---|---|---|---|---|
| L4 | 1347 / 1357 = **99.3%** | 10 | 8 | 0 | 5 | 3345 |
| adult | 1424 / 1454 = **97.9%** | 30 | 20 | 3 | 14 | 3735 |
| dauer (clean control) | 1853 / 1871 = **99.0%** | 18 | 10 | 1 | 8 | 4227 |

† context only — a different animal/stage; does **not** attribute the miss (see scope caveat).

**Interpretation.** The subset property holds to ~98–99% — CIRCE's synaptic data is consistent with the
independently-scored contactomes. Reading the misses:

1. **Detection-limit synapses:** ~55–80% of misses are single-section (weight 1) — the profile of
   borderline/false-positive synapses or sub-threshold contacts. Benign.
2. **Scope (White↔Brittin only):** the White datasets extend beyond the nerve ring, so some misses are
   synapses made outside Brittin's volume. This cannot be confirmed per-miss (no synapse coordinates),
   but it makes the two White↔Brittin fractions **lower bounds**. The **dauer same-scope control at
   99.0%** has no such confound and bounds the within-reconstruction inconsistency at ~1%.
3. **The one systematic cluster — PLN in the adult.** All three heavy adult misses involve **PLN**,
   bilaterally symmetric (`PLNL–SMBVL`/`SAADL`, `PLNR–SMBVR`/`SAADR`). Checking the data: these pairs
   **are scored as nerve-ring contacts in the L4 and dauer reconstructions, but not in the adult (N2U)**,
   and PLN has fewer adult contacts (PLNL 21 / PLNR 14) than in L4 (27 / 25) or dauer (31 / 23). So this
   is a **within-N2U discrepancy** (White's N2U synapse vs Brittin's N2U no-contact), *not* out-of-volume
   (the appositions sit in the nerve ring in other stages). The most likely cause is an **incomplete
   Brittin adult contact-scoring of PLN's process**; a White adult false-positive or real adult
   remodeling can't be excluded without the EM. Either way — a flagged upstream discrepancy, **not a
   CIRCE error.** *(An earlier draft called this out-of-nerve-ring; the data refute that.)*

**Calibration.** Robust: the subset agreement, and the absence of any cell with many heavy misses (a
genuine naming/merge bug would surface that way). **Not resolved:** the cause of individual White↔Brittin
misses — those fractions are lower bounds (scope confound); only the dauer same-scope control (99.0%) is
a clean within-reconstruction estimate. The "apposed in another stage" field is context, not a cause
attribution (cross-stage, non-transferable). The large **potential-synapse** counts (3–4k contacts per
volume without a synapse) simply restate that contact is **necessary but not sufficient** for a synapse —
the premise of the contactome literature,
not a finding of this check.

## Run

```
uv run cckg build                                               # writes outputs/connectome.json
uv run python -I analysis/contactome_synapse_consistency/consistency.py
```
