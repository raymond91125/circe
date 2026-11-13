# Ingestion assessment — Morillo et al. 2025 male neuronal expression

*Status: proposal / scoping (no schema or data changes yet). Drafted 2026-11-12; calibration
validated against the real supplements 2026-11-13.*

## Paper

Morillo KS, St Ange J, Weng Y, Kaletsky R, Murphy CT. **Single-nucleus neuronal
transcriptional profiling of male *C. elegans* uncovers regulators of sex-specific and
sex-shared behaviors.** *Cell Reports* 2025. PMID 40682776 · PMC12908727 ·
doi:10.1016/j.celrep.2025.116016.

Adult **male** single-nucleus RNA-seq (snRNA-seq) neuronal atlas, with a matched
hermaphrodite dataset from the same platform. Maps dimorphic GPCR, neuropeptide, and
ion-channel expression; the male companion to the St. Ange & Weng hermaphrodite atlas.

## Why it fits CIRCE

CIRCE's [sex extension](../PLANNING.md#sex-extension--male--hermaphrodite-unified-sex-aware-graph)
already ingests the Cook 2019 **male structural connectome** (datasets `cook_2019_male` /
`cook_2019_hermaphrodite`, tagged by `Dataset.sex`), and male-specific cells are grounded to
male WBbt terms. But every **expression**-driven layer is still hermaphrodite-only:

- Gene expression / `GeneExpression` is CeNGEN adult hermaphrodite (threshold-4).
- The neuropeptide connectome (Ripoll 2023) and monoamine connectome layers derive their
  edges from that hermaphrodite expression.

So CIRCE has male *wiring* but no male *molecular* layer. Morillo 2025 is the most direct
source to fill that gap — it is the **expression counterpart to the already-ingested male
structural connectome**, and would enable male neuropeptide / monoamine connectome
projections parallel to the existing hermaphrodite ones.

**Schema impact: none required.** `Dataset` already carries `sex` + `life_stage`, and
`GeneExpression` is reified against a `dataset`. Male expression ingests as `GeneExpression`
records tied to a new `morillo_2025_male` dataset (`sex: male`, `life_stage: adult`) — exactly
the Cook 2019 pattern.

## Supplement inventory (verified from the downloaded files)

| File | Content | Structure | Ingestibility |
|---|---|---|---|
| **mmc12** (103 MB) | Per-neuron expression matrix, **Male + Herm**, 4 thresholds (0.5/1/1.5/3%), genome-wide (~17k genes) | 145 **clusters** × genes; values = SCT-normalized mean | ⭐ richest, but cluster-level (see caveats) |
| **mmc8** | Sex-shared **neuropeptide (FLP/NLP/INS)** expression by neuron, male & herm (SCT + PCT) | **Per individual neuron class**; % + normalized | ⭐ cleanest for the neuropeptide layer |
| mmc7 | csGPCR (chemosensory GPCR) expression, male/herm/sex-shared | per-neuron | receptors, but *chemosensory* — not the NPP-GPCRs the Beets/Ripoll pairs use |
| mmc3 (39 MB) | Expression matrix, all thresholds (single sex) | cluster-level | superseded by mmc12 |
| mmc2 / mmc9 / mmc13 | Cluster annotations, per-cluster markers, male-specific neuron list | — | male-specific neuron molecular identity; WBBT grounding aid |
| mmc5 / mmc6 / mmc10 / mmc11 / mmc4 | Male/herm-enriched genes, lectins, sex-DE per neuron, transcriptomic distances, spicule markers | analysis outputs | reference, not KG layers |

## Caveats (decisions this forces)

1. **Cluster vs individual-cell resolution (mmc12).** Columns are 145 clusters annotated to
   neuron classes, but **many are merged** (`VA/VB/DA/DB`, `CEP.ADE.PDE`, `CEM/Ray1,5,7B`) and
   some are **non-neuronal** (intestine, sperm, hypodermis, pharynx, muscle). CIRCE keys on
   individual cells; merged clusters cannot be cleanly split, and non-neuronal clusters must be
   filtered. → Ingest only clusters that map 1:1 (or 1:n with explicit fan-out) to CIRCE cells,
   and record the cluster→cell mapping as curation.
2. **Continuous values, not binary.** The "threshold" sheets filter *genes* by % of cells
   expressing; cells still hold continuous normalized values. Defining "cell X expresses gene G"
   needs an explicit per-cell binarization rule. → Resolved by **CeNGEN calibration** (see below),
   not an arbitrary cutoff.
3. **Platform mismatch with CeNGEN.** Morillo is single-**nucleus**; CeNGEN is single-**cell**.
   Male-vs-herm comparisons must use Morillo's **own** herm sheets (mmc12 `HermThr*`, or mmc8
   herm), **not** cross-wired against CeNGEN, or sex and platform confound. → Tag the Morillo
   datasets distinctly; never silently pool with CeNGEN.
4. **Gene keys are symbols** → resolve to WBGene (CIRCE's existing Alliance/expression-alignment
   machinery; same approach as the Ripoll layer).

## Binarization by CeNGEN calibration (chosen approach)

Rather than pick a Morillo cutoff arbitrarily, calibrate it against CeNGEN's established
hermaphrodite calls — the same gold standard CIRCE's hermaphrodite layer already uses — then
transfer the calibrated rule to the male data:

1. **Gold standard:** CeNGEN **threshold-4** binary calls for hermaphrodite neurons (the level
   CIRCE's neuropeptide layer uses) — available as the already-vendored
   `data/.../30072020_CENGEN_threshold4_*_allneurons.csv` (NPP/NPR/MR/LGC gene set) or the full
   cengenDataSC `cengen_sc_4` matrix.
2. **Match space:** join on `(neuron class, WBGene)` over the neurons that map **1:1** between
   CeNGEN's classes and Morillo's clusters, restricted to the **connectome-relevant gene set**
   (neuropeptides + their receptors) so the cutoff is tuned on the expression regime it will
   actually be applied to.
3. **Fit:** find the Morillo-herm cutoff *T* on the continuous value (SCT mean, or the PCT
   detection fraction) that best reproduces the CeNGEN threshold-4 calls — maximizing F1 /
   Youden's J, with matched prevalence as a cross-check. (If PCT is used, this effectively asks
   *which Morillo % threshold equals CeNGEN threshold-4*.)
4. **Transfer:** apply the *same* *T* to Morillo **male**. Sound because male and herm Morillo data
   share one platform, pipeline, and normalization — *T* is a property of the assay, calibrated on
   the sex where a gold standard exists and used on the sex where none does.

## Calibration — validation results (run against the real supplements)

The calibration above was executed against the vendored CeNGEN threshold-4 CSV, the
`np_genes.csv` symbol→WBGene bridge, and the downloaded supplements (parsed in isolation —
Python `-I`, paths passed as args, scripts kept out of the download dir). SCT-normalized mean was
used as the continuous value; Morillo's PCT detection fraction is **not** carried in the mmc12
threshold sheets, so the SCT-vs-PCT sub-choice (open question 1) resolves to **SCT**.

**mmc12 full atlas (HermThr1 sheet, 1:1-mappable clusters only).** The SCT value separates
CeNGEN-positive from CeNGEN-negative calls well — but **ligands and receptors need separate
cutoffs**:

| gene class | AUC | best-F1 cutoff *T* |
|---|---|---|
| neuropeptide ligands | **0.911** | *T* ≈ 0.28 |
| neuropeptide receptors (GPCRs) | **0.864** | *T* ≈ 0.07 |

(mmc8's per-individual-neuron FLP/NLP/INS SCT gives a consistent ligand AUC ≈ 0.90, confirming
cluster resolution does not degrade the ligand calibration.) So cluster-level resolution does
**not** hurt ligands; the single-*T* assumption does — receptors sit an order of magnitude lower
and must be binarized at their own cutoff.

**Male transfer — dimorphic-case check.** Applying the paper's own detection definition
(`value > 0` in the threshold sheet) to the transferred male sheet reproduces **all four** tested
dimorphic calls exactly (herm/male cluster columns verified identical, 145 each):

| gene | paper | Morillo @ `value>0` |
|---|---|---|
| flp-34 | male-biased | male 97 > herm 84 ✓ |
| nlp-43 | male-biased | male 97 > herm 80 ✓ |
| flp-1 | herm-biased | herm 92 > male 69 ✓ |
| flp-8 | herm-biased | herm 95 > male 70 ✓ |

**Threshold semantics — a finding that constrains ingestion.** Sex-bias *direction* is
threshold-dependent. flp-34 is male-biased only in the broad tail of **low-level** detections; at
its few strong sites it is herm-biased (DVC herm 13.3 / male 7.3; PVW herm 2.3 / male 0.56). So at
the stringent CeNGEN-calibrated cutoff (*T* ≈ 0.28) the sex-bias of low-expression peptides
**weakens or flips** relative to the paper's permissive `value>0` calls. These are two different
questions, not a contradiction:

- **Connectome layer → use the calibrated *T*.** It puts the male expression layer at the *same
  detection stringency* as CIRCE's existing CeNGEN hermaphrodite layer, so the male and herm
  connectomes are legitimately comparable (CIRCE's core value). Cost: the fine sex-dimorphism of
  low-expression peptides falls below the bar and is not captured as edges.
- **Dimorphism annotation → use the paper's `value>0`.** Preserves the paper's dimorphic calls,
  but is a dense, permissive regime (flp-34 "expressed" in ~90 of ~130 neuronal clusters) that is
  not comparable to CeNGEN and must not drive connectome edges.

**Recommendation:** binarize connectome edges at the calibrated *T* (separate ligand / receptor
cutoffs); if the dimorphism signal is wanted, surface Morillo's `value>0` sex-bias as a
**separate, clearly tagged annotation**, never as connectome edges. Do not serve both from one
threshold — and record in the dataset metadata that CIRCE's calibrated calls can differ from (even
oppose) the paper's `value>0` calls for low-expression genes.

**Residual caveats:** calibration is valid only on the 1:1-mappable neuron subset; merged Morillo
clusters inherit *T* but their cell-level calls stay coarse. The transfer assumption
(herm-calibrated *T* applies to male) is sound in principle (shared platform, pipeline, and
normalization) and passes the dimorphic-case check above, but remains an assumption — male
prevalence stays a guardrail.

## Proposed plan (phased; mirrors the sex-extension milestones)

**Phase A — focused neuropeptide layer (mmc8), lowest risk**
- M1 Pin `mmc8` (+ paper metadata, PMID/DOI) under `data/morillo-2025-male/` with MANIFEST +
  sha256; downloaded files treated as untrusted during parsing.
- M2 Define datasets `morillo_2025_male` (`sex: male`) and `morillo_2025_hermaphrodite`
  (`sex: hermaphrodite`), `life_stage: adult`.
- M3 Parse per-neuron FLP/NLP/INS expression → `GeneExpression` records, binarized at the
  CeNGEN-calibrated cutoff *T* (see above); resolve gene symbols → WBGene.
- M4 Reconcile neuron names against the CIRCE registry; curation work-list for merged/ambiguous.
- M5 Build + export; SPARQL verify (male neuropeptide-ligand expression present; sex partitions
  cleanly; spot-check known dimorphic cases — e.g. `nlp-43`/`flp-34` male-biased,
  `flp-1`/`flp-8`/`ins-30` herm-biased, per the paper).

**Phase B — full male expression atlas (mmc12)** — adds the genome-wide male matrix (incl. the
NPP-receptor GPCRs the connectome needs), resolving the cluster→cell mapping and binarization
from Phase A. Use the **receptor-specific** cutoff (*T* ≈ 0.07) for GPCRs, not the ligand cutoff
(validated above).

**Phase C — male neuropeptide / monoamine connectome** — combine the male expression layer with
the Beets NPP-GPCR pairs (already vendored) to project a **male** neuropeptide connectome
mirroring the hermaphrodite Ripoll layer. *Requires a prior-art check first* (per
[CLAUDE.md](../../CLAUDE.md) analysis conventions) — confirm whether a male neuropeptide
connectome has been published before presenting a CIRCE-derived one as novel.

## Open questions

1. Binarization — **resolved**: CeNGEN-calibrated cutoff on the SCT value, with **separate ligand
   and receptor cutoffs** (validated above, AUC 0.91 / 0.86; PCT is not carried in the mmc12
   threshold sheets, so SCT is used).
2. Phase scope — stop at the neuropeptide layer (A), the full atlas (B), or go to a male
   connectome (C)?
3. For dimorphic/merged clusters that can't map to individual cells, drop them or record them as
   a coarse cluster-level annotation outside `GeneExpression`?
4. Keep Morillo's hermaphrodite sheets too (enables within-platform sex comparison) or male-only?
