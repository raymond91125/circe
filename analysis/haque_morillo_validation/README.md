# Validating the Morillo 2025 male neuropeptide layer against Haque 2025

CIRCE's male neuropeptide connectome (Phase C) rests entirely on the Morillo et al. 2025 male
single-**nucleus** expression layer. Haque et al. 2025 independently profiled the **same** sex-shared
nervous system by single-**cell** RNA-seq, in **both sexes**, and published per-neuron-class
expression — a platform-matched, independent yardstick. This analysis uses it two ways:

- **Part A — concordance:** does Haque corroborate the Morillo calls CIRCE ingested?
- **Part B — remodeling cross-check:** Haque is the platform-matched both-sexes dataset Morillo is
  not, so it can test the male-vs-hermaphrodite neuropeptide "remodeling" claim that Morillo's
  detection-depth confound forced us to drop from the connectome.

## Prior art (what Haque claims, which this reconciles with)

Haque et al. 2025 (*Sci Adv* eadv9106) state that "the male neuropeptide connectome undergoes
substantial remodeling, with most of the neuropeptides showing a male-bias expression" — an
**expression** claim; they build no connectome. Part B tests that claim directly from their tables;
Part A simply asks whether their expression agrees with Morillo's, gene by gene.

## Data & method

- **Haque tables** (not vendored — large untrusted third-party downloads): S7 `Herm.Exp.Mean` and S8
  `Male.Exp.Mean` (continuous mean expression per neuron class, per sex); S6 (binary calls, threshold
  0.04); S5 (the 62 high-confidence sex-shared clusters: ≥9 cells/sex, no mixed identities); S3
  (herm-vs-male differential expression per cluster, `avg_logFC > 0` = male-biased).
- **Genes:** the connectome's neuropeptide ligands (FLP/NLP/INS, `morillo_genes.csv`) and NPP-GPCR
  receptors (`morillo_receptor_genes.csv`), matched to Haque by symbol (two Morillo names via their
  canonical synonyms: `nlp-73`=`ctdh-1`, `nlp-77`=`srlf-1`).
- **Neuron classes:** Haque's class namespace; Morillo cells mapped to their CIRCE `cell_class`. The
  headline restricts to Haque's 62 high-confidence clusters; all annotated shared classes are a
  sensitivity check.
- **Metrics.** (A1) **AUC** — how well Haque's continuous per-class mean separates the classes
  Morillo called positive from those it called negative (Mann-Whitney; 0.5 = chance). AUC is the
  primary metric because it needs no shared threshold across two platforms. (A2) Cohen's **κ** and raw
  agreement between Morillo's binary calls and Haque's (S6). κ is reported but de-emphasised: most
  (gene, class) pairs are negatives, so the prevalence imbalance drives κ far below the raw agreement
  (the well-known κ paradox) — it is a floor, not the headline.

## Part A results — the Morillo layer is independently corroborated

Haque's 62 high-confidence sex-shared clusters:

| Layer | Sex | pooled AUC | mean per-gene AUC | binary agreement | κ |
|---|---|---|---|---|---|
| Ligands | male | **0.889** | 0.883 | 0.73 | 0.09 |
| Ligands | herm | **0.898** | 0.873 | 0.75 | 0.16 |
| Receptors | male | **0.791** | 0.831 | 0.76 | 0.34 |
| Receptors | herm | **0.774** | 0.795 | 0.77 | 0.35 |

- **Ligands are strongly validated** (AUC ≈ 0.89), consistent across both sexes — the layer the
  connectome's *sources* rest on holds up against an independent platform. This matches the original
  CeNGEN calibration (AUC 0.907), now reproduced by a third dataset.
- **Receptors are moderately validated** (AUC ≈ 0.78). Crucially, restricting to well-sampled
  clusters **raises** receptor AUC from ~0.70 (all shared classes, below) to ~0.78 — so the receptor
  layer's known softness **concentrates in poorly-sampled clusters**, and edges among well-sampled
  neurons are more trustworthy than the global figure implies.
- κ is higher for receptors (~0.34) than ligands (~0.09–0.16) precisely because receptors have more
  positive calls (less prevalence imbalance) — a property of the κ statistic, not of layer quality;
  the AUC ordering (ligands > receptors) is the real signal.

**Sensitivity — all annotated shared classes:** ligands AUC 0.865 (male) / 0.881 (herm); receptors
0.700 / 0.696. Ligands are stable; receptors drop on noisy clusters (see above).

**Discordant genes (AUC < 0.5, candidates for curation review, not assumed errors):** ligands
`nlp-52`, `flp-16`, `nlp-56` (male), `nlp-77`, `nlp-58`, `nlp-5`, `ins-16` (herm); receptors
`gnrr-6`, `dmsr-6` (male), `ntr-1` (herm). `nlp-77` (=`srlf-1`) is the strongest disagreement and
overlaps a Morillo symbol-synonym resolution — worth a closer look.

## Part B results — remodeling is real, but it is Haque's to claim, not Morillo's

Over the 62 high-confidence clusters:

- **Haque independently confirms net male-bias.** 86 of 133 neuropeptide ligands have higher mean
  expression in males than hermaphrodites (vs 47 herm-biased) — ~65%, i.e. "most," matching Haque's
  abstract. Their significant differential-expression events agree: **149 male-biased vs 108
  herm-biased** (padj < 0.05). The remodeling claim stands on platform-matched data.
- **Morillo cannot make that claim.** Morillo's per-gene sex-bias *direction* agrees with Haque's at
  only **49.2 %** (30 agree / 31 disagree) — chance. So Morillo's detection-depth confound does not
  merely inflate the hermaphrodite network globally; it **scrambles per-gene direction to noise**.

This is the quantitative vindication of the connectome's design decision: ship the male connectome as
an **absolute** resource (validated here), and source any **sex-difference / remodeling** statement
from Haque — never from a Morillo male-vs-herm comparison.

## Confidence calibration

- **Robust:** ligand-layer validation (AUC ≈ 0.89, three independent platforms agree); Haque's net
  male-bias (two independent measures, S7/S8 means and S3 DE); Morillo's unreliability for per-gene
  sex direction (chance-level agreement).
- **Method-dependent:** absolute receptor AUC depends on cluster-quality filtering (0.70 → 0.78); κ
  values depend on prevalence and are not comparable across layers. The discordant-gene list is a
  screen, not a verdict.
- **Out of scope:** male-specific neurons (CEM, HOA/HOB, spicule/ray) — Haque profiles sex-shared
  neurons only, so the connectome's male-specific edges are **not** validated here (they remain
  supported by Morillo + Cook alone).

## Files

| File | What |
|------|------|
| `validate.py` | Deterministic analysis (Parts A + B); reads Haque tables via `$HAQUE_DIR`, writes `results.json` |
| `results.json` | Committed output (so the numbers live in the repo without the large Haque files) |

## Reproduce

```
uv run cckg build      # produces outputs/connectome.json (used for the cell -> class map)
HAQUE_DIR=/path/to/PMID40644535 \
  uv run python -I analysis/haque_morillo_validation/validate.py
```

Haque supplementary tables S3/S5/S6/S7/S8 (`adv9106_table_s*.xlsx`) must be in `$HAQUE_DIR`; they are
intentionally **not** committed. The Morillo layer is documented in `data/morillo-2025-male/`, the
derived connectome in `analysis/morillo_male_neuropeptide/`.

## References

- Haque R, et al. *Decoding sexual dimorphism of the sex-shared nervous system at single-neuron
  resolution.* Sci Adv 11:eadv9106 (2025). PMID 40644535.
- Morillo KS, et al. *Single-nucleus neuronal transcriptional profiling of male C. elegans…* Cell
  Reports (2025). PMID 40682776.
