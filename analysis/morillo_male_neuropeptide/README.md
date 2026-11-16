# A predicted male long-range neuropeptide connectome (Phase C)

The first **male** neuropeptidergic connectome for *C. elegans*: a predicted "wireless" signaling
network derived in-build from the Morillo et al. 2025 male single-nucleus expression atlas, the male
counterpart to the Ripoll-Sánchez et al. 2023 *hermaphrodite* network. It is a KG dataset
(`morillo_2025_male_neuropeptide_lr`, `connection_type = neuropeptidergic`, `sex = male`), not a
standalone file — this doc records how it is built, why the method is trustworthy, what it covers,
and (importantly) what it deliberately does **not** claim.

## Headline

| | value |
|---|---|
| Edges (directed) | **11,586** |
| Neurons in network | 188 |
| Edge weight (NPP–GPCR pathways) | 1–9 (total 19,797) |
| Edges touching a male-specific neuron | **2,080** |
| Male-specific neurons participating | 17 (CEM×4, HOA, HOB, PHD×2, SPC/SPD/SPV×6, CA7, CP9, PVX) |
| NPP–GPCR pairs usable in male | 75 of 92 |

The 2,080 edges onto/from male-specific neurons (CEM, HOA/HOB, the spicule and ray-associated
neurons, …) are the part no hermaphrodite network can contain — the reason a male connectome is
worth deriving at all.

## Method — the *same* validated method as the hermaphrodite network

A directed edge *source → target* exists when, for one of the 92 biochemically deorphanized
NPP–GPCR pairs (Beets et al. 2023, EC50 ≤ 500 nM), the **source expresses the ligand** and the
**target expresses the cognate GPCR**; the edge weight is the number of such pairs. There is no
anatomical gating — every expressing source reaches every expressing target — which is exactly the
Ripoll-Sánchez **long-range** (LR) diffusion model.

The only thing swapped is the expression substrate: **Morillo 2025 male** (ligands from mmc8,
receptors from mmc12, binarized at the CeNGEN-calibrated cutoffs T = 0.2526 / 0.0491) in place of
Ripoll's CeNGEN hermaphrodite expression. The derivation runs in-build
(`assemble.py` → `ingest/neuropeptide.py:derive_neuropeptide_edges`) from the already-ingested male
`GeneExpression` records, so the connectome **regenerates from source** — there is no vendored male
adjacency matrix to drift out of sync.

### Why the engine is trustworthy (method validation)

`derive_neuropeptide_edges` is new code, so it is pinned against ground truth: fed CIRCE's **own**
Ripoll expression (`mechanistic/np_gene_expression.csv`) and the 92 pairs, it reproduces the
published Ripoll long-range network **exactly** — 53,558 edges, identical per-edge weights, zero
mismatches (unit test `test_derive_neuropeptide_edges_reproduces_published_herm_network`). So the
male connectome is the published, validated method applied to male data, not a bespoke algorithm.

## Prior art (checked first, per the repo's analysis conventions)

No published male neuropeptide connectome exists, so this is novel but not in conflict with anyone:

- **Ripoll-Sánchez et al. 2023** built the connectome, but **hermaphrodite only** (CeNGEN is
  hermaphrodite); it has no male expression and no male network.
- **Morillo et al. 2025** published the male expression **atlas** but **no connectome**.
- **Haque et al. 2025** (single-cell, sex-shared neurons) and **Haque et al. 2024** (developmental
  expression atlas) characterize male-biased neuropeptide *expression* but explicitly leave the male
  neuropeptide *connectome* to "future studies."

## External validation against Haque (expression ground truth)

The derivation can only be as good as the male expression layer feeding it, so that layer is
cross-checked against an independent dataset. Haque et al. 2024 (PMID 38769103) identify **`ins-39`**
as a prominent male-biased insulin-like peptide "expressed specifically in ciliated sensory
neurons," used as a male-isolation marker. CIRCE's ingested Morillo male layer recovers `ins-39` in
**AFDL/R, ASEL/R, ASJL/R, ASKL/R** (all ciliated amphid sensory) and in **zero** hermaphrodite
neurons — a direct, independent-platform match (Morillo snRNA vs. Haque bulk) on a gene with
published ground truth. (`ins-39` is not itself a deorphanized-pair ligand, so it validates the
expression substrate, not the edges.)

## Scope and caveats (calibrated confidence)

- **Long-range only.** The short- and mid-range Ripoll models gate edges by anatomical reach using a
  hand-curated hermaphrodite neuron × process-bundle matrix. No male equivalent exists, so SR/MR are
  **deferred**, not approximated. LR is the honest, fully-supported model for male.
- **Receptor precision.** Receptors come from snRNA at cluster level (mmc12), calibrated AUC 0.834;
  snRNA detects low-level receptor transcripts below CeNGEN's threshold-4, so individual **edges**
  (which require a receptor call on the target) are the least certain layer. Treat edge *existence*
  as a prediction, not an observation.
- **Binarized from continuous data** at a single calibrated cutoff; low-expression peptides are
  threshold-sensitive.
- **This is an *absolute* male resource, not a sex comparison.** Morillo's male sheet is sampled less
  deeply than its hermaphrodite sheet (~51 vs. ~75 ligands detected per neuron), a detection-depth
  confound that makes a male-vs-hermaphrodite edge-count comparison unreliable (it would spuriously
  make the male network look sparser). We therefore make **no remodeling / sex-difference claim**
  from this dataset. A defensible sex comparison needs a platform-matched both-sexes substrate (e.g.
  Haque 2025 single-cell) — a possible future "Phase D," kept out of scope here on purpose.

## Files

| File | What |
|------|------|
| `summarize.py` | Deterministic summary stats from `outputs/connectome.json` + source CSVs |
| `results.json` | Committed output of `summarize.py` |

## Reproduce

```
uv run cckg build                                   # derives the dataset into outputs/connectome.json
python analysis/morillo_male_neuropeptide/summarize.py
```

The derivation itself lives in the build (`src/celegans_connectome_kg/build/assemble.py`, the
"Male long-range neuropeptide connectome (DERIVED)" block) and its engine in
`src/celegans_connectome_kg/ingest/neuropeptide.py` (`derive_neuropeptide_edges`). Data provenance
for the male expression layer is in `data/morillo-2025-male/MANIFEST.md`.

## References

- Ripoll-Sánchez L, et al. *The neuropeptidergic connectome of C. elegans.* Neuron 111:3570–3589
  (2023). PMID 37935195.
- Beets I, et al. *System-wide mapping of peptide-GPCR interactions in C. elegans.* Cell Rep
  42:113058 (2023). PMID 37656621.
- Morillo KS, et al. *Single-nucleus neuronal transcriptional profiling of male C. elegans…* Cell
  Reports (2025). PMID 40682776.
- Haque R, et al. *Sex-specific developmental gene expression atlas unveils dimorphic gene networks
  in C. elegans.* Nat Commun 15:4273 (2024). PMID 38769103.
- Haque R, et al. *Decoding sexual dimorphism of the sex-shared nervous system at single-neuron
  resolution.* Sci Adv 11:eadv9106 (2025).
