# Morillo et al. 2025 — male neuropeptide expression

Per-neuron neuropeptide expression for **both sexes**, from the Morillo et al. 2025 adult male
single-nucleus RNA-seq atlas — the **male molecular counterpart** to the already-ingested Cook 2019
male *structural* connectome. Added to the KG as `Gene` + `GeneExpression` records, in two layers
kept at their native resolution (never pooled with CeNGEN) and scoped in
[`docs/proposals/morillo-2025-male-expression.md`](../../docs/proposals/morillo-2025-male-expression.md):

- **Phase A — ligands** (`morillo_2025_male`, `morillo_2025_hermaphrodite`): the FLP/NLP/INS
  neuropeptide ligands, **per neuron**, from supplement `mmc8`.
- **Phase B — receptors** (`morillo_2025_male_receptors`, `morillo_2025_hermaphrodite_receptors`):
  the NPP-GPCR receptors, **per cluster**, from the genome-wide supplement `mmc12`. Together the two
  layers let a male neuropeptide connectome be derived (Phase C, not yet built), paralleling the
  hermaphrodite Ripoll layer.

## Source

Morillo KS, St Ange J, Weng Y, Kaletsky R, Murphy CT. *Single-nucleus neuronal transcriptional
profiling of male* C. elegans *uncovers regulators of sex-specific and sex-shared behaviors.*
**Cell Reports** 2025. PMID 40682776 · PMC12908727 · doi:10.1016/j.celrep.2025.116016.

- **Ligands (Phase A):** `mmc8.xlsx`, sheets *"Herm/Male - Sexshared FlpNlpIns SCT"* —
  SCT-normalized mean expression, **neuron** column × gene row, for the sex-shared FLP/NLP/INS
  ligands.
- **Receptors (Phase B):** `mmc12.xlsx`, sheets *"HermThr1(0.5%)"* / *"MaleThr1(0.5%)"* —
  SCT-normalized mean expression, **cluster** column (row 1 = cluster number, row 2 = annotation) ×
  gene row, genome-wide; we read the NPP-GPCR receptor rows only.
- Neither supplement is **vendored** in this repo (large third-party downloads, treated as untrusted
  data). `build_expression.py` (mmc8) and `build_receptor_expression.py` (mmc12) read them from a
  path given on the command line or via `MORILLO_MMC8` / `MORILLO_MMC12`, run with `python -I`. The
  committed CSVs below are derived from them.

## Binarization cutoff (calibrated against CeNGEN, not arbitrary)

mmc8 holds **continuous** SCT values, so "cell X expresses gene G" needs an explicit cutoff. The
Morillo-herm value is thresholded at **T = 0.2526**, the cutoff that best reproduces (max F1)
**CeNGEN threshold-4** hermaphrodite calls — the same gold standard CIRCE's hermaphrodite
neuropeptide layer already uses.

- **Calibration set:** the 42 neuropeptide-ligand genes shared between `mmc8` and
  `../ripoll-2023-neuropeptide/mechanistic/np_genes.csv`, over the mmc8 columns that map 1:1 to a
  CeNGEN neuron/class (2,730 gene×neuron pairs, 12.2% CeNGEN-positive).
- **Fit:** AUC = **0.907**; best-F1 cutoff T = 0.2526 (precision 0.80, recall 0.64); prevalence
  CeNGEN 12.2% vs Morillo@T 9.7%.
- **Transfer:** the *same* T is applied to the male sheet — male and herm mmc8 share one platform,
  pipeline, and normalization, so T is a property of the assay, calibrated where a gold standard
  exists (herm) and transferred to the sex where none does (male). See the proposal doc for the
  mmc12 cross-check (ligand AUC 0.911 at T≈0.28) and the threshold-semantics caveat.

**Receptors (Phase B)** get their own cutoff — receptors are lower-expressed and noisier. The mmc12
herm value is thresholded at **T = 0.0491**, calibrated the same way over the 45 NPP-GPCRs and the
mmc12 clusters mapping 1:1 to a CeNGEN neuron/class (AUC = **0.834**; precision 0.47, recall 0.66;
CeNGEN prevalence 16.0% vs Morillo@T 22.8%) and transferred to male. The modest precision reflects
snRNA detecting low-level receptor transcripts below CeNGEN's stringent threshold-4 — a caveat for
any **edges** derived from this layer (Phase C), not for the expression records themselves.

## Neuron reconciliation (mmc8 column label → CIRCE cells)

mmc8 columns mix individual neurons (`HOA`), neuron classes (`AIY`, `CEM`), merged clusters the
snRNA could not resolve (`PHA/PHB/ASH`), generic groupings (`RnA`, `CP Neurons`, `Male 1`), and
non-neuronal cells (`Intestine`, `Sperm`). Policy (per the proposal's resolved curation choices):

- **Merged clusters dropped** — any label naming more than one class (`/` or `,`): the value is for
  the pooled group, not per neuron.
- **Non-neuronal / ungroundable groupings dropped** (`Intestine`, `Sperm`, `Hypodermis`, `RnA`,
  `RnB`, `CP Neurons`, `Male N`, `IL2DV`, `IL2LT`, `PVS`, …) — not individual CIRCE neurons.
- **Ventral-cord classes aliased** to CIRCE's `-n` class names (`VB`→`VBn`, `VC`→`VCn`, `VD`→`VDn`,
  `DA`→`DAn`, `DB`→`DBn`, `AS`→`ASn`).
- **Male-specific neurons included** where Morillo resolved them individually (17 cells: `CEM`×4,
  `HOA`, `HOB`, `PHD`×2, `SPC`/`SPD`/`SPV`×2, `CA7`, `CP9`, `PVX`).
- Kept labels are an individual cell name or a neuron class; the build (assemble) expands a class
  label to its member cells present in that sex (class-level calls, as with CeNGEN).

## Gene resolution

**Ligands.** `morillo_genes.csv` maps each mmc8 gene symbol to its WormBase gene id (`WB:WBGene…`),
resolved via MyGene.info (species 6239) and **cross-checked against `np_genes.csv` with zero
mismatches** for the 42 overlapping genes. All **133** mmc8 ligands are resolved. Two were not in
MyGene under their Morillo names and were resolved by WormBase curation to their canonical synonyms:
**`nlp-73` = `ctdh-1` (WBGene00013964)** and **`nlp-77` = `srlf-1` (WBGene00015514)**. The gene nodes
keep the Morillo symbol (the source's naming); the WBGene id is the authoritative key.

**Receptors.** `morillo_receptor_genes.csv` is the **45** NPP-GPCRs of the Beets/Ripoll pair set that
are present in mmc12, taken directly from `np_genes.csv` (already authoritative — no external lookup).
Six (`npr-37/39/40/41/42/43`) fall below mmc12's 0.5% gene-inclusion threshold and are absent from
the atlas, so they carry no Morillo receptor expression.

## Files

| File | What |
|------|------|
| `morillo_genes.csv` | symbol → WBGene for the 133 ligands (committed curation artifact) |
| `morillo_herm_expression.csv` / `morillo_male_expression.csv` | `cell,gene_symbol` ligand calls (mmc8) |
| `build_expression.py` | regenerates the ligand CSVs from `mmc8` + `morillo_genes.csv` |
| `morillo_receptor_genes.csv` | symbol → WBGene for the 45 NPP-GPCR receptors present in mmc12 |
| `morillo_herm_receptor_expression.csv` / `morillo_male_receptor_expression.csv` | `cell,gene_symbol` receptor calls (mmc12) |
| `build_receptor_expression.py` | regenerates the receptor CSVs from `mmc12` + `morillo_receptor_genes.csv` |

`cell` is a source column label (individual cell or neuron class, `-n`-aliased); the build expands
classes to member cells per sex. mmc12 is a combined atlas (all 145 clusters appear in every sheet),
so the build skips columns that belong to the other sex (e.g. male-specific `CEM` in the herm sheet).

## Reproduce

```
MORILLO_MMC8=/path/to/mmc8.xlsx   python -I data/morillo-2025-male/build_expression.py
MORILLO_MMC12=/path/to/mmc12.xlsx python -I data/morillo-2025-male/build_receptor_expression.py
uv run cckg build
```

Current output — **ligands:** 1,556 herm + 809 male per-cell records (172 / 132 cells);
**receptors:** 1,664 herm + 1,603 male per-cell records (169 / 178 cells; 17 male-specific each).

## Caveats

- **snRNA vs scRNA.** Morillo is single-*nucleus*; CeNGEN is single-*cell*. Male-vs-herm comparisons
  use Morillo's **own** herm sheet (never cross-wired against CeNGEN), so the two Morillo datasets
  are tagged distinctly and never pooled with the CeNGEN-derived layers.
- **Class-level calls.** A class column's call is assigned to every member cell; this is a class-
  level measurement, not a true per-cell one (the same convention as the CeNGEN-derived layers).
- **Binarized from continuous data at a single cutoff;** sex-bias of low-expression peptides is
  threshold-sensitive (see the proposal's threshold-semantics finding).
