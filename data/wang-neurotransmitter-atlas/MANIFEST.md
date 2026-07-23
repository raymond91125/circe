# Pinned reference — Wang et al. C. elegans neurotransmitter atlas

Per-neuron neurotransmitter assignments from the comprehensive CRISPR-reporter atlas — a
second **cited neurotransmitter source** (alongside [`../wormatlas-neurotransmitter/`]
(../wormatlas-neurotransmitter/MANIFEST.md)) toward the `nt` provenance gap noted in
[`../neuron-graph/MANIFEST.md`](../neuron-graph/MANIFEST.md).

## Source

- **Paper:** Wang C, Vidal B, Sural S, Loer C, Aguilar GR, Merritt DM, et al. *A neurotransmitter
  atlas of C. elegans males and hermaphrodites.* eLife 13:RP95402 (**version of record**), 2024.
- **DOI:** 10.7554/eLife.95402 · https://elifesciences.org/articles/95402
- **Crawled (UTC):** 2026-07-03T00:13:05Z
- **Note:** the published article page is bot-blocked (HTTP 403 to headless fetch); the
  supplementary files were retrieved from the eLife CDN (`cdn.elifesciences.org/articles/95402/`).
  The VOR **Table S2 is content-identical** to the earlier reviewed-preprint v2 supplement
  (302 assignments, no changes) — only the provenance here points to the version of record.
- **Lab:** Hobert lab. Assignments are based on CRISPR-engineered reporter alleles for the
  transmitter-pathway genes (eat-4/VGLUT → glutamate, unc-17/cha-1 → ACh, unc-25/unc-47 →
  GABA, cat-2/bas-1 → dopamine, tph-1/cat-1 → serotonin, tdc-1/tbh-1 → tyramine/octopamine),
  cross-checked against scRNA-seq (CeNGEN) and prior antibody data.

## Files

| File | What |
|------|------|
| `TableS2_expression.xlsx` | pinned supplement (`elife-95402-supp2-v1.xlsx`) — per-neuron **hermaphrodite** reporter table + consolidated "Neurotransmitter(s)" call |
| `parse_atlas.py` | reproducible parser (`uv run --with openpyxl python parse_atlas.py`) |
| `assignments.csv` | parsed **302** hermaphrodite per-neuron assignments (`class, neuron, neurotransmitters`) — a reconciliation reference for `Cell.neurotransmitter` |
| `Supp3_male_expression.xlsx` | pinned supplement (`elife-95402-supp3-v1.xlsx`) — the **male** reporter atlas (male-specific neurons) |
| `Supp4_dimorphic.xlsx` | pinned supplement (`elife-95402-supp4-v1.xlsx`) — **sexually-dimorphic** sex-shared neurons (hermaphrodite vs male) |
| `parse_male_atlas.py` | parser for Supp3 + Supp4 → `sex_neurotransmitters.csv` |
| `sex_neurotransmitters.csv` | **116** per-sex assignments (`cell, sex, neurotransmitter, confidence, note`) — the build input for the reified `NeurotransmitterAssignment` records |

## Male / per-sex assignments (`sex_neurotransmitters.csv`)

`parse_male_atlas.py` builds the per-sex assignments the KG ingests:
- **Supp3** → the 90 male-specific neurons (mechanical map of the consolidated call to the
  neuron-graph code scheme `a`/`l`/`g`/`d`/`s`/`u`; `*`/`?` → confidence `putative`). Grouped
  serial names (`DX1/2`, `EF3/4`) are expanded to member cells and de-duplicated.
- **Supp4** → 8 sex-shared neurons with a transmitter-**identity** difference by sex, curated to
  both a hermaphrodite and a male code (e.g. `AIM` Glu→ACh; several gain GABA via `unc-47`). AVG
  is excluded (its difference is `unc-17` expression level only, cholinergic in both sexes).
- `CP0`, `DX4`, `EF4` appear in Supp3 but are absent from the Cook connectome. They are minted as
  grounded, connectivity-free male neurons via `data/curation/atlas_only_cells.csv` so their atlas
  neurotransmitter attaches. The unqualified `Cell.neurotransmitter` (hermaphrodite/neuron-graph)
  is left untouched.

`TableS2_expression.xlsx` (= `elife-95402-supp2-v1.xlsx`) sha256:
`ad5175b5e53748acf77c6dfff1b5efcf15ed4364091239bd2def0095544eb18c`

Other VOR supplements (not vendored), at `https://cdn.elifesciences.org/articles/95402/elife-95402-<supp>-v1.xlsx`:
`supp1` (reporter-detail atlas) `bc472336…`, `supp3` `34e8fddd…`, `supp4` `d69a4d78…`,
`supp5` `fce87f65…`.

## Parse notes

- The consolidated **"Neurotransmitter(s)"** column is used (the per-reporter columns encode
  expression by cell fill colour, which isn't recoverable from cell values).
- Values: `ACh`, `Glu`, `GABA`, `DA`, `5-HT`, `octopamine`, `betaine`, and `unknown (orphan)`;
  `*` marks assignments updated by this study, `- NEW` newly identified, `(uptake)` = uptake
  (not synthesis). These markers are preserved verbatim in `assignments.csv`.
- Covers both sexes; per-neuron rows are the hermaphrodite unless noted. For neurons with a
  sex- or stage-specific transmitter switch (e.g. AIM), the paper's text is the authority.

Attribution: cite Wang et al. 2024 (eLife) when using these assignments.
