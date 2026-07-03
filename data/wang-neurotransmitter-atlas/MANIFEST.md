# Pinned reference — Wang et al. C. elegans neurotransmitter atlas

Per-neuron neurotransmitter assignments from the comprehensive CRISPR-reporter atlas — a
second **cited neurotransmitter source** (alongside [`../wormatlas-neurotransmitter/`]
(../wormatlas-neurotransmitter/MANIFEST.md)) toward the `nt` provenance gap noted in
[`../neuron-graph/MANIFEST.md`](../neuron-graph/MANIFEST.md).

## Source

- **Paper:** Wang C, Vidal B, Sural S, Loer C, Aguilar GR, Merritt DM, et al. *A neurotransmitter
  atlas of C. elegans males and hermaphrodites.* eLife (reviewed preprint), 2024.
- **DOI:** 10.7554/eLife.95402.2 (version 2) · https://elifesciences.org/reviewed-preprints/95402v2
- **Crawled (UTC):** 2026-07-03T00:00:22Z
- **Lab:** Hobert lab. Assignments are based on CRISPR-engineered reporter alleles for the
  transmitter-pathway genes (eat-4/VGLUT → glutamate, unc-17/cha-1 → ACh, unc-25/unc-47 →
  GABA, cat-2/bas-1 → dopamine, tph-1/cat-1 → serotonin, tdc-1/tbh-1 → tyramine/octopamine),
  cross-checked against scRNA-seq (CeNGEN) and prior antibody data.

## Files

| File | What |
|------|------|
| `TableS2_expression.xlsx` | pinned supplement (eLife `.../95402/v2/…/573258_file07.xlsx`) — per-neuron reporter table + consolidated "Neurotransmitter(s)" call |
| `parse_atlas.py` | reproducible parser (`uv run --with openpyxl python parse_atlas.py`) |
| `assignments.csv` | parsed **302** per-neuron assignments (`class, neuron, neurotransmitters`) |

`TableS2_expression.xlsx` sha256: `a2c211f763560c74fffdcb4881c8b704a242083d1e59f58a090d1b390fbe6f2d`

Other supplements (not vendored) from the same version, for reference:
`file06` (reporter-detail atlas) `4a39fe0b…`, `file08` `6cec5b91…`, `file09` `6d28a147…`,
`file10` `8f70df6c…` — at `https://prod--epp.elifesciences.org/api/files/95402/v2/content/supplements/573258_file<NN>.xlsx`.

## Parse notes

- The consolidated **"Neurotransmitter(s)"** column is used (the per-reporter columns encode
  expression by cell fill colour, which isn't recoverable from cell values).
- Values: `ACh`, `Glu`, `GABA`, `DA`, `5-HT`, `octopamine`, `betaine`, and `unknown (orphan)`;
  `*` marks assignments updated by this study, `- NEW` newly identified, `(uptake)` = uptake
  (not synthesis). These markers are preserved verbatim in `assignments.csv`.
- Covers both sexes; per-neuron rows are the hermaphrodite unless noted. For neurons with a
  sex- or stage-specific transmitter switch (e.g. AIM), the paper's text is the authority.

Attribution: cite Wang et al. 2024 (eLife) when using these assignments.
