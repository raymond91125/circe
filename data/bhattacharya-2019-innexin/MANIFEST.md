# Bhattacharya et al. 2019 — innexin (gap-junction) gene expression

## Source
Bhattacharya A, Aghayeva U, Berghoff EG, Hobert O (2019). **Plasticity of the Electrical Connectome
of *C. elegans*.** *Cell* 176(5):1174–1189.e16. PMID 30686580 · doi:10.1016/j.cell.2018.12.024

**Figure 1B** — a color-coded matrix of neuronal innexin gene expression: 17 neuronally-expressed
innexins × 121 neuron classes, scored as expressed in **both** non-dauer and dauer (yellow),
**non-dauer only** (red), **dauer only** (green), or not expressed (white); `*` marks developmental
plasticity within non-dauer.

## Vendored data
`innexin_expression.csv` (`neuron, innexin, expression, developmental_plasticity`) and
`innexin_genes.csv` (`fig_label, symbol, isoform, wbgene, category, systematic_name`).

**Figure 1B has no machine-readable supplement** — it is published only as an image. This table was
**extracted programmatically** from the high-resolution figure (grid detection + per-cell color
classification), then **manually verified**. The extraction was cross-checked three ways:
- **Panel 1C** (the paper's own per-innexin neuron-class counts) — matched within ~1–2 for all 17.
- **Table S2** (developmental changes) — 23/24 listed pairs recovered as `*`.
- **Table S1** (non-neuronal expression) — the 10 innexins S1 reports as gut/gonad/muscle-only are
  correctly absent from the neuronal Fig 1B; the 17 present are the neuronal subset.

Genes are keyed to persistent WormBase gene IDs (Alliance-resolved); isoforms (inx-1 a/b, inx-10 a,
inx-18 a/b) are kept as qualifiers on the same gene, matching the Cook 2020 SI6 convention.

## Build
Ingested as two hermaphrodite gene-expression datasets so the dauer plasticity is explicit:
`bhattacharya_2019_innexin` (non-dauer) and `bhattacharya_2019_innexin_dauer` (life_stage = dauer).
An expression scored "both" yields a record in **both** datasets; "non-dauer only" only in the
former; "dauer only" only in the latter. Category `innexin`, over cells already in the KG.

## Not ingested / notes
- The `*` developmental-plasticity flag is retained in the CSV; recorded on the non-dauer record.
- Non-neuronal innexin expression (Table S1) and developmental detail (Table S2) are not ingested.
