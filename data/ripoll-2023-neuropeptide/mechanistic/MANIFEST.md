# Ripoll-Sánchez et al. 2023 — neuropeptide connectome *mechanistic layer*

The base [../MANIFEST.md](../MANIFEST.md) ingests the neuropeptide connectome as a network: a
directed edge source→target weighted by the **number** of NPP–GPCR pathways linking the pair. This
folder adds the mechanistic decomposition — **which** peptide→receptor pairs mediate each edge.

## Source
Same publication (Ripoll-Sánchez et al. 2023, *Neuron* 111:3570–3589; PMID 37935195) and repository
(github.com/LidiaRipollSanchez/Neuropeptide-Connectome). NPP–GPCR pairs are the biochemically
deorphanized interactions of Beets et al. 2023 (*Cell Rep* 42:113058; PMID 37656621;
doi:10.1016/j.celrep.2023.113058), filtered to **EC50 ≤ 500 nM**.

## Vendored files

### `npp_gpcr_pairs.csv` — the 92 validated NPP–GPCR pairs
Columns: `pair_index` (1–92, matching the paper's per-pair network order), `ligand` (NPP gene,
lowercase), `gpcr` (receptor gene, lowercase), `ec50_nm` (EC50 in nanomolar), `npp_family`,
`gpcr_class`. Built by joining the canonical pair order in
`Adjacency matrices for networks/neuropeptide_pairs (network identities for Individual_net folders).csv`
with EC50 / family / class from `Supplementary Tables/Table S5. NPP-GPCR pairs information.xlsx`.
All 92 pairs matched an EC50 row.

### `edge_pairs.csv` — per-edge pair attribution
Columns: `source`, `target`, `pair_index`. One row per (edge, mediating pair): edge source→target
is mediated by pair `pair_index` (i.e. source expresses that pair's ligand and target expresses its
GPCR). 145,834 rows. Cell names normalized to CIRCE (`DA01`→`DA1`).

Derived from the paper's 92 **individual** per-pair adjacency matrices
(`Adjacency matrices for networks/Individual NPP-GPCR networks LR/`), which decompose the aggregate
network: file *i* is pair *i*'s source×target binary matrix. This is the **long-range** superset;
the mid- and short-range models are edge subsets of it with **identical** per-edge pair sets (the
range models gate edge *existence* by anatomical reach, not the pathway composition of an edge), so
one attribution serves all three range databases.

## Validation
Summing the 92 individual matrices reproduces the published aggregate long-range network exactly:
**53,558 edges, and for every edge the pair count equals the published weight — 0 mismatches.**
So `edge_pairs.csv` is a faithful decomposition, not an estimate. (Cross-check: group by
`source,target` in `edge_pairs.csv` and the row count per edge equals that edge's `weight` in
`../long_range_network.csv`.)

### `np_genes.csv` — the 100 NPP + GPCR genes, resolved to WormBase
Columns: `symbol`, `wbgene` (WB:WBGene… curie), `systematic_name`, `category`
(`neuropeptide` for the 49 ligands, `neuropeptide_receptor` for the 51 GPCRs), `wbgene_source`.
WBGene resolution provenance (`wbgene_source`): 12 genes already in CIRCE (Cook 2020 / innexin
maps, reused for gene sharing); 84 by **expression-pattern alignment** — the gene's recovered
expression uniquely matches one WBGene row of the CeNGEN threshold-4 matrix
(`Scripts & data/30072020_CENGEN_threshold4_..._allneurons.csv`); this method was validated as
exact on all 12 pre-existing genes. The remaining 4 (`nlp-2`, `nlp-23` — collide on a sparse
2-neuron pattern; `dmsr-5`, `npr-34` — the two genes whose scRNA-seq the 2024 network corrected,
so they no longer match the older matrix) were resolved authoritatively via the Alliance / NCBI
gene APIs.

### `np_gene_expression.csv` — per-neuron NPP / GPCR expression
Columns: `cell`, `gene_symbol`. 3,514 (cell, gene) records over 302 neurons. **Recovered from the
92 individual per-pair matrices**: each matrix is the outer product (ligand expression) ×
(GPCR expression) binarized, so a matrix's row-support is exactly the neurons expressing its
ligand and its column-support the neurons expressing its GPCR. This recovers the *corrected*
expression the 2024 network was built on (so `dmsr-5`/`npr-34` are current), and every gene's
expression is identical across all pairs it appears in (verified).

## Validation
Beyond the per-edge check above, the expression + pairs together **derive** the network from first
principles — an edge source→target exists (weight = pair count) iff the source expresses a
ligand and the target its cognate GPCR for a validated pair. Joining `np_gene_expression.csv`
with `npp_gpcr_pairs.csv` this way reproduces the published long-range network **exactly**: 53,558
edges, identical to `edge_pairs.csv` — 0 differences. So the ingested expression is faithful and
the KG can regenerate the network without the pre-computed edge list.

## Notes / not included
- **Sensitivity variants** (other CeNGEN thresholds / EC50 cutoffs) are not ingested — only the
  reference threshold-4 / 500 nM expression and pairs.
- The expression covers the **NPP/GPCR genes of the 92 pairs**, the substrate that derives the
  network; the fuller CeNGEN matrix (NPR/MR/LGC genes without a ≤500 nM pair) is not ingested.
