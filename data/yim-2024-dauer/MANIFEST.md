# Yim, Choe, Bae et al. 2024 — dauer nerve-ring connectome

## Source
Yim H, Choe DT, Bae JA, Choi M-K, Kang H-M, Nguyen KCQ, et al. (2024).
**Comparative connectomics of dauer reveals developmental plasticity.**
*Nat Commun* 15:1546. PMID 38413604 · doi:10.1038/s41467-024-45943-3

A complete electron-microscopy reconstruction of the *C. elegans* **dauer** nerve ring, with
deep-learning synapse detection. Dauer is the stress-induced alternative third larval stage; the
animal is a **hermaphrodite**. Scope is the **nerve ring only** (not the whole nervous system).

EM volume & segmentation: BossDB (`yim_choe_bae2023`). Analysis code: github.com/jabae/Yim-Choe-Bae-et_al-2024.

## Vendored data
`dauer_connections.csv` — the aggregated **chemical** connectome (`pre,post,synapses,size_sum`),
derived deterministically from Supplementary Data (`MOESM5`, sheet `dauer_221021`), the per-synapse
table (`syn_id, pre_neuron, post_neuron, size`).

- **6,371 synapses → 2,200 directed edges** over 221 partners (181 neurons + muscle / other cells).
- **`synapses`** — number of synapses per ordered pair. Used as the connection **weight**, matching
  the synapse-count convention of the White/Witvliet datasets already in the KG.
- **`size_sum`** — summed active-zone size (voxels) per edge, the paper's own primary weight metric
  (active-zone size as a proxy for synaptic strength). Retained for reference; not currently loaded.

## Build
Ingested as the hermaphrodite, **dauer**-stage dataset **`yim_2024_dauer`** over cells that already
exist in the KG. All edges are **chemical**: the study did not reconstruct gap junctions ("Gap
junctions … are not included in this study").

Life stage is recorded via the `Dataset.life_stage` slot (see `data/curation/dataset_life_stage.csv`).
The one partner not already in neuron-graph's neuron list, **`exc_duct`** (the excretory duct cell,
postsynaptic in 2 edges), is added as a curated endpoint cell grounded to WBbt:0004540.

## Not ingested (scope)
- **Gap junctions** — not reconstructed by the study.
- The paper's re-tabulated **L1–L4 / adult** comparison connectomes (`MOESM6`/`MOESM9` adjacency
  matrices): those stages are already represented in the KG by Witvliet 2020 and White 1986.
- The **active-zone size** weighting (kept in `size_sum` for reference, but weight = synapse count).
