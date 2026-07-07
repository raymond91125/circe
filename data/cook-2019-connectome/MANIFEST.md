# Pinned reference — Cook et al. 2019 whole-animal connectomes (both sexes)

The adult **hermaphrodite and male** whole-animal connectomes, as adjacency matrices. This is the
canonical source for the **male** connectome and is the input for extending the KG to a
**sex-aware** graph (male + hermaphrodite). See `docs/PLANNING.md` for the sex-extension plan;
this directory is phase **M1** (pin data + provenance) only — no parser/schema changes yet.

## Source

- **Paper:** Cook SJ, Jarrell TA, Brittin CA, Wang Y, Bloniarz AE, Yakovlev MA, Nguyen KCQ,
  Tang LT-H, Bayer EA, Duerr JS, Bülow HE, Hobert O, Hall DH, Emmons SW. *Whole-animal
  connectomes of both Caenorhabditis elegans sexes.* Nature 571:63–71 (2019).
- **DOI:** 10.1038/s41586-019-1352-7 · https://www.nature.com/articles/s41586-019-1352-7
  (open full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC6889226/)
- **Data host:** WormWiring (Emmons Lab), https://wormwiring.org/pages/adjacency.html
- **File URL:** `https://wormwiring.org/si/SI 5 Connectome adjacency matrices, corrected July 2020.xlsx`
- **Crawled (UTC):** 2026-07-07T22:19:03Z
- **Version:** the **"corrected July 2020"** revision (supersedes the original `SI 5 Connectome
  adjacency matrices.xlsx`; the correction reconciles gap-junction diagonal/symmetry
  inconsistencies against the asymmetric tables).

## Files

| File | What |
|------|------|
| `SI5_connectome_adjacency_matrices_corrected_2020.xlsx` | pinned supplement, verbatim (original name `SI 5 Connectome adjacency matrices, corrected July 2020.xlsx`) |

`SI5_connectome_adjacency_matrices_corrected_2020.xlsx` sha256:
`1f4fdbf84746b69b49a8da0816f52787860ce349b638dce37924ba80f90c70c9`

## Workbook structure

Seven sheets — a legend plus one adjacency matrix per (sex × connection-type):

| Sheet | Dimensions (rows × cols) | Contents |
|------|------|------|
| `TITLE AND LEGEND` | — | provenance + weight semantics (see below) |
| `hermaphrodite chemical` | 304 × 458 | directed chemical adjacency (row = pre, col = post) |
| `hermaphrodite gap jn asymmetric` | 475 × 476 | gap junctions, per-direction |
| `hermaphrodite gap jn symmetric` | 475 × 476 | gap junctions, symmetrized |
| `male chemical` | 386 × 579 | directed chemical adjacency (row = pre, col = post) |
| `male gap jn asymmetric` | 590 × 590 | gap junctions, per-direction |
| `male gap jn symmetric` | 590 × 590 | gap junctions, symmetrized |

Each matrix has a region-grouping header (e.g. `PHARYNX`) above the cell-name header row/column;
cell names begin ~row 3 / col 3 and weights fill the body. Nodes include **neurons, muscles, and
non-muscle end organs** (gonad, gut, hypodermis, etc.). The male graph is ~579 nodes (385 neurons,
155 muscles, 39 end organs); the hermaphrodite is smaller (no male-specific neurons).

## Weight semantics (from the legend, verbatim intent)

> "For chemical connections, rows give the pre-synaptic cell, columns give the post-synaptic cell.
> **Weights in the body of the matrices are the total number of EM serial sections of connectivity,
> taking into account both the number of synapses and the sizes of synapses.** To provide complete
> coverage of the entire nervous system, the data are assembled from multiple animals and include
> connections added by extrapolation in gaps where no data were available."

**Important:** these weights are **serial-section counts (synapse number × size)**, *not* raw
synapse counts like the neuron-graph datasets (White 1986, Witvliet). Downstream ingest (phase M3)
must record this so weights are not silently compared across sources. Gap junctions: prefer the
**symmetric** sheet (matches the KG's existing undirected/merged gap-junction handling).

## Licensing / attribution

The WormWiring SI files carry "Emmons Lab Copyright (c) 2020" with **no explicit reuse license**.
Standard practice is to **cite Cook et al. 2019 (Nature)** when using these data. The same data are
also redistributed by the MIT-licensed OpenWorm ConnectomeToolbox
(https://github.com/openworm/ConnectomeToolbox, `Cook2019MaleReader`) if an explicitly-licensed
copy is later preferred. Cite Cook et al. 2019 in any KG release derived from these matrices.
