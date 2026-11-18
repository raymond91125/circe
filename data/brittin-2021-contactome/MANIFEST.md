# Brittin et al. 2021 — L4 + adult nerve-ring contactome (membrane adjacency)

## Source
Brittin CA, Cook SJ, Hall DH, Emmons SW, Cohen N. (2021).
**A multiscale brain map derived from whole-brain volumetric reconstructions.**
*Nature* 591:105. PMID 33627874 · doi:10.1038/s41586-021-03284-x

Volumetric EM reconstruction of the *C. elegans* **nerve ring** in two animals — the **L4 (JSH)** and
**adult (N2U)** series — the SAME two electron-micrograph series **White et al. 1986** scored synapses
in. Both animals are **hermaphrodites**. Scope is the **nerve ring** (the restricted bilateral-cell
dataset of the paper). The "contactome" term is coined in this paper.

Segmentation software: github.com/cabrittin/parsetrakem2 · analysis: github.com/cabrittin/elegansbrainmap

## Contactome (membrane adjacency) — `brittin_2021_l4_contactome`, `brittin_2021_adult_contactome`

Physical membrane contacts (appositions) between neurons — the substrate that bounds the synaptic
connectome (a synapse can only form where membranes touch). Ingested as the **`adjacency`** connection
type (undirected; stored in both orientations like gap junctions; weight = **contact area in µm²**,
EM-measured; sex hermaphrodite; life stage **L4** / **adult**).

- **Source (not vendored — large third-party file, ~8.9 MB):** Supplementary Information 1
  (`NIHMS2028296-supplement-Supp_Info_1.csv`): a header-less CSV, **one row per membrane-contact
  instance scored in one EM section**, columns
  `[dataset, cell1, cell2, cell1_index, cell2_index, EM_section_number, membrane_contact_in_pixels]`,
  `dataset` ∈ {`l4`, `adult`}. 285,214 rows total (156,941 L4 + 128,273 adult).
- **`build_contactome.py`** reads it and emits the two committed edge lists by **summing the pixel
  column over all EM sections and neurite-segment index pairs per unordered cell pair**, then scaling
  **pixels → µm² by ×4.5e-4** (same scale for both series, per the Supp Info):
  - `l4_contactome.csv` — **4,763 undirected contacts over 181 cells**.
  - `adult_contactome.csv` — **5,245 undirected contacts over 185 cells**.
  - Format `cell_a,cell_b,contact_area` (µm²), `cell_a < cell_b`. Brittin zero-pads single-digit
    ventral-cord indices (`VB01`→`VB1`, `VC01`→`VC1`, `VD01`→`VD1`); the converter unpads them so all
    185 labels ground to existing CIRCE cells (0 new cells minted).
- **Validation (the point of the layer):** because each contactome and its synapses come from the
  SAME EM series, White synapse pairs are almost entirely a subset of contact pairs
  (synapses ⊆ contacts) — the within-animal substrate check:
  - **L4 vs `white_1986_jsh`: 1,347 / 1,357 = 99.3%**
  - **adult vs `white_1986_n2u`: 1,424 / 1,454 = 97.9%**
  Pinned in `tests/test_sex_build.py::test_brittin_contactome_ingested`.
- **Reproduce:** `python -I data/brittin-2021-contactome/build_contactome.py /path/to/Supp_Info_1.csv`
  (or set `BRITTIN_CONTACTOME_CSV`).

## Not ingested (scope)
- **Supplementary Information 3** (`Mᵟ/Cᵟ/Gᵟ` reference graphs, also on GitHub as `reference_graph_*`):
  these are the **thresholded** contact graphs (only contacts reproducible in ≥4 datasets, restricted
  to 173 cells), which are far sparser (~900–1,046 contacts) and fail the synapse⊆contact check — so
  CIRCE ingests the **raw SI 1 scoring** instead, not these derived subsets.
- **Supplementary Information 2 / 4 / 5** — reproducibility rasters, cluster assignments, and the
  White-vs-Brittin contact cross-check for 3 cells (an independent validation; not needed for ingest).
- **Synapses** — already in CIRCE as `white_1986_jsh` / `white_1986_n2u`; this layer adds only the
  membrane-contact substrate for those same two animals.
