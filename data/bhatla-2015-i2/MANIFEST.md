# Bhatla et al. 2015 — pharyngeal I2 neuron EM synapses

## Source
Bhatla N, Droste R, Sando SR, Huang A, Horvitz HR (2015).
**Distinct neural circuits control rhythm inhibition and spitting by the myogenic pharynx of
*C. elegans*.** *Curr Biol* 25(16):2075–2089. PMID 26212880 · doi:10.1016/j.cub.2015.06.052

Data: Supplemental Data 1 (`mmc1.pdf`), **Figure S6 panel E** — two tables ("I2L chemical synapses",
"I2R chemical synapses") quantifying, per postsynaptic partner, the chemical synapses of the
pharyngeal I2 neurons reconstructed by EM (anterior neurites from worm #1, posterior from worm #4).

## Vendored data
`i2_synapses.csv` — the transcribed edge list (`pre,post,sections,synapses`). The source table is
**image-only in the PDF**, so it was transcribed by hand from a high-resolution render.

- **`pre`** is `I2L` / `I2R`; **`post`** is the postsynaptic partner, normalized to CIRCE cell names
  (SI6 `PM*` → `pm*`, `BM` → `bm`; `e3VL`/`e3VR` and neuron names kept).
- **`sections`** — the number of EM sections over which the synapses onto the partner are
  distributed (Fig S6E column 4). Used as the connection **weight**, comparable to the Cook
  EM-serial-section metric.
- **`synapses`** — the synapse count "in this work" (Fig S6E column 3), retained for reference.

## Build
Ingested as the hermaphrodite dataset **`bhatla_2015_i2`** (26 chemical connections) over cells that
already exist in the KG. Notably includes **I2 → pharyngeal-muscle synapses** (pm1, pm3VL/VR,
pm4/pm4VR, pm5, bm) that are **absent from the White/Cook connectomes**.

## Not ingested (scope)
- The **Albertson & Thomson 1976** comparison column and the **DP-volume / vesicle-area**
  morphometrics (Fig S6E columns 2, 5, 6) — this pass is connectivity only.
- Rows with **0 synapses in this work** (reported by A&T 1976 but not confirmed here: `MCL` for
  both cells, the posterior `M1` row for I2L) are omitted.
