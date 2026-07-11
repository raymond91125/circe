# CIRCE

**C**onnectome **I**ntegration & **R**easoning for **C**. **E**legans — the
knowledge-graph core.

[![CI](https://github.com/raymond91125/circe/actions/workflows/ci.yml/badge.svg)](https://github.com/raymond91125/circe/actions/workflows/ci.yml)

CIRCE integrates published _C. elegans_ connectome datasets (which vary by age,
sex, and reconstruction) into a knowledge graph for reasoning over chained
connections.

**This repository is the knowledge-graph core** — a re-runnable pipeline that
converts the source connectivity data into a **LinkML-modeled knowledge graph**,
with every cell grounded in the **_C. elegans_ Gross Anatomy Ontology (WBBT)** and
every connection traceable to its source evidence. The LinkML schema is the single
source of truth; RDF/OWL (for SPARQL + research) and a neuron-graph-shaped JSON
projection (consumed by the [CIRCE visualization](https://github.com/raymond91125/circe-viz))
are generated from it.

See [`docs/PLANNING.md`](docs/PLANNING.md) for the full design and roadmap.

## Status

All roadmap phases implemented (ingest → match → build → export → verify). Anatomy curation
of the work-list tail (114 cells) is ongoing.

## Layout

```
src/celegans_connectome_kg/
  schema/   LinkML .yaml (single source of truth)   [Phase 1]
  ingest/   neuron-graph file readers               [Phase 2]
  match/    anatomy name->WBbt matcher + report      [Phase 2]
  build/    LinkML data assembly                      [Phase 3]
  export/   RDF + neuron-graph JSON emitters          [Phase 3]
  verify/   Oxigraph load + sample SPARQL             [Phase 4]
data/
  neuron-graph/   pinned source snapshot (see MANIFEST)
  wbbt/           pinned WBBT release (see MANIFEST)
outputs/          generated artifacts (gitignored)
```

## Setup

```bash
uv sync --extra dev      # create .venv and install
uv run cckg --help       # pipeline CLI
```

## Pipeline

```bash
uv run cckg ingest       # read pinned neuron-graph files -> normalized records
uv run cckg match        # resolve cell names to WBbt; write match report + work-list
uv run cckg build        # assemble the LinkML graph -> outputs/connectome.json
uv run cckg export       # serialize RDF/OWL + neuron-graph JSON projection
uv run cckg verify       # load the RDF into Oxigraph; check counts; run sample SPARQL
```

Each stage is independently runnable and re-runnable against updated inputs. Sample SPARQL
queries live in `src/celegans_connectome_kg/verify/queries/`.

## Inputs

Source data is pinned under `data/` for reproducibility; provenance and checksums are
recorded in each subdirectory's `MANIFEST.md`. neuron-graph and WBBT remain external
upstreams — this repo does not fork them.

## Related repositories

- [CIRCE](https://github.com/raymond91125/circe) (knowledge-graph core + project home — this repo)
- [CIRCE visualization](https://github.com/raymond91125/circe-viz) (interactive browser)
- [NemaNode](https://github.com/zhenlab-ltri/NemaNode) (upstream connectome viz, via FunCoNN)

## Data sources & credits

CIRCE integrates published data: chemical-synapse and gap-junction wiring from John White and
the [Zhen](https://www.zhenlab.com), [Samuel](https://scholar.harvard.edu/aravisamuel), and
[Lichtman](https://lichtmanlab.fas.harvard.edu) labs
([NemaNode](https://nemanode.org/) / Witvliet et al., 2021; Cook et al., 2019), and functional
connectivity from the [Leifer Lab](http://leiferlab.princeton.edu)
([Randi et al., 2023](https://www.nature.com/articles/s41586-023-06683-4)) — all grounded in the
_C. elegans_ Gross Anatomy Ontology (WBBT).

## License & contact

This project is in beta and has not yet been peer reviewed. Code is released under the
[MIT License](LICENSE).

CIRCE is developed by [WormBase](https://wormbase.org/) and the
[Alliance of Genome Resources](https://www.alliancegenome.org/).

Contact: raymond+github@caltech.edu
