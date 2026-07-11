# C. elegans Connectome Knowledge Graph

[![CI](https://github.com/raymond91125/circe/actions/workflows/ci.yml/badge.svg)](https://github.com/raymond91125/circe/actions/workflows/ci.yml)

A re-runnable pipeline that converts **neuron-graph** (FunCoNN) connectivity data into a
**LinkML-modeled knowledge graph**, with every cell grounded in the
**C. elegans Gross Anatomy Ontology (WBBT)** and every connection traceable to its source
evidence.

The LinkML schema is the single source of truth; RDF/OWL (for SPARQL + research) and a
neuron-graph-shaped JSON projection (for the viz) are generated from it.

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
upstreams — this repo does not vendor or fork them.
