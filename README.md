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

## Schema

The LinkML schema is the single source of truth
([`src/celegans_connectome_kg/schema/connectome.yaml`](src/celegans_connectome_kg/schema/connectome.yaml)).
Browsable docs and machine-readable exports are generated from it:

- **Docs site:** <https://raymond91125.github.io/circe/> (classes, slots, enums; published to
  GitHub Pages by [`.github/workflows/schema-docs.yml`](.github/workflows/schema-docs.yml)).
- **OWL:** [`docs/schema/connectome.owl.ttl`](docs/schema/connectome.owl.ttl) — the schema/TBox
  for triplestores and reasoners (pairs with the RDF data from `cckg export`).
- **JSON Schema:** [`docs/schema/connectome.schema.json`](docs/schema/connectome.schema.json) —
  for validating data.

Regenerate all three with `sh scripts/gen-schema-docs.sh` (preview the site:
`uv run --extra docs mkdocs serve`).

## Query it (SPARQL)

The full graph is available as RDF/Turtle — download `connectome.ttl` from the
[latest release](https://github.com/raymond91125/circe/releases/latest) and load it into any
triplestore or `rdflib`/[Oxigraph](https://github.com/oxigraph/oxigraph).

To explore interactively, run a local read-only SPARQL endpoint over your build:

```bash
uv run cckg build && uv run cckg export   # produces outputs/connectome.ttl
sh scripts/serve-sparql.sh                # -> http://127.0.0.1:7878/  (query console)
```

It serves a browser query console plus a `/query` HTTP endpoint (SPARQL-JSON / CSV / TSV, and
Turtle for CONSTRUCT), auto-prepends the common prefixes (`cckg:`, `WBbt:`, …), and refuses
updates. The shipped example queries in
[`src/celegans_connectome_kg/verify/queries/`](src/celegans_connectome_kg/verify/queries/) are
loadable from the console. Pass `--host 0.0.0.0 --port <p>` to expose it beyond localhost.

## Related repositories

- [CIRCE](https://github.com/raymond91125/circe) (knowledge-graph core + project home — this repo)
- [CIRCE visualization](https://github.com/raymond91125/circe-viz) (interactive browser)
- [NemaNode](https://github.com/zhenlab-ltri/NemaNode) (upstream connectome viz, via FunCoNN)

## Data sources & credits

CIRCE integrates published data across several modalities. Each source is pinned under `data/`
with provenance and checksums in its `MANIFEST.md`.

**Connectomes — chemical synapses & gap junctions**

- **White et al., 1986** — the original serial-section EM reconstruction of the adult
  hermaphrodite ([10.1098/rstb.1986.0056](https://doi.org/10.1098/rstb.1986.0056)).
- **Witvliet et al., 2021** — connectomes across development (Zhen / Samuel / Lichtman labs)
  ([10.1038/s41586-021-03778-8](https://doi.org/10.1038/s41586-021-03778-8)).
- **Cook et al., 2019** — whole-animal connectomes of both sexes (Emmons Lab / WormWiring)
  ([10.1038/s41586-019-1352-7](https://doi.org/10.1038/s41586-019-1352-7)).
- **Cook et al., 2020** — the pharyngeal connectome (Emmons Lab)
  ([10.1002/cne.24932](https://doi.org/10.1002/cne.24932)).

**Functional connectivity**

- **Randi et al., 2023** — a neural signal-propagation atlas from optogenetic activation +
  calcium imaging (Leifer Lab)
  ([10.1038/s41586-023-06683-4](https://doi.org/10.1038/s41586-023-06683-4)).

**Neurotransmitter identity**

- **Wang et al., 2024** — a *C. elegans* neurotransmitter atlas
  ([10.7554/eLife.95402](https://doi.org/10.7554/eLife.95402)).
- **Hardege et al., 2022** — betaine neurotransmission at RIM
  ([10.1073/pnas.2201783119](https://doi.org/10.1073/pnas.2201783119)).
- [WormAtlas Neurotransmitters Table](https://www.wormatlas.org/neurotransmitterstable.htm).

**Anatomy ontology**

- Every cell is grounded in the *C. elegans* Gross Anatomy Ontology (**WBBT**), maintained by
  [WormBase](https://wormbase.org/).

The connectome and functional data reach CIRCE via a pinned snapshot of
[NemaNode](https://nemanode.org/) / [FunCoNN](https://funconn.princeton.edu/) (White, Witvliet,
Randi); Cook, the neurotransmitter sources, and WBBT are pinned directly.

## License & contact

This project is in beta and has not yet been peer reviewed. Code is released under the
[MIT License](LICENSE).

CIRCE is developed by [WormBase](https://wormbase.org/) and the
[Alliance of Genome Resources](https://www.alliancegenome.org/).

Contact: raymond+github@caltech.edu
