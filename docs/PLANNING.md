# C. elegans Connectome Knowledge Graph — Planning

> Status: draft from planning conversation, 2026-06-23. This is the agreed scope before implementation.

## Goal

Build a re-runnable pipeline that converts **neuron-graph** (FunCoNN) connectivity data
into a **LinkML-modeled knowledge graph**, with every cell grounded in the
**C. elegans Gross Anatomy Ontology (WBBT)** and every connection traceable to its
source evidence.

The LinkML schema is the **single source of truth**; all outputs are generated from it.

## Consumers (what the graph is for)

1. **Triplestore + SPARQL** — canonical RDF/OWL store for querying.
2. **neuron-graph visualization** — enriched, ontology-linked data fed back into the
   existing viz; requires a JSON projection matching neuron-graph's shape.
3. **Direct research exploration** — a browsable artifact for the user / collaborators.

Implication: **RDF/OWL is the primary serialization** (consumers 1 & 3), with a
**neuron-graph-shaped JSON projection** generated alongside it (consumer 2). Because both
neuron-graph and WBBT evolve, this is a **pipeline**, not a one-shot conversion.

## Inputs (pinned, reproducible)

- **neuron-graph** — ingest from the **raw files in its repo** (the JSON loaded into MySQL
  + the `functional/randi_funconn_*.tsv` files), pinned as a versioned snapshot. No live
  API/DB dependency.
- **WBBT** — a **pinned release** of `wbbt.json` / `wbbt.owl` (via PURL or committed
  artifact), for reproducible name→URI mappings.
- **linkml-model / LinkML toolchain** — consumed as a normal Python dependency.

These three stay as external upstreams; this project does not vendor or fork them.

## Scope (v1)

- **Connection types:** all three — chemical synapse, gap junction, functional.
- **Datasets:** all datasets/specimens. The same neuron pair therefore yields **multiple
  connection records** (one per observing dataset, each with its own weight) — matching
  neuron-graph's own per-dataset structure, so it round-trips back to the viz.
- **Cells:** neurons **and** muscles (and any other non-neuronal cells neuron-graph carries).
- **Provenance:** **first-class** — every connection traces to its source data.

## Data model (LinkML)

Core classes (names indicative):

- **`Cell`** — a neuron or muscle. Carries its neuron-graph name, type/attributes, and a
  link to its WBBT anatomy term (`anatomy: WBbt:NNNNNNN`).
- **`Connection`** — a **reified node** (not a plain triple), because it is an n-ary
  relation: `pre`, `post`, `type` (chemical/gap/functional), `weight`, `dataset`, and a
  link to evidence.
- **`Dataset`** — a neuron-graph specimen/dataset, with its metadata.
- **`Evidence`** — provenance layer, **plain LinkML slots** (no PROV-O). Designed to be
  **additive**: today it holds "aggregated weight from dataset X"; later it can hold
  "individual synaptic contact" (section number, volume, image) **without schema churn**.
  Open-world: `Connection → has many Evidence`.

### Anatomy grounding (the name→URI bridge)

Mapping is **lexical**: normalize each neuron-graph cell name and match it against WBBT
`rdfs:label` + synonyms (`hasExactSynonym` etc.). Most resolve automatically. Results fall
into three buckets, emitted as a **first-class match report**:

- **Matched** — auto-linked to a `WBbt_*` URI.
- **Ambiguous** — multiple hits, or only a related/broad synonym (lower confidence).
- **Unmatched** — no hit.

The ambiguous + unmatched tail becomes a flat **work-list (CSV/YAML)** a human resolves by
hand. (No AI curation tool in scope — pipeline only.)

## Outputs

- **RDF/OWL** (Turtle) — primary, for triplestore/SPARQL/research.
- **neuron-graph JSON** — projection matching `/api/cells`, `/api/connections` shapes,
  losslessly preserving pre/post/type/weight/dataset for the viz. **Note:** this matches the
  API *output* shape, which is **not** how data is loaded into the viz — see the finding
  below.
- **Match report** — matched / ambiguous / unmatched, doubling as QA signal + curation
  work-list.

### Finding: feeding the viz needs the *raw populate* format, not the API projection

(Recorded 2026-06-29, after a live load of the KG data into neuron-graph against a remote
MariaDB.) The `/api/cells` and `/api/connections` projection is correct as an **API-output**
artifact, but loading it straight into neuron-graph's `connections` table breaks the viz:

- The client's default (class-level) view queries connections by **cell class** (e.g.
  `cells=AVA`), relying on **class-expansion rows** (`AVA→X`, `pre→postClass`, …) that
  neuron-graph's `populate-database` generates from cell-level edges. The API-shaped
  projection has no class rows, so any multi-cell class (AVA, DA, VB, …) shows as
  unconnected. Singleton classes (DVB, AVL) and explicitly-input individual cells happen to
  work, which masks the problem.
- Gap junctions: the API projection is already reverse-merged into one sorted orientation,
  whereas the DB layer expects raw per-direction edges and merges at populate/query time.
- The viz also needs each dataset's `datatypes` (`cs`/`gj`/`fc`); it is **absent from
  neuron-graph's own `datasets.json`** for the EM datasets and must be derived from the
  connection types actually present per dataset.

**Implication:** to feed the viz from the KG, emit the **raw populate format** (per-dataset
`{pre, post, typ, syn}` cell-level files + `neurons.json` + `datasets.json` *with*
`datatypes`) and run neuron-graph's `populate-database` (which owns class expansion + gap
dedup + the `synapses` table). The current API-shaped projection stays useful for shape
validation and for an API-level consumer — it is just the wrong layer to bulk-load. A
viz-loading helper (KG → raw populate files) is a candidate follow-up, out of v1 scope.

## Pipeline stages

1. **ingest** — read pinned neuron-graph files → normalized intermediate records.
2. **match** — lexical name→WBbt resolution; emit match report + work-list.
3. **build** — instantiate LinkML data (cells, connections, datasets, evidence).
4. **export** — serialize RDF/OWL + neuron-graph JSON projection.

Each stage is independently runnable and re-runnable against updated inputs.

## Stack

- **Python**, LinkML toolchain (schema → RDF/JSON/Python dataclasses).
- Pure pipeline — **no TypeScript, no UI, no interview feature**. (The interview feature
  from course-video-manager was used only as a *planning device* — this conversation.)

## Proposed repo layout

```
celegans-connectome-kg/
├── src/celegans_connectome_kg/
│   ├── schema/                # LinkML .yaml (single source of truth)
│   ├── ingest/                # neuron-graph file readers
│   ├── match/                 # anatomy name→WBbt matcher + report
│   ├── build/                 # LinkML data assembly
│   └── export/                # RDF + neuron-graph JSON emitters
├── data/
│   ├── neuron-graph/          # pinned source snapshot
│   └── wbbt/                  # pinned WBBT release
├── outputs/                   # generated RDF, JSON, match report
├── docs/                      # this plan, ADRs
├── tests/
├── pyproject.toml
└── README.md
```

## Phased roadmap

- **Phase 0 — scaffold:** ✅ repo structure, pyproject, LinkML installed, pin inputs.
- **Phase 1 — schema:** ✅ author the LinkML schema (Cell, Connection, Dataset, Evidence);
  generate RDF/JSON-Schema/dataclasses; validate with a tiny hand-made sample.
- **Phase 2 — ingest + match:** ✅ read neuron-graph files; build the WBBT label/synonym
  index; produce the match report and work-list.
- **Phase 3 — build + export:** ✅ assemble full LinkML data; emit RDF + neuron-graph JSON;
  round-trip-check the JSON against neuron-graph's expected shape.
- **Phase 4 — load + verify:** ✅ load RDF into a triplestore (Oxigraph); write sample SPARQL
  queries; sanity-check counts vs. source.

## Open questions / to confirm later

- Exact WBBT release/version to pin.
- Which synonym types count as a confident match (exact only, vs. related/broad). **Resolved:**
  only `label` + `hasExactSynonym` are confident; `related`/`broad`/`narrow` → ambiguous.
- ~~Triplestore choice for Phase 4~~ **Resolved: Oxigraph** (embedded, pip-installable, lets
  CI run real SPARQL).
- ~~GitHub project creation~~ **Resolved:** public repo under `raymond91125`.
- ~~Curating the 114-cell work-list and the 14 class-level/fragment connection endpoints~~
  **Resolved:** all 114 cells curated to WBbt terms (BWM → individual cell terms, etc.); the
  11 mappable endpoints minted as stub cells, 3 EM fragments left unmapped. Every cell is now
  WBBT-grounded.
- **Follow-up (out of v1 scope):** a viz-loading helper that emits the raw populate format
  (KG → per-dataset `{pre,post,typ,syn}` + `datasets.json` with `datatypes`) so the KG can
  feed neuron-graph's `populate-database` directly. See the Outputs finding above.
```
