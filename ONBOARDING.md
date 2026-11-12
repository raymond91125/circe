# CIRCE — contributor onboarding

**CIRCE** (Connectome Integration & Reasoning for *C. elegans*) is a LinkML-modeled,
WBBT-grounded knowledge-graph pipeline. It integrates published *C. elegans* connectome
datasets (varying by age, sex, reconstruction) into RDF/OWL plus a neuron-graph JSON
projection consumed by the [circe-viz](https://github.com/raymond91125/circe-viz) viewer.

**The LinkML schema is the single source of truth:**
`src/celegans_connectome_kg/schema/connectome.yaml`. RDF/OWL and the JSON projection are
*generated* from it — never hand-edit the outputs.

## Setup & pipeline

Uses [`uv`](https://docs.astral.sh/uv/). The pipeline is three stages:

```
uv run cckg build      # assemble from pinned sources -> outputs/connectome.json
uv run cckg export     # RDF/OWL + neuron-graph JSON projections -> outputs/
uv run cckg verify     # SPARQL cross-checks over the built graph
```

## The gate (run before every commit)

```
uv run --extra dev python -m pytest
uv run ruff check .            # CI lints the WHOLE repo, not just src/tests
uv run ruff format --check src/ tests/
```

CI runs the test suite on Python 3.10 / 3.11 / 3.12; all three are required to merge.

## Repo layout

| Path | What's there |
|------|--------------|
| `src/celegans_connectome_kg/schema/connectome.yaml` | **Source of truth** — the LinkML schema |
| `src/celegans_connectome_kg/` | build / export / verify pipeline, ingest adapters |
| `data/` | pinned source datasets, each with a `MANIFEST.md` documenting provenance |
| `outputs/` | generated graph (`connectome.json`, `connectome.ttl`) + JSON projections |
| `docs/` | generated schema docs; `docs/PLANNING.md` is the design + roadmap |
| `tests/` | the test suite |
| `analysis/` | **worked analyses over the built graph** (see below) |

## Analyses

Worked, reproducible analyses live in [`analysis/`](analysis/) — each self-contained and
deterministic (fixed seeds; independent of `PYTHONHASHSEED`), reading the built graph. Start
from [`analysis/README.md`](analysis/README.md). Before adding one, read the analysis
conventions in [`CLAUDE.md`](CLAUDE.md) (check prior art, justify method choices,
sensitivity-check the headline, calibrate confidence).

Current analyses include:

- [`analysis/can_neuron/`](analysis/can_neuron/) — **Is the CAN pair a chemical-transmitting
  neuron?** Combines the CIRCE connectome (zero chemical outputs), CeNGEN scRNA-seq (CAN lacks
  `snt-1`/`unc-10`), and Hobert 2015 reporter data (CAN is the sole neuron missing the
  pan-neuronal `rab-3` reporter). Write-up in
  [`article.md`](analysis/can_neuron/article.md); regenerate CIRCE-side facts with
  `uv run python analysis/can_neuron/can_circe.py`.
- [`analysis/connectivity_paths/`](analysis/connectivity_paths/) — whether/how two cells
  connect (reachability, shortest path, routes ≤ k hops), SPARQL + a `paths.py` tool.
- [`analysis/feeding_pathways/`](analysis/feeding_pathways/) — how extra-pharyngeal feeding
  neurons could reach the pharynx, answered in SPARQL.
- [`analysis/neurotransmitter_reconciliation.md`](analysis/neurotransmitter_reconciliation.md)
  — reconciles neurotransmitter assignments across three sources.
- `analysis/pharynx_motif_analysis.py` / `pharynx_communities.py` — triad-motif and community
  analysis of the Cook 2020 pharyngeal connectome.

## Contribution workflow

Branch off `main`, run the gate, open a PR (CI must be green), squash-merge, delete the branch.
See recent PRs (#76–#80) for the house style of commits and PR descriptions.
