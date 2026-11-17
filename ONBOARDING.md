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

## Deployment

CIRCE runs on two hosts: a dev box and the **public production** host `wobrbig.caltech.edu`
(viz at **:3000**, SPARQL endpoint at **:3001** — **both HTTPS**).

### Public SPARQL endpoint (`:3001`), after a KG change

`circe-sparql` is a **system** service that serves `outputs/connectome.ttl`, loaded into memory at
startup. After the KG changes (data, schema), regenerate the TTL and reload it:

```
ssh <you>@wobrbig.caltech.edu
cd ~/local/src/git/circe
git pull --ff-only
PATH="$HOME/.local/bin:$PATH" uv run cckg build && uv run cckg export   # uv isn't on the non-interactive PATH
sudo systemctl restart circe-sparql     # system service -> needs sudo; it reloads the ~150 MB TTL (~1 min)
```

Query it **over HTTPS** — the cert is self-signed, so pass `-k`; a plain-HTTP request fails with an
SSL error and returns nothing:

```
curl -sk https://wobrbig.caltech.edu:3001/query \
  -H 'Content-Type: application/sparql-query' \
  --data-binary 'PREFIX cckg: <https://wormbase.org/resources/connectome/>
SELECT (COUNT(?c) AS ?n) WHERE { ?c a cckg:Connection }'
# GET https://…:3001/query?query=<urlencoded> also works; https://…:3001/queries lists sample queries
```

### Viz (`circe-viz`, `:3000`)

The viz reads its own MySQL DB, populated from circe's `outputs/`. After a KG change, in the
`circe-viz` repo run `npm run refresh-kg` (regenerates the bundled client maps + server raw-data from
circe), commit, then on the host: `npm run build-prod` (client bundle, served from `dist/` — live on
rebuild, no restart) and `npm run populate-database` (truncates + reloads the serving DB; brief blip).
`:3000` redirects HTTP→HTTPS. Browsers may cache the old bundle — hard-refresh to pick up a new one.
