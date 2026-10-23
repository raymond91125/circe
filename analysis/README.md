# Analyses

Reproducible analyses over the CIRCE knowledge graph. Each is self-contained and deterministic
(fixed seeds; results independent of `PYTHONHASHSEED`), and reads the built graph
(`outputs/connectome.ttl`, from `uv run cckg build && uv run cckg export`) or its JSON projection.
See the conventions in [`../CLAUDE.md`](../CLAUDE.md) before adding one.

| Analysis | What it answers |
|----------|-----------------|
| [`connectivity_paths/`](connectivity_paths/) | **Are two cells connected — directly or indirectly — and by what paths?** Reachability, shortest path, and all routes ≤ k hops, per connection type or over the combined `synaptic` graph, with strength filtering. SPARQL (`queries/*.rq`) + a `paths.py` tool. |
| [`feeding_pathways/`](feeding_pathways/) | How the extra-pharyngeal "feeding" neurons AIN / ASI / AVK could reach the pharynx — a whole-brain imaging correlation turned into testable molecular hypotheses, answered entirely in SPARQL. |
| [`neurotransmitter_reconciliation.md`](neurotransmitter_reconciliation.md) | Reconciles neuron neurotransmitter assignments across three sources (`reconcile_neurotransmitters.py`). |
| `pharynx_motif_analysis.py` → [`pharynx-motifs.html`](pharynx-motifs.html) | Network-motif (triad) analysis of the Cook 2020 pharyngeal connectome. |
| `pharynx_communities.py` → [`pharynx-communities.html`](pharynx-communities.html) | Louvain community detection on the Cook 2020 pharyngeal connectome. |

Most scripts run with `uv run python analysis/<name>`; those needing extra libraries say so in their
docstring (e.g. `uv run --with networkx python analysis/connectivity_paths/paths.py …`).
