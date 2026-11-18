# Analyses

Reproducible analyses over the CIRCE knowledge graph. Each is self-contained and deterministic
(fixed seeds; results independent of `PYTHONHASHSEED`), and reads the built graph
(`outputs/connectome.ttl`, from `uv run cckg build && uv run cckg export`) or its JSON projection.
See the conventions in [`../CLAUDE.md`](../CLAUDE.md) before adding one.

| Analysis | What it answers |
|----------|-----------------|
| [`connectivity_paths/`](connectivity_paths/) | **Are two cells connected — directly or indirectly — and by what paths?** Reachability, shortest path, and all routes ≤ k hops, per connection type or over the combined `synaptic` graph, with strength filtering. SPARQL (`queries/*.rq`) + a `paths.py` tool. |
| [`feeding_pathways/`](feeding_pathways/) | How the extra-pharyngeal "feeding" neurons AIN / ASI / AVK could reach the pharynx — a whole-brain imaging correlation turned into testable molecular hypotheses, answered entirely in SPARQL. |
| [`can_neuron/`](can_neuron/) | **Is the CAN pair a chemical-transmitting neuron?** Connectome (zero chemical outputs) + CeNGEN scRNA-seq (lacks `snt-1`/`unc-10`) + Hobert 2015 reporter (sole neuron missing the `rab-3` pan-neuronal reporter) converge: no conventional fast synaptic transmission. |
| [`morillo_male_neuropeptide/`](morillo_male_neuropeptide/) | **The first predicted *male* neuropeptide connectome** (Phase C): the Ripoll-Sánchez long-range method (validated to reproduce the published 53,558-edge herm network exactly) applied to the Morillo 2025 male expression layer — 11,586 edges, 2,080 onto/from male-specific neurons. Absolute male resource, not a sex comparison (documented confound). |
| [`haque_morillo_validation/`](haque_morillo_validation/) | **Does Haque et al. 2025 corroborate the Morillo male layer?** Independent single-cell platform: ligands validate strongly (AUC ≈ 0.89), receptors moderately (≈ 0.78 on well-sampled neurons). Haque's own data independently confirms net neuropeptide male-bias (86 vs 47 genes), while Morillo's per-gene sex *direction* matches Haque only at chance (49%) — so remodeling claims come from Haque, the absolute connectome from Morillo. |
| [`contactome_synapse_consistency/`](contactome_synapse_consistency/) | **Do CIRCE's synapses fall within the contactome (same volume)?** A synapse needs a membrane contact, so within one EM reconstruction synapses ⊆ contacts. Holds at 97.9–99.3%; the **dauer same-scope control (99.0%)** is the one clean estimate (the White↔Brittin figures are lower bounds — contactomes are nerve-ring-only). Misses are mostly detection-limit synapses; the bilateral PLN adult cluster is a within-N2U discrepancy (likely incomplete Brittin adult contact-scoring), not a CIRCE error. The reverse, all-contact, cross-scoring QC direction Brittin 2021 does not report. |
| [`neurotransmitter_reconciliation.md`](neurotransmitter_reconciliation.md) | Reconciles neuron neurotransmitter assignments across three sources (`reconcile_neurotransmitters.py`). |
| `pharynx_motif_analysis.py` → [`pharynx-motifs.html`](pharynx-motifs.html) | Network-motif (triad) analysis of the Cook 2020 pharyngeal connectome. |
| `pharynx_communities.py` → [`pharynx-communities.html`](pharynx-communities.html) | Louvain community detection on the Cook 2020 pharyngeal connectome. |

Most scripts run with `uv run python analysis/<name>`; those needing extra libraries say so in their
docstring (e.g. `uv run --with networkx python analysis/connectivity_paths/paths.py …`).
