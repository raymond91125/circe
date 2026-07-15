# CIRCE — Connectome Integration & Reasoning for *C. elegans*

A LinkML-modeled, WBBT-grounded knowledge-graph pipeline that converts neuron-graph
connectivity data into RDF/OWL plus neuron-graph JSON projections (consumed by the
`circe-vis` viz). The LinkML schema is the source of truth:
`src/celegans_connectome_kg/schema/connectome.yaml`.

## Dev

```
uv run cckg build      # assemble the connectome from pinned sources -> outputs/connectome.json
uv run cckg export     # RDF/OWL + neuron-graph JSON projections -> outputs/
uv run cckg verify     # SPARQL cross-checks over the built graph
```

Gate before committing:

```
uv run --extra dev python -m pytest
uv run ruff check src/ tests/
uv run ruff format --check src/ tests/
```

## Analysis conventions

Before presenting any quantitative or statistical analysis (e.g. work under `analysis/`):

1. **Check prior art first.** If the dataset has a source publication, or there are standard
   references for the method, see whether they ran a comparable analysis and reconcile any
   differences *before* drawing conclusions.
2. **Justify method choices.** For each pivotal choice (null model, threshold, normalization,
   statistical test), state what was chosen, why, the standard alternative, and whether the
   conclusion survives that alternative. For enrichment / motif tests specifically, the null
   must control for all lower-order structure (e.g. preserve reciprocity when testing 3-node
   motifs — otherwise a dyadic property masquerades as a triadic one).
3. **Sensitivity-check the headline.** Re-run the main result under at least one reasonable
   alternative assumption; report if it flips.
4. **Treat surprising results as red flags.** A finding that contradicts established results is
   probably an artifact until shown otherwise — say so, and dig in before presenting it.
5. **Calibrate confidence.** Separate robust findings from ones that hinge on a method choice.

Make analyses reproducible and deterministic (fixed seeds; results independent of
`PYTHONHASHSEED`), and commit them under `analysis/`.
