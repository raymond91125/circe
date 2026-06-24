# Pinned input — neuron-graph

Source: `mattpocock`-style sibling repo at `/home/raymond/local/src/git/neuron-graph`
(FunCoNN, a NemaNode fork).

- Upstream commit: `6ea5f9b55715e68eee4089dec2b3a4c2574d23e2`
- Pinned: 2026-06-23

## Files pinned

| File | sha256 | notes |
|------|--------|-------|
| `randi_funconn_unc31.tsv`  | `21a626af9dae7795ff16761196312fa32cb060b834819ded8decaf877f27431d` | functional connectivity, unc-31 mutant |
| `randi_funconn_wildcp.tsv` | `cb89227b9e9b5b6f8f277dae556b544b0fbc0108efe287c2807fb290ffc8a174` | functional connectivity, wild-type (see note) |
| `randi_funconn_wildty.tsv` | `cb89227b9e9b5b6f8f277dae556b544b0fbc0108efe287c2807fb290ffc8a174` | functional connectivity, wild-type (see note) |

**Note:** `wildcp` and `wildty` have identical checksums in the upstream repo — confirm in
Phase 2 whether this is intentional (same data under two dataset labels) or an upstream
duplication before treating them as distinct datasets.

## Still to pin (Phase 2)

The **chemical synapse** and **gap junction** (EM) connectivity and the **cell** list are
loaded into neuron-graph's MySQL via its `populate-database` step; the source JSON for that
was not located during Phase 0 (only test fixtures under `test/` were found). Phase 2
(ingest) must locate the real populate-database inputs in neuron-graph and pin them here.
