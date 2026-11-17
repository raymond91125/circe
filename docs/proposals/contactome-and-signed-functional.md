# Proposal: a contactome (adjacency) layer + explicit functional sign

Two enrichments to CIRCE's **connection model**, scoped together because both add edge semantics the
current graph lacks:

- **A — Contactome (new data):** ingest *physical membrane adjacency* (which neurons touch),
  from Brittin et al. 2021 (nerve ring) and Yim et al. 2024 (dauer), as a new `adjacency`
  connection type.
- **B — Explicit functional sign (model enhancement, no new data):** make the
  excitatory/inhibitory **sign** of the Randi 2023 functional connections a first-class field.
  CIRCE already carries the sign implicitly (in the functional weight); this exposes it.

Status: scoping/feasibility — not yet built. Follows the Morillo proposal's structure.

## A. Contactome

### Why (biological significance)

Physical adjacency is the **substrate that bounds the synaptic connectome**: a chemical synapse or
gap junction can only form where two membranes touch, so synapses are a *subset* of contacts. Brittin
et al. 2021 (*Nature*) showed the nerve ring's local contact neighbourhoods **constrain** which
synapses form, and that contact-profile clustering reveals a layered brain architecture invisible in
the synapse graph (Moyle et al. 2021). Functionally, adjacency is the physical basis for
**extrasynaptic/volume transmission** (neuropeptides, monoamines act on nearby receptor-expressing
cells — exactly CIRCE's predicted NP/monoamine layers) and for **ephaptic coupling**; it also helps
rationalize the function-vs-anatomy gap Randi 2023 attributes to extrasynaptic signalling. For CIRCE
specifically, a contactome is the **physical bridge between the observed synaptic layer and the
predicted extrasynaptic layers**, and an empirical anatomical gate for short-range NP models (the
gate CIRCE otherwise lacks, especially for the male).

### Prior art (checked)

- **Brittin et al. 2021**, *Nature* 591:105 — the "contactome": all neuronal membrane contacts in the
  nerve ring, two volumetric reconstructions (**adult + L4**, each left/right), hermaphrodite.
  Data: Zenodo [`4383277`](https://zenodo.org/records/4383277) + GitHub `cabrittin/elegansbrainmap`
  (SQL/GraphML, ~17 GB) + **ASCII adjacency in WormAtlas SI3** (small).
- **Moyle et al. 2021** — contact-profile clustering (analysis, not a new dataset).
- **Yim et al. 2024**, *Nat Commun* — dauer **contactome** (contact-*area* matrix, symmetric, ~363
  pairs), alongside the dauer synaptic connectome CIRCE already ingests.
- No **whole-animal** or **male** contactome exists. CIRCE has **no adjacency layer** today.

### Feasibility

- **Ingestion route (recommended):** the OpenWorm **`cect`** package (ConnectomeToolbox) already
  packages both Brittin and Yim contactomes as graph datasets — export a clean edge list
  (`node, node, weight`) from it, avoiding Brittin's 17 GB raw SQL. **Fallbacks:** Brittin ASCII
  (WormAtlas SI3); Yim supplementary matrices.
- **Name reconciliation:** Brittin nerve-ring neurons and Yim dauer cells are already CIRCE cells
  (dauer via the Yim synaptic layer), so grounding is mostly a lookup — a short reconciliation like
  the Morillo one, not a new cell inventory.
- **Symmetric data:** physical apposition has no pre/post. Store both directions, reusing CIRCE's
  existing **gap-junction convention** (gap junctions are already stored in both orientations).

### Schema / build design

- New `ConnectionType` value **`adjacency`** — "observed physical membrane contact (neurons whose
  processes are apposed); not a synapse." (Enum currently: chemical, gap_junction, functional,
  neuropeptidergic, monoaminergic.)
- New datasets, each with distinct provenance and **units kept separate** (do NOT pool):
  `brittin_2021_contactome_adult`, `brittin_2021_contactome_l4` (herm; weight = δ, the number of the
  EM reconstructions a contact appears in, 1–4, and/or contact length) and `yim_2024_contactome`
  (dauer; weight = summed **contact area**). Document the weight unit in each `Dataset.description`.
- `sex` = hermaphrodite (Brittin) / dauer is a life stage (Yim, hermaphrodite); `life_stage` =
  adult / L4 / dauer.

### Caveats

- **Nerve ring only** (Brittin) and **dauer** (Yim) — not whole-animal, not male. Weight units are
  **not comparable across sources** (δ/length vs contact area). EM-adjacency has resolution limits.

## B. Explicit functional sign

> **Status: IMPLEMENTED** (option (a) — explicit `sign` + non-negative `weight`). Schema adds a
> `ConnectionSign` enum (`excitatory`/`inhibitory`) and a `sign` slot on `Connection`; the build
> splits the signed Randi weight into magnitude + sign (1,239 excitatory / 244 inhibitory); the viz
> projections reconstruct the signed value, so `kg_connections.json` is byte-identical (viz
> unaffected). Remaining parts A (contactome) and a possible part C (predicted chemical-synapse sign,
> Fenyves 2020 / derived from NT + ionotropic receptors) are still proposals.

### Finding (feasibility: the data is already here)

Randi 2023 functional connectivity is ingested from the neuron-graph `randi_funconn_*` files, whose
`syn` values are **signed** (`[22]` excitatory, `[-12]` inhibitory). CIRCE's functional weight =
`round(sum(syn))`, so the sign already flows through: the built graph has **1,239 excitatory (+) and
244 inhibitory (−)** functional edges. It is just **not modelled explicitly** — a consumer must
notice the weight can be negative.

### Proposal

- Add an explicit **`sign`** slot to `Connection` (enum `excitatory` / `inhibitory`; optional — only
  the functional layer populates it), derived from the functional weight's sign.
- Decide how to carry magnitude (open decision below). Enables clean queries ("inhibitory functional
  partners of AVA"), and reconciliation against **neurotransmitter identity** (e.g. a GABAergic
  source predicting an inhibitory edge) and against sign-prediction work (Fenyves 2020).

## Open decisions (for review)

1. **Adjacency directedness** — store both directions (gap-junction convention). *(recommended)*
2. **Adjacency weight** — keep raw per-source (δ/length for Brittin, contact area for Yim),
   non-pooled; expose normalization as a separate concern. *(recommended)*
3. **Functional sign representation** — (a) add `sign` enum **and** make `weight` the non-negative
   **magnitude** (cleanest: weight ≥ 0 for every type, sign separate); or (b) add `sign` enum but
   leave `weight` signed for backward-compatibility. *(a) is cleaner but touches existing functional
   weights; needs a call.*
4. **Scope** — Brittin + Yim now; no whole-animal/male contactome exists to add.

## Validation plan (per CLAUDE.md analysis conventions)

- **Contactome — synapses ⊆ contacts.** Within the reconstructed region, every synaptic (chemical /
  gap-junction) pair should also be a contact pair. Measuring that overlap is a strong correctness
  check (adjacency is the substrate of synapses) and quantifies the "potential synapse" set
  (contacts without synapses). Cross-check contact counts against Brittin's reported totals.
- **Functional sign — vs neurotransmitter identity.** Where the source neuron's transmitter is known
  (Wang 2024), check that inhibitory transmitters (GABA) trend to `inhibitory` functional edges;
  report concordance honestly (functional sign is net circuit effect, not a monosynaptic guarantee).

## References

- Brittin CA, et al. *Volumetric reconstruction of main Caenorhabditis elegans neuropil at two
  different time points.* Nature 591:105–110 (2021).
- Moyle MW, et al. *Structural and developmental principles of neuropil assembly in C. elegans.*
  Nature 591:99–104 (2021).
- Yim H, et al. *Comparative connectomics of dauer reveals developmental plasticity.* Nat Commun
  (2024).
- Randi F, et al. *Neural signal propagation atlas of C. elegans.* Nature 623:406–414 (2023).
