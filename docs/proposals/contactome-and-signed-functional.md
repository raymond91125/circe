# Proposal: a contactome (adjacency) layer + explicit functional sign

Enrichments to CIRCE's **connection model**, scoped together because all add edge semantics the
current graph lacks:

- **A — Contactome (new data):** ingest *physical membrane adjacency* (which neurons touch) as a new
  `adjacency` connection type. *(IMPLEMENTED — three contactomes, each scored on the same EM volume as
  a synaptic connectome already in CIRCE, so synapses ⊆ contacts validates within-animal: Yim 2024
  **dauer** (~99%), Brittin 2021 **L4** (99.3%) and **adult** (97.9%). Brittin came from the raw
  per-section scoring in **Supplementary Information 1**, not the 17 GB Zenodo MySQL and not the
  repo's committed graphs — those are a thresholded analysis subset that fails the subset check.)*
- **B — Explicit functional sign (model enhancement, no new data):** make the
  excitatory/inhibitory **sign** of the Randi 2023 functional connections a first-class field.
  *(IMPLEMENTED — see below.)*
- **C — Predicted chemical-synapse sign:** assign an excitatory/inhibitory **sign** to *chemical*
  synapses (far broader than Randi's measured functional subset), from presynaptic neurotransmitter
  + postsynaptic ionotropic receptor — reusing B's `sign` slot. *(data + build IMPLEMENTED,
  increments 1–3; validation analysis (increment 4) pending — see below)*

Status: A and B are built; C's data and derivation build are built (increments 1–3), with the
validation analysis (increment 4) pending. Follows the Morillo proposal's structure.

## A. Contactome

> **Status: IMPLEMENTED** — three hermaphrodite contactomes as `adjacency` edges, each validated
> within-animal against the matching synaptic dataset (synapses ⊆ contacts):
> - `yim_2024_dauer_contactome` — 6,080 undirected → 12,160 edges; vs `yim_2024_dauer` = ~99%.
> - `brittin_2021_l4_contactome` — 4,763 undirected → 9,526 edges; vs `white_1986_jsh` = 99.3%.
> - `brittin_2021_adult_contactome` — 5,245 undirected → 10,490 edges; vs `white_1986_n2u` = 97.9%.

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

- **Ingestion route (used):** read each publication's **own supplementary file** directly. Brittin =
  Supplementary Information 1 (raw per-section membrane-contact scoring, ~8.9 MB ASCII CSV), aggregated
  per unordered cell pair; Yim = Supplementary Data 6 (contact-area matrix). This avoids Brittin's
  17 GB Zenodo MySQL *and* the OpenWorm `cect` package entirely. Note the repo's committed graphs
  (`reference_graph_*` = Supp Info 3) are a **thresholded** subset (~900–1,046 contacts, ≥4-dataset
  reproducibility, 173 cells) that fails the synapse⊆contact check — SI 1 is the full raw data.
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

## C. Predicted chemical-synapse sign

> **Status: data + derivation built (increments 1–3); validation (increment 4) pending. Decision:
> derive natively; validate against Fenyves 2020.**

### Why

Part B gives *measured* functional sign, but only for Randi's small head subset. The much broader
signed layer is the **predicted polarity of chemical synapses** — whether each EM synapse is
excitatory or inhibitory — inferred from the presynaptic neurotransmitter and the postsynaptic
ionotropic receptor. This puts a sign on most of the ionotropic chemical connectome (thousands of
edges), complementing B's handful of measured ones, and enables signed-circuit analysis (feedforward
excitation, feedback inhibition, E/I balance). The hermaphrodite layer is also the **validated engine**
for the genuinely novel follow-on — a male signed connectome built by swapping the expression
substrate (see "Follow-on" below).

### Prior art (checked) — this is a *reproduction*, not new science

The hermaphrodite prediction is **already published**, so Part C's value is integration into CIRCE
(sign queryable alongside measured functional sign and NT identity) and a trustworthy, reusable
derivation — not a new biological result.

- **Fenyves et al. 2020** (*PLoS Comput Biol* 16:e1007974) / **EleganSign** — the method and the
  hermaphrodite result. For each chemical synapse, pairs presynaptic NT with postsynaptic ionotropic
  receptor expression: **cation-channel receptor → excitatory, anion-channel → inhibitory**; a synapse
  is excitatory if **only cation** receptor genes match, inhibitory if **only anion**, **"complex"** if
  both, unpredicted if none. Covers ~**73%** of the ionotropic chemical synapses; ~**4:1** E:I.
  - Connectome used: **WormWiring** (3,638 connections / 20,589 synapses) = Cook 2019 hermaphrodite.
  - Receptors: **62** ionotropic postsynaptic genes (Glu/ACh/GABA); **42** found expressed; `acc-4`
    and `lgc-46` excluded (presynaptic localization). Per-gene cation/anion calls are the paper's
    **S1 Table**.
  - Inputs: older NT compilations + CeNGEN/Altun receptor expression, binarized.
- **Data availability (confirmed downloadable):** **S1 Data** (per-connection predicted sign — the
  comparison ground truth) at `linkgroup.hu/docs/S1_Data.xlsx` (xlsx, ~1.0 MB; also on PLoS SI), and
  **S1 Table** (receptor → cation/anion mapping) on PLoS SI (`.docx`).
- **Signed-motif analysis** (Csermely group, bioRxiv 2025, PMC13374236) relabels Fenyves: of 3,638
  connections, 1,555 excitatory / 553 inhibitory / 1,530 unknown. Downstream, not new primary data.
- Direct single-synapse polarity measurements remain sparse, so these are *predictions*.

### Approach — derive natively (decided)

Mint a per-synapse polarity exactly as CIRCE derives the neuropeptide / monoamine connectomes from
expression × pairs, rather than ingesting EleganSign's output. Native derivation is the CIRCE pattern,
is reusable for the male swap, and keeps EleganSign as an **independent** validation target (ingesting
its table would make the "validation" circular). The engine:

> For each chemical synapse `pre → post`: let `T = NT(pre)` (Wang 2024). Collect the ionotropic
> receptor genes for transmitter `T` that `post` expresses (CeNGEN). The matched genes being
> **cation-only → `excitatory`**, **anion-only → `inhibitory`**, **mixed → `complex` (sign unset)**,
> **none → unpredicted (sign unset)**.

### Inputs

| Need | Status |
|---|---|
| Chemical connectome to sign | ✅ `cook_2019_hermaphrodite` (= WormWiring; validation target); White / Witvliet as bonus coverage |
| Presynaptic NT identity | ✅ Wang 2024 `neurotransmitter_assignments` (ACh / Glu / GABA + monoamines) |
| Postsynaptic ionotropic receptor expression | ➕ add CeNGEN expression for the fast-transmitter ionotropic receptors |
| Receptor → cation/anion mapping | ➕ curate from Fenyves **S1 Table** (cited) |
| NT → receptor-family map | ➕ small curated table (ACh→nAChR/ACC; Glu→iGluR; GABA→GABA-A-type) |

### New data to add (clone the `cengen-monoamine-ionotropic` layer — exact template)

1. `data/cengen-ionotropic-receptors/ionotropic_receptor_genes.csv` — the fast-transmitter ionotropic
   receptor genes (ACh `acr-*`/`acc-*`/`lgc-*`; Glu `glr-*`/`glc-*`/`avr-*`; GABA `unc-49`/`exp-1`/
   `lgc-35`), each with **transmitter**, **cation/anion** (reproducing Fenyves S1 Table, cited), and
   WormBase gene id; `acc-4`/`lgc-46` excluded with a noted reason.
2. CeNGEN expression for those genes → `ionotropic_receptor_expression.csv` (same ingestion as the
   monoamine ionotropic layer: `build_expression.py` + `MANIFEST.md`).
3. `nt_receptor_family.csv` — maps each presynaptic transmitter to the receptor families that gate it.

### Schema / build design

- **Reuse the `sign` slot** (added in B) on the chemical connections. Distinguish predicted-chemical
  from measured-functional sign by `connection_type` (chemical vs functional) and by a derived
  provenance dataset (e.g. `fenyves_derived_synapse_sign`). "complex" / unpredicted leave `sign` unset
  (a provenance note records *which* of the two it was).
- Measured (Randi functional) and predicted (chemical) sign coexist cleanly; no viz change required
  unless we choose to surface it later.

### Validation plan (the point of doing the hermaphrodite layer)

Compare our native predictions to Fenyves **S1 Data**, decomposing any divergence:

1. **Logic check (same inputs).** Feed Fenyves' own NT + receptor calls through our engine; we must
   reproduce their per-connection labels (near-exact). A mismatch here is a bug in our rule.
2. **Data-vintage axis (our inputs).** Swap in Wang 2024 NT + current CeNGEN; report where and why a
   sign changes. These differences are *data*, not errors — attribute them.
3. **Headline sanity (per CLAUDE.md):** coverage ≈ 73%, E:I ≈ 4:1, "complex" rate in range. A result
   outside these is a red flag to chase before trusting the layer.
4. **Internal cross-checks:** concordance with Part B measured functional sign where they overlap, and
   with NT identity (GABAergic source → inhibitory trend). Report honestly (functional sign is a net
   circuit effect, not a monosynaptic guarantee).

Downloading S1 Data (~1 MB, `linkgroup.hu`) happens at this step, with explicit go-ahead; it is a
**validation reference only** — not vendored, not the source of CIRCE's signs.

### Design decisions (recommended defaults)

1. **Validation target connectome:** `cook_2019_hermaphrodite` (= WormWiring). Sign White / Witvliet
   too for coverage, but validate against Cook only.
2. **Cation/anion source:** reproduce Fenyves S1 Table as a curated, cited CIRCE file — not ingest
   EleganSign's output (keeps validation independent).
3. **"complex" / unpredicted:** `sign` unset, the two distinguished by a provenance note.
4. **NT source:** Wang 2024 (newer than Fenyves'); treat the difference as a measured axis, not an
   error.

### Caveats

- Predictions, not measurements; the cation=excitatory / anion=inhibitory mapping is an assumption;
  "complex" synapses are unresolved; gene-expression binarization thresholds apply; metabotropic
  (GPCR) transmission and monoamine / neuropeptide sign are out of scope here (fast ionotropic only).

### Deliverables / increments

1. ✅ Curate the ionotropic receptor gene list + cation/anion mapping (data only). *(PR #97 — 64
   S1-Table rows = 62 WBGene, grounded offline against WS298. No separate `nt_receptor_family.csv`
   was needed: each gene already carries its `transmitter`, so the NT→family relation is the gene
   table itself.)*
2. ✅ Ingest CeNGEN ionotropic receptor expression (new layer; template = monoamine ionotropic).
   *(PR #98 — published 021821 threshold-4 matrix, per-class → per-cell expansion; 58/62 genes with
   neuronal expression, 290 neurons.)*
3. ✅ Build step: derive per-synapse sign, populate `sign` + provenance; tests. *(This increment —
   the `fenyves_derived_synapse_sign` provenance dataset writes `sign` + a `predicted_sign` evidence
   note on the hermaphrodite, non-dauer chemical connections; pure engine unit-tested.)*
4. 🔶 Validation analysis under `analysis/partc_sign_validation/` (mostly done). Against two verified
   references: **Fenyves 2020** (per-synapse) and **Szilágyi et al., *BMC Biol* 2026;24:167** (DOI
   10.1186/s12915-026-02641-4; Fenyves' own group), which relabels the 2020 complex synapses to one
   sign **by the anion/cation receptor-count ratio** (ties → unknown; 3,638 → 1,555 / 553 / 1,530).
   **Findings:** the engine is self-consistent (4,294/4,294); under the authors' own complex
   resolution our result reproduces the hermaphrodite balance (E:I ≈ 3.08 vs their 2.81–3.52,
   definite coverage ≈ 65%) — the raw E:I ≈ 1.5 was a complex-convention artifact. The **biological
   reality check vs measured functional sign (Randi 2023)** passes at 64% (resolved) but exposes a
   real limitation: **glutamatergic synapses are over-called inhibitory**, driven by the broadly
   threshold-4-expressed GluCls `avr-14`/`glc-4` — so ACh/GABA predictions are trustworthy,
   glutamatergic ones are low-confidence. **Still pending (needs a download):** the exact per-synapse
   logic check + data-vintage decomposition vs Fenyves S1 Data, and a CeNGEN-threshold sensitivity
   test for the GluCl over-call.

> **Increment-3 build result, and what increment 4 must actually explain.** On Cook 2019
> hermaphrodite (Fenyves' target connectome), over the in-scope chemical connections (presynaptic
> fast transmitter present): 704 excitatory, 465 inhibitory, 1,842 complex, 1,283 unpredicted. The
> engine is correct (unit-tested; textbook AVA→VA excitatory, RIS→AVE inhibitory reproduce); the
> differences from Fenyves are **data**, not a rule bug. Note the comparison pitfalls before reading
> any single number:
> - **Count the same way.** Fenyves' headline "73% predicted" counts *complex as predicted*
>   (excitatory + inhibitory + complex). Counted that way our predicted fraction is ≈70% — i.e.
>   **coverage is NOT the main divergence**. The ≈27% figure is the *definite-sign-only* fraction, a
>   different quantity. (Fenyves also counts **synapses**, we count **connections**, and per-synapse
>   vs per-connection weighting shifts the split — reconciled in increment 4.)
> - **Fenyves is neuron→neuron too.** Their network is 20,589 synapses among 297 neurons; **muscle is
>   not scored**. So the ~1,000 body-wall-muscle NMJs that fall to `unpredicted` here (post absent
>   from CeNGEN) are an *extension beyond* Fenyves, not a gap against it — they appear only because we
>   sign all hermaphrodite chemical synapses rather than restricting to the neuron→neuron scope.
>   Adding a muscle receptor-expression substrate (Altun / WormAtlas) is an optional extension, not a
>   fix needed to match Fenyves.
> - **The genuine divergence is the complex rate and E:I ratio, from expression density.** Fenyves
>   2020 (per-synapse) is ≈3.5:1 E:I with 16.6% complex; ours is ≈1.5:1 with far more complex. This
>   is the "data-vintage axis": our raw CeNGEN threshold-4 calls the anion channels broadly (raw
>   matrix: avr-14 68/128 classes, unc-38 92/128, lgc-47 45/128, unc-49 38/128) and finds 58/62 genes
>   expressed, versus Fenyves' hand-curated binary matrix (WormBase + WormAtlas + CeNGEN, 42/62
>   expressed). Fenyves themselves report that *adding* CeNGEN lowered their E:I (3.7–4.1 → 3.1–4.1) —
>   we see a stronger dose of the same effect. There is **no independent "~4:1" prior**; ~4:1 is
>   Fenyves' own measured value, and the per-connection split depends on a complex-resolution choice
>   (see below).
> - **"Complex" resolution is a documented, author-endorsed step.** Fenyves' own group (Szilágyi
>   2026) resolves complex to a definite sign by receptor-count ratio. CIRCE therefore keeps the
>   **stored** chemical `sign` raw (complex unset); a receptor-ratio-resolved *view* is an
>   increment-4 derived comparison, not baked into the KG.

### Follow-on (out of scope here): a male signed connectome

The genuinely novel extension — no male signed connectome exists. Same engine, male substrate:
`cook_2019_male` synaptic connectome × male ionotropic receptor expression. The male ionotropic
receptors sit in Morillo 2025's genome-wide `mmc12` (CIRCE currently ingests only its NPP-GPCR rows),
so the substrate is obtainable by the same pipeline. Deferred until the hermaphrodite engine is
validated.

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
- Fenyves BG, et al. *Synaptic polarity and sign-balance prediction using gene expression data in
  the C. elegans chemical synapse neuronal connectome network.* PLoS Comput Biol 16(12):e1007974
  (2020). (EleganSign; S1 Data = per-connection predicted sign, S1 Table = receptor cation/anion.)
- Szilágyi GS, et al. *Signed motif analysis of the C. elegans neuronal network reveals positive
  feedforward and negative feedback loops.* bioRxiv 2025.01.09.632090 (2025; PMC13374236). (Relabels
  the Fenyves polarity dataset.)
- Wang J, et al. *A neurotransmitter atlas of C. elegans males and hermaphrodites.* (Wang 2024;
  presynaptic transmitter identity.)
