# Validating the Part C predicted chemical-synapse sign

CIRCE's Part C layer (PRs #97–#99) predicts an excitatory/inhibitory sign for each hermaphrodite
chemical synapse from the presynaptic fast transmitter × postsynaptic ionotropic receptor expression
(Fenyves 2020 method). This analysis asks three things, all on `cook_2019_hermaphrodite` (= WormWiring,
Fenyves' own connectome):

1. **Does it reproduce the published hermaphrodite result** (Fenyves 2020; Szilágyi/Fenyves 2026)?
2. **Does it pass a biological reality check** against *measured* functional sign (Randi 2023)?
3. **Where is it trustworthy, and where not** (confidence calibration)?

Run (after `uv run cckg build`): `uv run python -I analysis/partc_sign_validation/validate.py`.
Deterministic (verified identical across `PYTHONHASHSEED`); numbers are captured in `results.json`.
The threshold-sensitivity follow-on (Result 5) is a second script that needs the CeNGEN threshold
1–3 matrices (not vendored): `CENGEN_MATRICES_DIR=/path/to/matrices uv run python -I
analysis/partc_sign_validation/threshold_sensitivity.py` → `threshold_sensitivity_results.json`.

## Prior art (the references this reconciles with)

- **Fenyves et al. 2020** (*PLoS Comput Biol* 16:e1007974) — the method and the hermaphrodite result:
  **per synapse**, neuron→neuron, 20,589 synapses → 9,070 excitatory / 2,580 inhibitory / 3,413
  complex / 5,526 unpredicted; *"we predicted synaptic polarity for 73% of chemical synapses"*;
  E:I ≈ 3.5:1. Receptor expression was a **hand-curated binary matrix** (WormBase + WormAtlas +
  CeNGEN), with 42/62 receptor genes found expressed.
- **Szilágyi et al. 2026** (*BMC Biol* 24:167, DOI 10.1186/s12915-026-02641-4 — **Fenyves' own
  group**) — relabels the 2020 complex synapses to one sign *"based on the ratio of anion- and
  cation-gated postsynaptic receptor expression"* (ties → unknown), reporting **per connection**:
  3,638 → 1,555 excitatory / 553 inhibitory / 1,530 unknown (E:I ≈ 2.8:1).
- **Randi et al. 2023** (*Nature* 623:406) — *measured* functional sign (signed signal propagation),
  CIRCE's Part B. The only independent measured sign in the graph; used as the reality-check yardstick.

## Method choices (and why)

- **Two views of our sign, reported side by side.** *Raw* = the sign as stored in the KG (complex
  left unset — CIRCE's deliberate choice). *Receptor-ratio-resolved* = complex relabelled by the
  majority of matched cation vs anion receptor genes (ties → unknown), **replicating Szilágyi 2026**.
  The resolved view is the only like-for-like comparison to the published connection-level numbers;
  it is an analysis lens only and does **not** change the KG. Reporting both is the sensitivity check
  on the single most consequential method choice (how complex is handled).
- **Functional sign as the reality check.** It is the one *measured* sign available. Caveat, stated
  up front: functional sign is a **net signal-propagation effect** (possibly polysynaptic /
  extrasynaptic), not a monosynaptic ground truth — so concordance *bounds*, it does not *prove*,
  correctness. All 89→134 overlap pairs are nonetheless direct chemical synapses with a measured
  response, which is the closest monosynaptic-ish test we have.
- **Self-consistency guard.** The script recomputes each synapse's raw label from the receptor
  catalogue and asserts it matches the sign stored in the KG (**4,294/4,294 match**) — so the analysis
  is validating the shipped engine, not a re-implementation that could drift.

## Result 1 — once complex is resolved the same way, we ~reproduce the published result

| view | exc | inh | complex/unknown | unpredicted | definite % | E:I |
|---|---|---|---|---|---|---|
| **ours, raw** | 704 | 465 | 1,842 complex | 1,283 | 27.2% | 1.51 |
| **ours, receptor-ratio-resolved** | 2,113 | 686 | 212 unknown | 1,283 | 65.2% | **3.08** |
| Fenyves 2020 (per **synapse**) | 9,070 | 2,580 | 3,413 complex | 5,526 | 56.6% | 3.52 |
| Szilágyi 2026 (per **connection**) | 1,555 | 553 | 1,530 unknown | — | — | 2.81 |

The raw E:I (1.51) looks alarmingly flat, but that is almost entirely the **complex-handling
convention**, not a broken layer: when complex is resolved by receptor ratio — *the authors' own
method* — E:I jumps to **3.08** (between Szilágyi's 2.81 and Fenyves' 3.52) and definite coverage to
**65%** (vs Fenyves' 56.6%). So CIRCE's engine on current CeNGEN, run through the same complex
resolution, roughly reproduces the published hermaphrodite polarity balance. (Counts differ in
magnitude because Fenyves tallies *synapses* and we tally *connections*; the exact per-synapse
reconciliation needs Fenyves' S1 Data — see "Not done".)

## Result 2 — biological reality check: glutamatergic synapses are over-called inhibitory

Predicted chemical sign vs measured functional sign (Randi 2023), on the pairs that have both:

| view | overlap | agreement | predicted-exc / measured-exc | **predicted-inh / measured-exc** | pred-exc / meas-inh |
|---|---|---|---|---|---|
| raw | 55 | 47.3% | 26 | **26** | 3 |
| resolved | 134 | 64.2% | 86 | **36** | 12 |

Resolving complex lifts agreement to 64%. But a clear residual error survives both views: we predict
**inhibitory for synapses the functional data measures as excitatory** (26 raw / 36 resolved). The
anion receptors driving these false-inhibitory calls are overwhelmingly the **glutamate-gated chloride
channels `avr-14` (21) and `glc-4` (16)** — exactly the genes CeNGEN threshold-4 calls most broadly
(`avr-14` in 68/128 classes). Caveat: Randi's *measured-inhibitory* set barely overlaps our signed
pairs (3 raw / 12 resolved of its 244), so this analysis cannot assess inhibitory **recall** — the
well-powered finding is the false-inhibitory **precision** problem, not a recall claim.

## Result 3 — confidence calibration (raw outcome by presynaptic transmitter)

| presynaptic NT | excitatory | inhibitory | complex | unpredicted | reads as |
|---|---|---|---|---|---|
| acetylcholine | 566 | 31 | 1,257 | 918 | **sane** (excitatory/complex, almost no inhibitory) |
| GABA | 1 | 96 | 65 | 176 | **sane** (inhibitory/complex, almost no excitatory) |
| glutamate | 147 | **340** | 654 | 215 | **suspect** (skews inhibitory via broad GluCls) |

Cholinergic and GABAergic sign predictions are biologically well-behaved. **Glutamatergic predictions
are the unreliable ones**: glutamate is genuinely mixed in *C. elegans*, but current threshold-4
over-represents the inhibitory GluCls (`avr-14`, `glc-4`), pulling glutamatergic synapses toward
inhibitory/complex and depressing the raw E:I. This is the same effect Fenyves themselves noted
(adding CeNGEN lowered their E:I 3.7–4.1 → 3.1); CIRCE's raw threshold-4 is a stronger dose.

## Result 4 — external validation vs Fenyves 2020 S1 Data (logic check + data-vintage)

Fenyves' S1 Data (xlsx, ~1 MB, `linkgroup.hu`, sha256 in `results.json`; a validation reference, **not
vendored** — give its path via `$FENYVES_S1`) contains their inputs *and* outputs: NT assignments
(sheet 1), the cation/anion table (sheet 2), the receptor-expression matrix (sheet 3), and the
per-connection sign prediction (sheet 5). That lets us separate rule from data cleanly.

**Logic check — the rule is identical.** Running CIRCE's sign rule on Fenyves' *own* receptor-match
booleans (S1 sheet 5) reproduces their per-connection label **3,638 / 3,638 = 100%**. So the engine
logic is exactly Fenyves'; every divergence below is **data**, not rule.

**Output comparison — divergence is expression density, not sign direction.** On the 2,435 connections
both connectomes share and we both make a sign call for, raw agreement is only **29.2%** — but that is
category reshuffling, not sign error:

- When **both make a definite (+/−) call, they agree on direction 87.7%** of the time (316 such pairs).
- The disagreements are dominated by *complex/unpredicted* reshuffling: we call **complex** where
  Fenyves calls a clean sign (562 excitatory + 281 inhibitory), and we make a definite or complex call
  where Fenyves had **no prediction** (298 + 269 + 202). Current CeNGEN threshold-4 simply detects far
  more receptor expression than Fenyves' 2020 curated matrix, which *raises* coverage (fewer "no pred":
  theirs 1,422/3,638 = 39%) *and* inflates complex.
- **Attribution:** of 1,725 disagreements, **1,513 (88%) are expression-driven** and only **212 (12%)
  are NT-source-driven** — confirming the data-vintage axis quantitatively. (Fenyves' own
  per-connection balance: + 1,314 / − 425 / complex 477 / no pred 1,422; definite E:I = 3.09.)

So the engine is Fenyves' engine; feeding it today's denser CeNGEN is what moves the numbers, in the
direction the earlier results already showed (more complex, more inhibitory glutamate calls).

## Result 5 — threshold sensitivity: the GluCl over-call is irreducible, not a tunable artifact

The layer uses CeNGEN's **threshold 4** matrix — CeNGEN's *most stringent* (fewest genes called
expressed). CeNGEN publishes four nested binary matrices (1 liberal → 4 stringent;
expr(t4) ⊆ expr(t3) ⊆ expr(t2) ⊆ expr(t1), **verified** here). The obvious question Result 2/3
raises: does a different threshold tame the broad `avr-14`/`glc-4` expression driving the
inhibitory over-call? Since threshold 4 is already the floor, the only directions available are
*less* stringent. Reclassifying the same 4,294 in-scope hermaphrodite chemical connections at each
threshold (swapping **only** the postsynaptic expression) gives
(`threshold_sensitivity.py` → `threshold_sensitivity_results.json`):

| CeNGEN threshold | expr edges | raw exc/inh/cpx/unp | def % | raw E:I | resolved E:I | glu exc/inh/cpx | **glu E:I** | `avr-14` cells | `glc-4` cells |
|---|---|---|---|---|---|---|---|---|---|
| 1 liberal | 6,086 | 55 / 119 / 3,021 / 1,099 | 4.1% | 0.46 | 5.45 | 2 / 33 / 1,213 | **0.06** | 291 | 168 |
| 2 medium | 4,727 | 128 / 243 / 2,812 / 1,111 | 8.6% | 0.53 | 3.18 | 9 / 148 / 1,091 | **0.06** | 274 | 148 |
| 3 conservative | 3,445 | 491 / 413 / 2,205 / 1,185 | 21.1% | 1.19 | 3.20 | 96 / 312 / 799 | **0.31** | 210 | 110 |
| **4 stringent (shipped)** | 2,709 | 704 / 465 / 1,842 / 1,283 | 27.2% | 1.51 | 3.08 | 147 / 340 / 654 | **0.43** | 170 | 86 |

Every over-call metric is **monotonically best at threshold 4** and worsens as the threshold
loosens: the GluCls spread (`avr-14` 170 → 291 cells, `glc-4` 86 → 168), glutamate E:I collapses
toward zero (0.43 → 0.06), complex inflates (42.9% → 70.4%) and raw definite coverage craters
(27.2% → 4.1%). So the over-call is **not a tunable artifact** — it is intrinsic to how broadly
CeNGEN's single-cell data detect these glutamate-gated chloride channels, and it is already as small
as CeNGEN allows. There is no cut more stringent than threshold 4, so the layer cannot tune it away;
the shipped choice of threshold 4 is the right (most conservative) operating point, vindicated.

The lever that *does* help is **complex resolution, not threshold**: the receptor-ratio-resolved
E:I stays ≈3.1–3.2 (t1's 5.45 is an inflation artifact of near-universal expression) and the
resolved reality-check agreement is flat-to-slightly-higher at looser thresholds (64.2% → 72.2%,
but only because ratio resolution leans on the cation/anion *balance* in denser data while raw
definite coverage is destroyed). This confirms increment 4's conclusion from the other direction:
**how complex is handled dominates the threshold choice.** (Raw reality-check agreement swings
41–50% across thresholds but on 8–55 overlap pairs — too few to read into.)

## Result 6 — per-synapse reconciliation: the counting unit is not the gap

Results 1 and 4 compare **per connection**; Fenyves' headline is **per synapse** (20,589 synapses).
This weights each connection by its synapse count (the connectome `weight` = Cook's aggregated
synapse number) and re-tallies per synapse, restricted to the **neuron→neuron** scope Fenyves
scored (`per_synapse_reconciliation.py` → `per_synapse_reconciliation_results.json`).

First, the scope lines up: our neuron→neuron chemical total is **3,707 connections / 20,905
synapses**, within ~2% of Fenyves' **3,638 / 20,589** — confirming `weight` is a synapse count and
that the muscle NMJs (1,051 conns / 6,193 synapses) really are a separate extension, not overlap.

| view (neuron→neuron) | exc | inh | complex | unknown | unpred | **E:I** | predicted incl-complex |
|---|---|---|---|---|---|---|---|
| ours, per **connection** | 704 | 465 | 1,842 | — | 696 | 1.51 | 81.2% |
| ours, per **synapse** (raw) | 3,061 | 2,301 | 12,921 | — | 2,622 | **1.33** | 87.5% |
| ours, per **synapse** (resolved) | 12,663 | 4,044 | — | 1,576 | 2,622 | **3.13** | 79.9% |
| Fenyves 2020, per **synapse** | 9,070 | 2,580 | 3,413 | — | 5,526 | **3.52** | 73.2% |

The counting unit is **not** the source of the raw E:I gap: weighting per synapse moves raw E:I the
*wrong* way (1.51 → **1.33**), because the strongest (high-synapse-count) edges are
disproportionately **complex** — complex is 42.9% of connections but **61.8% of synapses**
(12,921/20,905). So per synapse the complex-handling convention matters *more*, not less. Once
complex is receptor-ratio-resolved, per-synapse **E:I = 3.13 lands beside Fenyves' 3.52**, now in his
counting unit *and* scope; our higher absolute counts (and lower unpredicted: 2,622 vs 5,526) are the
denser-CeNGEN coverage of Result 4, not a sign-direction difference. This is the cleanest
apples-to-apples statement of the reconciliation: same connectome, same neuron→neuron scope
(~20,900 synapses), same per-synapse unit, same complex resolution → **E:I 3.13 vs 3.52.**

## Bottom line

- The **engine is Fenyves' engine** — running our rule on Fenyves' own data reproduces their
  per-connection labels **100%** (3,638/3,638), and the analysis is self-consistent with the shipped
  build (4,294/4,294). Under the published complex-resolution it **reproduces the hermaphrodite
  polarity balance** (E:I ≈ 3.1, coverage ≈ 65%). The raw divergence is a counting/convention +
  expression-density effect (88% of connection-level disagreements are expression-driven, 12% NT),
  not a bug — and **not a counting-unit artifact either**: per synapse, in Fenyves' own scope, the
  resolved E:I is 3.13 vs his 3.52 (Result 6).
- **Trust ACh and GABA sign predictions; treat glutamatergic predictions as low-confidence**, because
  threshold-4 over-calls the inhibitory GluCls — a limitation the functional reality check exposes
  independently of Fenyves, and one that **no CeNGEN threshold removes** (Result 5: threshold 4 is
  already the most stringent, hence the least over-called, cut available).
- CIRCE keeps the stored sign **raw** (complex unset); the receptor-ratio-resolved labels are provided
  here as a derived view, not written back to the KG.

## Beyond this analysis

Both increment-4 follow-ons are now done (Results 5 and 6). The remaining Part C direction is a
separate layer, not a validation of this one: a **male signed connectome** (`cook_2019_male` ×
male-specific ionotropic receptor expression), which would need a male receptor-expression substrate
(CeNGEN is hermaphrodite) before the same engine could be applied.
