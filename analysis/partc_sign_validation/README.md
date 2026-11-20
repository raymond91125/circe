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

## Bottom line

- The **engine is correct** (4,294/4,294 self-consistent) and, under the published complex-resolution,
  **reproduces the hermaphrodite polarity balance** (E:I ≈ 3.1, coverage ≈ 65%). The headline raw
  divergence was a counting/convention artifact, not a bug — consistent with the #99 write-up.
- **Trust ACh and GABA sign predictions; treat glutamatergic predictions as low-confidence**, because
  threshold-4 over-calls the inhibitory GluCls — a limitation the functional reality check exposes
  independently of Fenyves.
- CIRCE keeps the stored sign **raw** (complex unset); the receptor-ratio-resolved labels are provided
  here as a derived view, not written back to the KG.

## Not done (needs data not yet fetched)

- **Exact logic check** (feed Fenyves' own NT + receptor calls through our engine; reproduce their
  per-connection labels) and the **per-synapse data-vintage decomposition** — both need Fenyves' S1
  Data / S6 Table / S7 Data (≈1 MB, `linkgroup.hu` / PLoS SI), a validation reference only, to be
  downloaded with explicit go-ahead.
- **Threshold sensitivity** (does a less stringent/more stringent CeNGEN threshold change the GluCl
  over-call?) — needs the other CeNGEN threshold matrices.
