# Neurotransmitter reconciliation: neurons.json vs WormAtlas vs Wang et al.

Cross-check of the neuron-graph `nt` field (carried into the KG, provenance-less — see
`data/neuron-graph/MANIFEST.md`) against two cited sources:

- **WormAtlas** neurotransmitter table — curated "By Class" summaries (`data/wormatlas-neurotransmitter/`)
- **Wang et al. 2024** eLife atlas — Table S2 consolidated calls (`data/wang-neurotransmitter-atlas/`)

Reproduce: `uv run --with 'beautifulsoup4,lxml,openpyxl' python analysis/reconcile_neurotransmitters.py`
→ writes `neurotransmitter_reconciliation.csv` (per-class, all three sources).

**117 neuron classes compared; the three sources agree on ~101.** The rest split into three
kinds — only one is a real `neurons.json` error.

## 1. Robust discrepancy — a likely `neurons.json` error
| class | neurons.json | WormAtlas | Wang | verdict |
|---|---|---|---|---|
| **HSN** | glutamate, serotonin | **ACh**, serotonin | **ACh** | Both curated sources say ACh (+serotonin); neither supports **glutamate**. `neurons.json`'s glutamate is almost certainly wrong — HSN is serotonergic + **cholinergic**. |

## 2. `neurons.json` gaps that the atlas fills (`neurons.json` has no transmitter; Wang assigns one)
`ALA`→GABA, `AVF`→GABA, `AVJ`→GABA, `AWA`→ACh, `I4`→glutamate, `RIP`→ACh, `ASI`→betaine (uptake).
These are cells neuron-graph left blank (`nt` empty / `u`); Wang provides an identity (several
are recent/orphan calls). Worth adopting if the KG wants transmitter coverage for them.

## 3. Not real conflicts — source granularity / known caveats
- **Wang lists the primary (fast) transmitter**, so co-transmitter neurons look like they
  "lost" a transmitter vs the dual-listing sources: `ADF`, `RIH`, `VCn` (ACh **+ serotonin**
  in neurons.json/WormAtlas; Wang shows ACh), `RIM` (glutamate **+ tyramine**; Wang glutamate),
  `RIB` (ACh **+ GABA**; Wang GABA). The monoamine/second transmitter is recorded on other
  sheets of the Wang supplement, not the consolidated column parsed here.
- **AIM** — WormAtlas puts AIM in the cholinergic summary, but that ACh is the **male, post-L3**
  phenotype (WormAtlas footnote 4; Pereira 2015). In the hermaphrodite AIM is
  glutamatergic/serotonergic, matching `neurons.json`. Not a discrepancy.
- **I5** — `neurons.json` adds serotonin; both WormAtlas summary and Wang give glutamate only
  (I5 serotonin is *weak* antibody signal, WormAtlas footnote 35). `neurons.json` likely
  over-asserts; low-confidence.
- **ASG** — serotonin appears only as a footnote-qualified WormAtlas mention; Wang says
  glutamate. Minor.

## Caveats (why this is a cross-check, not ground truth)
- **Wang consolidated column = primary/fast transmitter**; co-transmitters (serotonin, other
  monoamines) are under-counted here → category 3 above.
- **WormAtlas "By Class" mixes hermaphrodite + male** and encodes sex / developmental-stage /
  evidence-strength in `(m)`/`(h)` markers and footnotes; the matcher drops `(m)` male and
  non-neuronal entries but does **not** parse footnotes (hence the AIM caveat).
- Matching is class-level; ventral-cord motor classes (`DAn`…) matched by stripping the `n`,
  ranges (`DD1-6`) at the class prefix.

**Bottom line:** neuron-graph's `nt` agrees with the modern atlases across the vast majority of
classes. The single clear correction is **HSN (glutamate → ACh)**; secondarily, `neurons.json`
omits transmitters for a handful of cells the atlas now identifies (category 2).
