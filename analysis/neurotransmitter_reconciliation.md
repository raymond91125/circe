# Neurotransmitter reconciliation: neurons.json vs WormAtlas vs Wang et al.

Cross-check of the neuron-graph `nt` field (carried into the KG, provenance-less — see
`data/neuron-graph/MANIFEST.md`) against two cited sources:

- **WormAtlas** neurotransmitter table — curated "By Class" summaries (`data/wormatlas-neurotransmitter/`)
- **Wang et al. 2024** eLife atlas — Table S2 consolidated calls (`data/wang-neurotransmitter-atlas/`)

Reproduce: `uv run --with 'openpyxl' python data/wang-neurotransmitter-atlas/parse_atlas.py` then
`uv run --with 'beautifulsoup4,lxml,openpyxl' python analysis/reconcile_neurotransmitters.py`
→ writes `neurotransmitter_reconciliation.csv` (per-class, all three sources).

**117 neuron classes compared; 36 flagged as not-all-agree.** The full Wang `Neurotransmitter(s)`
block (cols U/V/W — see §2b) raised this from ~16: reading V/W surfaces the co-transmitters an
earlier column-U-only parse dropped. Most of the 36 are **Wang adding a co-transmitter**, not a
conflict — they split into the kinds below; only one is a real `neurons.json` error.

## 1. Robust discrepancy — a likely `neurons.json` error
| class | neurons.json | WormAtlas | Wang | verdict |
|---|---|---|---|---|
| **HSN** | glutamate, serotonin | **ACh**, serotonin | **ACh** | Both curated sources say ACh (+serotonin); neither supports **glutamate**. `neurons.json`'s glutamate is almost certainly wrong — HSN is serotonergic + **cholinergic**. |

## 2. `neurons.json` gaps that the atlas fills (`neurons.json` has no transmitter; Wang assigns one)
`ALA`→GABA, `AVF`→GABA, `AVJ`→GABA, `AWA`→ACh, `I4`→glutamate, `RIP`→ACh, `ASI`→betaine (uptake).
These are cells neuron-graph left blank (`nt` empty / `u`); Wang provides an identity (several
are recent/orphan calls). Worth adopting if the KG wants transmitter coverage for them.

**Applied to the curation overlay** (`data/curation/neurotransmitter_curation.csv`): `ALA`→GABA,
`RIP`→ACh (clean, non-hedged staining/reporter calls). `AVF`/`AVJ`/`AWA`/`I4` were **held** — either
atlas-hedged ("potential") or supported only by an uptake transporter with **no vesicular
transporter** to package the transmitter for release (e.g. AVF has GABA-uptake machinery but no
`unc-47`/VGAT), so import ≠ transmitter identity.

## 2b. Wang Table S2 `Neurotransmitter(s)` block spans cols U/V/W (primary/secondary/tertiary)
The `Neurotransmitter(s)` heading is **three merged sub-columns**, so co-transmitters live in V/W.
`parse_atlas.py` now reads all three (it previously read only U). Reading the full block adds a
transmitter `neurons.json` lacks for ~23 classes. These were triaged by evidence quality and the
**release-mechanism criterion** (a transmitter needs synthesis and/or a *vesicular* transporter —
uptake alone into the cytosol is not enough).

**Applied to the overlay (firmer = un-hedged, minus two hand-exclusions):**

| class | applied `nt` | Wang U/V/W | basis |
|---|---|---|---|
| AFD | `l`→`al` | Glu + ACh (NEW) | newly detected `unc-17`/VAChT |
| DVA | `a`→`al` | ACh + Glu (NEW) | newly detected `eat-4`/VGLUT |
| M5 | `a`→`al` | ACh + Glu | — |
| PDE | `d`→`dl` | DA + Glu (NEW) | newly detected `eat-4`/VGLUT |
| RIC | `o`→`lo` | octopamine + Glu | — |
| MI | `l`→`ls` | Glu + 5-HTP (NEW) | serotonin-precursor synthesis (not uptake) |
| AIB | `l`→`bl` | Glu + betaine (uptake, NEW) | `cat-1`/VMAT + `snf-3` |
| PHC | `l`→`bl` | Glu + betaine (uptake, NEW) | `cat-1`/VMAT + `snf-3` |
| RIS | `g`→`bg` | GABA + betaine (uptake, NEW) | `cat-1`/VMAT + `snf-3` |
| PVN | `a`→`abl` | ACh + Glu + betaine (uptake, NEW) | `eat-4`; `cat-1`/VMAT + `snf-3` |
| **CAN** | `n`→`b` | betaine (uptake) [primary] | un-hedged; CAN was a classical-transmitter orphan |
| **AUA** | `l`→`bl` | Glu + betaine (uptake) | un-hedged |
| **ASI** | `u`→`b` | *betaine (uptake) [primary] | `*`; snf-3 "very dim and variable" |
| **RIR** | `a`→`ab` | ACh + betaine (uptake) | `*`; snf-3 "dim and variable" |

(ASI/RIR carry Wang's `*` but *do* meet the vesicular-transporter bar; the four betaine-primary/
secondary classes were approved as a set.) **Betaine here is uptake-based**, distinct from RIM's
synthesis (Hardege 2022, curated `blt`) — the same pattern Wang encodes for RIH (`ACh` +
`5-HT (uptake)`, treated as serotonergic). Uptake + `cat-1`/VMAT = release-competent.

**Held (not applied):**
- **`*`-hedged additions** — ASG, AVJ, AWA, DVB, I4, M3, NSM, PDA, RMH, SMD, URX, and the
  ventral-cord motor classes `DAn`/`VAn`/`VBn`/`VCn` (all `*betaine (uptake) - NEW`).
- **AVF** (`+GABA`) — uptake transporter but **no `unc-47`/VGAT** to package for release; import ≠
  identity. Un-hedged but excluded on mechanism.
- **PVW** (`+5-HT`) — the 5-HT is the **male** phenotype ("male - 5-HT (alternative)"); the
  hermaphrodite call is orphan. Sex-specific, like the AIM caveat.
- **PDC** (also `cat-1`+`snf-3`) — **male-specific**, absent from the hermaphrodite connectome.

## 3. Not real conflicts — source granularity / known caveats
- **Wang under-counts vs the dual-listing sources for a few classes** even after reading U/V/W:
  `RIB` (neurons.json/WormAtlas ACh **+ GABA**; Wang U/V/W = GABA only) and `ADF`/`VCn` serotonin
  co-transmission are recorded elsewhere in the atlas, not in this block. Here `neurons.json` has
  *more*, so no action (we do not remove transmitters). Note `parse_atlas.py` now reads the full
  U/V/W block — RIH correctly shows `ACh + 5-HT (uptake)`, RIM `Glu + tyramine + betaine`.
- **AIM** — WormAtlas puts AIM in the cholinergic summary, but that ACh is the **male, post-L3**
  phenotype (WormAtlas footnote 4; Pereira 2015). In the hermaphrodite AIM is
  glutamatergic/serotonergic, matching `neurons.json`. Not a discrepancy.
- **I5** — `neurons.json` adds serotonin; both WormAtlas summary and Wang give glutamate only
  (I5 serotonin is *weak* antibody signal, WormAtlas footnote 35). `neurons.json` likely
  over-asserts; low-confidence.
- **ASG** — serotonin appears only as a footnote-qualified WormAtlas mention; Wang says
  glutamate. Minor.

## Caveats (why this is a cross-check, not ground truth)
- **Wang calls are reporter/staining-based**; `*` marks the atlas's own tentative calls and many
  new co-transmitters are uptake-based — the overlay applies these only when a vesicular
  transporter supports release (§2b), so some Wang additions are recorded here but not applied.
- **WormAtlas "By Class" mixes hermaphrodite + male** and encodes sex / developmental-stage /
  evidence-strength in `(m)`/`(h)` markers and footnotes; the matcher drops `(m)` male and
  non-neuronal entries but does **not** parse footnotes (hence the AIM caveat).
- Matching is class-level; ventral-cord motor classes (`DAn`…) matched by stripping the `n`,
  ranges (`DD1-6`) at the class prefix.

**Bottom line:** neuron-graph's `nt` agrees with the modern atlases across the vast majority of
classes. The single clear correction is **HSN (glutamate → ACh)**; secondarily, `neurons.json`
omits transmitters for a handful of cells the atlas now identifies (category 2).
