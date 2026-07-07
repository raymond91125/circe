# Anatomy curation

## `neurotransmitter_curation.csv` — corrections to neuron-graph's `nt`

Evidence-based overrides of the (provenance-less) neuron-graph `neurotransmitter` code,
keyed by cell name (`cell_name, neurotransmitter, note`). Applied in `build`
(`cckg build --nt-curation`), taking precedence over the ingested `nt`. Derived from the
three-way reconciliation in [`../../analysis/neurotransmitter_reconciliation.md`]
(../../analysis/neurotransmitter_reconciliation.md).

### Neurotransmitter code vocabulary

neuron-graph's `nt` codes: `a` ACh, `d` dopamine, `g` GABA, `l` glutamate, `o` octopamine,
`s` serotonin, `t` tyramine, `u` unknown, `n` none. This curation adds **`b` = betaine**
(RIM, per Hardege et al. 2022 — see [`../hardege-betaine-rim/`](../hardege-betaine-rim/MANIFEST.md)),
which also required adding betaine to the neuron-graph viz vocabulary (legend + colour + display).

### Criterion for applying a correction

A transmitter identity requires a **release mechanism**, not just presence of the transmitter:
biosynthesis and/or the **vesicular transporter** — `unc-17`/VAChT (ACh), `unc-47`/VGAT (GABA),
`cat-1`/VMAT (monoamines), `eat-4`/VGLUT (glutamate). Uptake alone (a plasma-membrane
transporter, e.g. `snf-11` GAT) can be clearance. We therefore apply only the Wang et al. 2024
atlas's **non-hedged** calls, and hold its own hedged ones. (Verified from Table S2 reporter
cells, which are graded — filled ≠ uniformly robust; many are "very dim and variable".)

**Applied** (atlas clean call):
- **Correction** — `HSNL, HSNR` `ls → as` (Glu+5HT → **ACh**+5HT): unc-17/VAChT + 5-HT
  *synthesis* (Wang comment: "5-HT synthesis (no uptake machinery detected)"; WormAtlas).
  neuron-graph's glutamate is unsupported.
- **Gap-fills** (`nt = u` → assignment): `ALA`→GABA (clean call, anti-GABA staining);
  `RIPL/RIPR`→ACh (unc-17/VAChT newly detected, "identity added").

**Held — documented but NOT applied** (the atlas's own hedged calls; left `unknown`):
- *"potential" identity on very-dim-and-variable reporter*: `AVJL/AVJR` (unc-47 very dim/variable
  → "potential GABAergic"), `AWAL/AWAR` (unc-17 very dim/variable → "potentially… ACh"),
  `I4` (eat-4 dim/variable → "may potentially… glutamate").
- *uptake, no releasing machinery*: `AVF` — Wang `GABA (uptake)`, **no** unc-47/VGAT detected
  (no reporter expression at all) → cannot vesicularly release GABA. `ASI` — Wang
  `*betaine (uptake)`, cat-1/VMAT+ but snf-3 "dim and variable" → "potential"; betaine also has
  no code in neuron-graph's `nt` vocabulary.

Consistency check: serotonin **uptake** neurons like `RIH` *are* accepted as serotonergic —
because RIH expresses cat-1/VMAT (verified: unc-17 + cat-1 + mod-5, no tph-1), i.e. it has the
vesicular machinery to release the 5-HT it takes up. AVF lacks the equivalent (VGAT), which is
why it is not applied.

## `class_anatomy_curation.csv` — cell-class → WBbt (for the viz cell-info link)

neuron-graph's cell-info panel links to WormBase by cell **class**, so the viz needs a
class → WBbt map. `cckg export` builds it: manual curation here (highest precedence), then a
unique strong (label/exact) WBBT match on the class name, then single-cell-class reuse of the
cell's own anatomy. This file curates the 28 classes that don't auto-resolve — the generic
ventral-cord motor-neuron classes (`ASn`/`DAn`/`DBn`/`DDn`/`VAn`/`VBn`/`VCn`/`VDn` → the
`* neuron` class terms), pharyngeal/labial (`I1`,`I2`,`M2`,`M3`,`MC`,`IL1`,`IL2`), `CEPsh`
(cephalic sheath), pharyngeal glands (`g1`/`g2`), `DefecationMuscles` (enteric muscle), and
the positional BWM classes (→ generic `body wall muscle cell`, since they span all 4
quadrants). With this, all 147 classes map. Emitted to `outputs/neuron-graph/anatomy_terms.json`.

## `connection_endpoint_cells.csv` — stub cells for class-level endpoints

11 names appear as connection pre/post but are not in neuron-graph's cell list (class-level
aggregates / pharyngeal structures). Each is mapped to a WBBT **class** term and minted as a
stub `Cell` so the references resolve and carry grounding. These are KG-only (not neuron-
graph cells), so they are **excluded from the viz `/api/cells` projection**.

Columns: `cell_name, wbbt_id, wbbt_label, cell_type, note`. Covers pharyngeal muscles
(`pm2/3/5/6/7`), marginal cells (`mc2/3`), pharyngeal glands (`G1→g1`, `G2→g2`),
`DefecationMuscles→enteric muscle`, and `VAn→VA neuron`. The 3 reconstruction fragments
(`Fragment`, `NR_fragment`, `vncfrag`) are **left unmapped** — they are unidentified EM
processes with no cell/class term, and remain dangling references (acceptable, open-world).

---



`anatomy_curation.csv` — human-resolved cell → WBbt mappings for the cells the lexical
matcher leaves in the work-list (ambiguous + unmatched). Curator-authored, not an upstream
artifact.

- Curated: 2026-06-25, against the pinned WBBT (`data/wbbt/`).
- Method: each term confirmed against the pinned `wbbt.json` (label/synonym), per the
  `external-term-lookup` skill — no guessed IDs.
- 114 rows = the full Phase 2 work-list (95 BWM muscles + 19 others).

## Columns

`cell_name, wbbt_id, wbbt_label, confidence, note`

## Notes on resolutions

- **BWM-\*** (95) — grounded to the **individual** body wall muscle cell term (under
  `WBbt:0006804`), matched via each cell's `DL`/`DR`/`VL`/`VR` synonym: neuron-graph
  `BWM-DL01` → synonym `DL1` → `WBbt:0006235` (lineage name `MSapappp`). All 95 map 1:1.
- **M1 / M4 / M5** — pharyngeal **neurons** (`WBbt:0004488 / 0004467 / 0004465`), correcting
  the lexical false-match to the pm1/pm4/pm5 **muscle** terms (synonyms "m1"/"m4"/"m5").
- **pm2D / pm3D / pm5D** — the dorsal pharyngeal-muscle pair terms (`pm?DL-pm?DR`).
- Cell-vs-tissue: `hyp` grounds to `WBbt:0007846` **hypodermal cell** (not the `hypodermis`
  tissue term `WBbt:0005733`).
- Lower-confidence rows are flagged in the `confidence` column: `LegacyBodyWallMuscles`
  (legacy aggregate → generic body wall muscle cell) and `excgl`.
