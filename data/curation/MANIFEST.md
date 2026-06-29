# Anatomy curation

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
