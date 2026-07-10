# Cook et al. 2020 — C. elegans pharyngeal connectome

## Source
Cook SJ, Crouse CM, Yakovlev MA, Nguyen KCQ, Hall DH, Emmons SW (2020).
**The connectome of the *Caenorhabditis elegans* pharynx.** *J Comp Neurol* 528:2767-2784.
PMC7601127 · doi:10.1002/cne.24932 · data also at https://www.wormwiring.org

## Vendored data
`SI3_combined_synapse_list.csv` — **Supplemental Data 3**, the per-synapse list exported from the
Elegance reconstruction databases. Each row is one combined synapse: `pre` -> `post1..post4`
(polyadic), `type` (chemical/electrical), `series` (N2W/JSA/N2T), and `sections` (EM
serial-section count). sha256: `8ddc0e2a35a2c50ea4ab61a162b6aeffa7929d2f5d4f9295cb2bd9db3796a933`

`edges.csv` — the weighted edge list produced by `aggregate_edges.py` (regenerate with
`python aggregate_edges.py`). 332 edges (276 chemical, 56 electrical).

## Why SI 3 (not SI 4 or the PDF)
- **SI 4** (published edge list CSV) is an internally inconsistent export: 36 directed pairs are
  listed twice with per-series weights the file never averages, and some `0.5` values don't match
  the source DB. Not trustworthy for weights.
- The **supplement PDF** renders the same synapse list with glued / `########` weight fields.
- **SI 3** is the clean raw synapse list (matches the Elegance DBs) — the authoritative source.

## Weight semantics
weight = **total EM serial sections** summed over every synapse connecting a pair (each polyadic
pre->post_i pairing counts the synapse's `sections`), across all series. This is the SAME
definition as Cook 2019 (synapse number x size), so the pharynx is directly comparable to the 2019
pharyngeal subset. We deliberately do NOT reproduce SI 4's per-series average.

## Raw provenance / cross-check (not vendored — large)
Elegance MySQL dumps, pinned by sha256:
- SI 1 `cne24932-sup-0001-supinfo1.sql` (database `n2w`): `abc82e86615a2d7f429e27300c0fd4285e43aca073d179e728caeca9c4ec46b4`
- SI 2 `cne24932-sup-0002-supinfo2.sql` (database `jsa`): `fb3c005d02a90e76b066596ea9e63e684963182372489d5ac3882b6f575ab991`

## Curation notes
- Non-cell endpoints dropped: `obj560962`, `obj586937`, `unk`, `unk1`.
- Lowercase name variants reconciled via `cook_name_aliases.csv` (pm4d->pm4D, pm5d->pm5D, mc3v->mc3V, ...).
- `g1vl`/`g1vr` (non-standard names for the g1 gland ventral processes) are compiled into the
  `g1` gland cell via `cook_name_aliases.csv` (g1vl/g1vr -> g1); `g1` is grounded to the
  pharyngeal g1 gland term WBbt:0003712 in `cook_anatomy_curation.csv`.
