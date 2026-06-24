# Pinned input — WBBT (C. elegans Gross Anatomy Ontology)

Source: sibling repo at `/home/raymond/local/src/git/c-elegans-gross-anatomy-ontology`.

- Upstream commit: `1d905a5312239c40db938171fffc60bf8f3fa041`
- Artifact build date (upstream): 2025-08-18
- Pinned: 2026-06-23

## Files pinned

| File | sha256 | notes |
|------|--------|-------|
| `wbbt.json` | `b85f49140542c6570dbfbaa861e44f39422c88c6b71e93e89c89d8fdca8099e3` | standard release, obographs JSON — carries `lbl` (rdfs:label) + `meta.synonyms`, the basis for name->WBbt matching |

## Why this artifact

`wbbt.json` (obographs format) is the matcher input: it exposes each term's primary label
and synonyms in a single parse-friendly file. `wbbt.owl` / `wbbt.obo` are also available
upstream if OWL-level reasoning is needed later. Confirm in Phase 2 which synonym scopes
(`hasExactSynonym` vs related/broad) count as a confident match.
