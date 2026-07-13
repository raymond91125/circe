# Pinned input — neuron-graph

Source: `mattpocock`-style sibling repo at `/home/raymond/local/src/git/neuron-graph`
(FunCoNN, a NemaNode fork).

- Upstream commit: `6ea5f9b55715e68eee4089dec2b3a4c2574d23e2`
- Pinned: 2026-06-23 (functional TSVs); re-pinned 2026-06-25 (populate-database raw data)

These are the actual inputs the neuron-graph `populate-database` step loads into MySQL
(`src/server/populate-db/raw-data/`). neuron-graph and WBBT remain external upstreams; this
repo pins a versioned snapshot for reproducibility but does not fork them.

## Cell list and datasets

| File | sha256 | notes |
|------|--------|-------|
| `neurons.json`  | `c780e0b7cab28b311a7e0916a7b8e37b9cf72de152712bd4a5eb29fd9d595492` | 447 cells (neurons + 95 BWM muscles). Fields: `name`, `classes`, `nt`, `typ`, `emb`, `inhead`, `intail`. `typ` is a NemaNode neuron-class code (s/i/m + combos; `b`=muscle; `u`/`sim`/… = pharyngeal/other), not a clean neuron/muscle flag — classification is deferred to the match/build stage. |
| `datasets.json` | `516b564f48f634a440d8e53917649c92dae1a9e97d72be6a9ca6f5308a66dfe5` | 15 datasets. Fields: `id`, `name`, `type`, `time`, `visualTime`, `description`, `datatypes` (`cs`/`gj`/`fc`). |

## Connections

Each file is a JSON array of `{pre, post, typ, syn, ids?, pre_tid?, post_tid?}`. The
`datasetId` is **not** in the records — it is the file stem (e.g. `white_1986_whole`).

Connection `typ` encoding (from `populate-connections.js`): **0 = chemical**,
**2 = electrical / gap junction**, **4 = functional**. Weight derivation mirrors
neuron-graph: functional → `round(sum(syn))`; chemical/electrical → `len(syn)` (contact
count). Gap junctions are symmetric and reverse-deduplicated by the viz — a build-stage
concern, not applied at ingest.

| File | sha256 | chem / gap / func | n |
|------|--------|-------------------|---|
| `connections/white_1986_whole.json` | `fe7ff4a6529ccb9d13183e665644e30ad8aded887c64b24732808757bddae41d` | 2405 / 1144 / 0 | 3549 |
| `connections/white_1986_jsh.json`   | `ac251c32cfbc3f3ee5add0f4f566c3006b715553e822210cc4ec9f1fc9816dd4` | 1480 / 576 / 0 | 2056 |
| `connections/white_1986_n2u.json`   | `ba949e1e1766020009d0db38aa3182b245aa6242bb0359a253cd65d0a6b9a5d4` | 1629 / 538 / 0 | 2167 |
| `connections/white_1986_jse.json`   | `4a77985e4038eb0ddf06f2eab9561c83be3d9a74cd07eb4328ee6023dba8fa15` | 104 / 0 / 0 | 104 |
| `connections/witvliet_2020_1.json`  | `50b3e052e55dbd4578f8faf8f171f9e586bda9f183702d53228cacd856c06d78` | 775 / 164 / 0 | 939 |
| `connections/witvliet_2020_2.json`  | `dd0f45f0762eb921fdfbc2d4a6e35f9abaeb9aac636665fd57d655bc865922c3` | 986 / 245 / 0 | 1231 |
| `connections/witvliet_2020_3.json`  | `3ba7f6619c1693633ae05697b31912c0da3a40cd8b4e27bb2406a04a2fa0adfd` | 1012 / 186 / 0 | 1198 |
| `connections/witvliet_2020_4.json`  | `4394d09c16318be6bcf503ef3d76a9f030e13ea691ca2e4af4c488901c727f7f` | 1136 / 415 / 0 | 1551 |
| `connections/witvliet_2020_5.json`  | `54957f08f53b48c1cf6f11502b40b76d6dc2f8ce84edd1af012c9a7c54dfab69` | 1515 / 577 / 0 | 2092 |
| `connections/witvliet_2020_6.json`  | `9068c14e4c40b689329d5eabfdb53fd46b374d4606ddb67a250bf7798e49dd4d` | 1525 / 427 / 0 | 1952 |
| `connections/witvliet_2020_7.json`  | `350b247214763fe6b5a89b76a95cc018f7551c0446ad6501bd9c4774fb05fbae` | 2202 / 579 / 0 | 2781 |
| `connections/witvliet_2020_8.json`  | `465aeb3c1f42d84273023cef8e5f5d8ebbcf63248550b2219db68fccb9787033` | 2186 / 616 / 0 | 2802 |
| `connections/randi_funconn_wildty.json` | `6f3c8d6e40411c9617c79542b1ba251258541b0bedb45aaca7f254905fb72c4f` | 0 / 0 / 1131 | 1131 |
| `connections/randi_funconn_wildcp.json` | `6f3c8d6e40411c9617c79542b1ba251258541b0bedb45aaca7f254905fb72c4f` | 0 / 0 / 1131 | 1131 |
| `connections/randi_funconn_unc31.json`  | `a52d23f9e7dd27a0ec96250e77ba7bc41ea062ca6334612a7267ee94d797e812` | 0 / 0 / 352 | 352 |

**`wildcp` ≡ `wildty`:** identical sha256 — the upstream ships the same wild-type functional
data under two dataset labels (they map to neuron-graph's "complete" vs "head" database
collections). Ingest reads both faithfully; the **build drops `randi_funconn_wildcp`** as a
duplicate (see `assemble._REDUNDANT_NG_DATASETS`) so KG aggregates don't double-count wild-type
functional weights. `randi_funconn_wildty` is kept.

## Superseded

The Phase 0 `randi_funconn_*.tsv` files (Randi correlation matrices including zero entries,
from neuron-graph's `functional/` dir) were removed: the canonical populate-database inputs
are the `connections/randi_funconn_*.json` files above, which carry only realized
connections and share the uniform record shape used for all three connection types.

## Provenance gap: neuron `nt` (neurotransmitter) and `typ` fields

The `neurotransmitter` (`nt`) and cell-`typ` values in `neurons.json` are inherited verbatim
with **no cited per-neuron source** in this data's lineage. Traced 2026-07-02:

- **Lineage:** this pin ← neuron-graph (FunCoNN fork) ← **NemaNode** (`zhenlab-ltri/NemaNode`,
  Zhen Lab). `neurons.json` is byte-identical across all three.
- **Origin commit:** in NemaNode, `neurons.json` was first committed **2020-05-14** by Daniel
  Witvliet as the initial commit ("after removing unpublished datasets") — no pre-history, no
  provenance recorded. NemaNode's in-app "Data sources" and docs cite only the *connectivity*
  reconstructions, never the neurotransmitter assignments.
- **Publication (Witvliet et al., 2020):** the Methods "Classification of neuron types"
  section sources only the **monoaminergic** identities — serotonin (AIM, HSN), dopamine
  (ADE, CEP), octopamine (RIC) — cited to **Sulston, Dew & Brenner 1975** and **Duerr et al.
  1999**, and only to justify the *modulatory* `typ`. Neuron **type** classification follows
  **White et al. 1986** (Table S1). The majority `nt` assignments — **acetylcholine (`a`),
  GABA (`g`), glutamate (`l`)** — are **not cited anywhere** in the lineage or the paper.
- **Likely (uncited) origin:** the field-standard C. elegans neurotransmitter atlas (Hobert
  lab; Pereira 2015 cholinergic, Gendrel 2016 GABAergic, Serrano-Saiz 2013 glutamatergic) via
  WormAtlas/WormBase — but this cannot be confirmed from any artifact in the lineage.

Treat `nt` (and `typ`) as unsourced upstream annotations; to establish provenance, ask the
NemaNode maintainers (Zhen Lab).

## Not pinned (out of v1 scope)

- `annotations/*.json` — connection annotations (head annotations exist; not modeled in v1).
- `trajectories/*.json` — 3D morphology (multi-MB, disabled in neuron-graph production).
