# Ionotropic fast-transmitter receptors — Fenyves 2020 polarity reference

The curated set of **ionotropic (ligand-gated) postsynaptic receptor genes** for the three
classical fast transmitters — **acetylcholine, glutamate, GABA** — each tagged with its
**ion selectivity** (cation vs anion), which fixes the **sign** a synapse using that
transmitter onto a cell expressing that receptor would carry.

This is the input gene table for **Part C — predicted chemical-synapse sign**
(see `../../docs/proposals/contactome-and-signed-functional.md`). It is **data only**: no
schema, build, or expression here. Expression (CeNGEN) and the per-synapse derivation are
later increments.

## Source: Fenyves et al. 2020, S1 Table

Fenyves BG, Qin P, Ch'ng Q, et al. *Transcriptional profiling of the signaling network of
the C. elegans connectome.* PLoS Comput Biol 16(12):e1007974 (2020).
doi:10.1371/journal.pcbi.1007974. The polarity assignments and per-gene references here are
transcribed **verbatim** from that paper's **S1 Table** ("62 ionotropic postsynaptic
receptor genes for glutamate, acetylcholine, and GABA … + for cation channel, – for anion
channel"). The S1 Table `.docx` is **not vendored** (copyright); the full content is
reproduced below and in `ionotropic_receptor_genes.csv`.

### Polarity → ion → sign mapping

Fenyves' "Polarity" column uses `+` / `−`. The mapping to CIRCE's functional sign:

| Fenyves polarity | ion selectivity | sign of the synapse |
|------------------|-----------------|---------------------|
| `+`              | cation          | **excitatory**      |
| `−`              | anion           | **inhibitory**      |

(The *gene* carries the ion; the *synapse* sign is derived later by matching a presynaptic
transmitter to the receptor the postsynaptic cell expresses. A cell expressing both a cation
and an anion receptor for the same transmitter → **complex** / sign unset; a transmitter with
no matched receptor → **unpredicted** / sign unset. See the proposal.)

## The genes (64 rows = 62 distinct genes)

Counts over the 64 S1-Table rows:

| transmitter   | cation (excitatory) | anion (inhibitory) |
|---------------|---------------------|--------------------|
| acetylcholine | 30                  | 8                  |
| glutamate     | 10                  | 6                  |
| GABA          | 2 (`exp-1`, `lgc-35`) | 8                |

**Acetylcholine — cation (+):** acr-1, acr-2, acr-3, acr-4, acr-5, acr-6, acr-7, acr-8, acr-9,
acr-10, acr-11, acr-12, acr-13, acr-14, acr-15, acr-16, acr-17, acr-18, acr-19, acr-20, acr-21,
acr-23, acr-25, deg-3, des-2, eat-2, lev-8, unc-29, unc-38, unc-63.
**Acetylcholine — anion (−):** acc-1, acc-2, acc-3, acc-4, lgc-46, lgc-47, lgc-48, lgc-49.
**Glutamate — cation (+):** glr-1, glr-2, glr-3, glr-4, glr-5, glr-6, glr-7, glr-8, nmr-1, nmr-2.
**Glutamate — anion (−):** avr-14, avr-15, glc-1, glc-2, glc-3, glc-4.
**GABA — cation (+):** exp-1, lgc-35.
**GABA — anion (−):** gab-1, ggr-1, ggr-2, ggr-3, lgc-36, lgc-37, lgc-38, unc-49.

### Excluded — presynaptic (not usable as postsynaptic sign predictors)

Fenyves flags two genes as presynaptic; they are kept as rows (`excluded = presynaptic`) for
faithful reproduction but must be dropped before predicting postsynaptic sign:

- **acc-4** (ACh, anion) — presynaptic (Fenyves ref [3])
- **lgc-46** (ACh, anion) — presynaptic (Fenyves ref [24])

## WBGene grounding

Every gene is grounded to a WormBase gene id, resolved **offline and deterministically**
against the pinned WormBase release **WS298** (not a live API — same spirit as the
monoamine-ionotropic layer). Resolution route per symbol: `geneIDs` (current CGC/sequence
name) → `geneOtherIDs` (aliases) → `functional_descriptions` (names used only in prose). All
64 resolve; **all are `Live`**; 0 unresolved.

### Validation checkpoint: 64 rows → 62 distinct WBGene

The 64 S1-Table name-rows collapse to exactly **62 distinct genes**, matching Fenyves'
own stated count ("62 ionotropic postsynaptic receptor genes"). Two gene pairs are listed in
the table under two names each (both members agree on transmitter and ion — self-consistent):

| co-referent pair | WBGene | transmitter / ion |
|------------------|--------|-------------------|
| `acr-4` = `des-2`  | WBGene00000955 | ACh / cation |
| `acr-13` = `lev-8` | WBGene00002975 | ACh / cation |

Downstream builds must key on **WBGene**, not symbol, so these de-duplicate naturally.

### Fenyves name → current WormBase CGC name

Six S1-Table names are now non-current (aliases or prose-only names). The CSV keeps Fenyves'
name in `symbol` (faithful to the table) and records the present CGC name in
`current_cgc_name`:

| Fenyves `symbol` | `current_cgc_name` | WBGene | resolved via |
|------------------|--------------------|--------|--------------|
| acr-1  | lev-1  | WBGene00002974 | geneOtherIDs (alias) |
| acr-4  | des-2  | WBGene00000955 | geneOtherIDs (alias) |
| acr-13 | lev-8  | WBGene00002975 | geneOtherIDs (alias) |
| ggr-1  | lgc-57 | WBGene00001586 | functional_descriptions (prose) |
| ggr-2  | lgc-58 | WBGene00001587 | functional_descriptions (prose) |
| ggr-3  | lgc-56 | WBGene00001588 | functional_descriptions (prose) |

### Cross-layer overlap

**`lgc-56` (WBGene00001588)** — Fenyves' `ggr-3` — also appears in the monoamine-ionotropic
layer (`../cengen-monoamine-ionotropic/`), where it is a deorphanized **amine**-gated channel
(tyramine/dopamine/octopamine; Morud et al. 2021). Both layers classify it as an **anion**
(inhibitory) channel, so its sign contribution is robust; Fenyves treats it as a GABA receptor
for the GABA-sign prediction. This is the only gene shared between the two layers.

## WS298 source files (not vendored)

Resolved against WormBase release WS298 (EBI mirror), downloaded transiently; **not vendored**
(large, and the grounding result is baked into the CSV):

| file | sha256 |
|------|--------|
| `c_elegans.PRJNA13758.WS298.geneIDs.txt.gz`      | `d87c69d1cabc03dc3e0d76b98944fecd952734f1010073b8d84c0970cda19346` |
| `c_elegans.PRJNA13758.WS298.geneOtherIDs.txt.gz` | `d549f0f41bc8cd301ddac1ac91f96fcdbb3815dc38ae29cce040bd758f447608` |
| `c_elegans.PRJNA13758.WS298.functional_descriptions.txt.gz` | `ea6acf482a4cffbd5723daab2da304b43c9b1063a913505fab6d77576a92fae1` |

Base URL: `https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS298/species/c_elegans/PRJNA13758/annotation/`

## Vendored data

- **`ionotropic_receptor_genes.csv`** — 64 rows. Columns:
  `symbol` (Fenyves S1-Table name), `wbgene` (`WB:WBGene…`), `current_cgc_name`,
  `category` (`ionotropic_receptor`), `transmitter` (`acetylcholine`/`glutamate`/`GABA`),
  `ion` (`cation`/`anion`), `excluded` (`presynaptic` or empty), `note` (S1-Table comment /
  Fenyves reference markers), `source` (`Fenyves 2020 S1 Table`), `wbgene_source`
  (`WormBase WS298`).
- **`ionotropic_receptor_expression.csv`** — 2,738 rows, columns `cell`, `gene_symbol`, `wbgene`;
  one row per (neuron, receptor) called expressed at CeNGEN threshold 4 (58 genes × 290 cells).
- **`cengen_class_map.csv`** — the auditable CeNGEN-class → CIRCE-member-cells expansion (128 rows).
- **`build_expression.py`** — regenerates the two CSVs above from the gene table, the CeNGEN
  matrix (env `CENGEN_T4`), and `../neuron-graph/neurons.json`. Deterministic.

## Expression (CeNGEN threshold 4)

Which neuron expresses each receptor, needed to predict a postsynaptic cell's response to its
presynaptic transmitter.

### Source: CeNGEN published threshold-4 matrix

`021821_stringent_threshold4.csv` — CeNGEN's **published** (Taylor et al. 2021, *Cell*
184:4329, doi:10.1016/j.cell.2021.06.023) single-cell expression matrix at the most stringent
threshold (4): 13,669 genes × **128 neuron classes**, keyed by `Wormbase_ID`; value > 0 means
the gene is called expressed in that class. **Not vendored** (9.8 MB); the projection below is
baked into the committed CSV.

- URL: `https://cengen.org/storage/021821_stringent_threshold4.csv`
- sha256: `478791be7f43d46e11cb751cf12016ae517d59b34db3325f29703646853df875`

**Vintage note.** This is a *different* CeNGEN vintage than the monoamine-ionotropic layer
(which uses Ripoll's July-2020 `30072020` threshold-4 subset to reproduce Ripoll). This layer
deliberately uses the **published final release** as a single consistent source for all 62
genes — the 2020 subset omits 6 of them, and mixing vintages within one receptor set would be
a confound for the Fenyves validation. The two layers' vintages are independent and each
documented.

### Genes with no CeNGEN neuronal expression (4 of 62)

`acr-10`, `acr-25`, `eat-2`, `lgc-48` are **absent from the matrix entirely** — CeNGEN lists
only genes detected above background in neurons, and these were not (`eat-2` is pharyngeal-muscle;
the others undetected in neurons). They therefore contribute no postsynaptic sign prediction.
58 of 62 genes have neuronal expression.

### Class → cell expansion

The matrix is per **neuron class**; CIRCE cells are individual neurons, so each class is expanded
to its member cells (a class's members share the class's expression — the standard assumption).
`cengen_class_map.csv` records the expansion for **all 128** classes, auditable. Of these, 105
match a CIRCE class token (e.g. `ADA` → ADAL, ADAR), 4 a cell name (`ASEL`, `ASER`, `DA9`,
`VA12`), and 19 are CeNGEN aggregates mapped explicitly (the ventral-cord motor-neuron groupings
`AS`/`DA`/`DB`/`VA`/`VB`/`VC` with their split-out singletons, and `IL2_DV/LR`, `RMD_DV/LR`,
`RME_DV/LR`). Two collapses are irreducibly ambiguous and resolved **conservatively (union)**:

- **`VD_DD`** → all 13 VD + 6 DD neurons. CeNGEN cannot separate the VD and DD GABAergic motor
  neurons, so expression is assigned to all 19.
- **`AWC_ON` / `AWC_OFF`** → both AWCL and AWCR. The ON/OFF identity is stochastic per animal,
  not fixed to left/right, so a gene called in either is assigned to both AWC neurons.

## Per-synapse sign derivation (build; Part C increment 3)

The build (`celegans_connectome_kg.ingest.synapse_sign` + `build.assemble`) consumes the two CSVs
above to populate the `sign` slot on the **hermaphrodite, non-dauer chemical connections**, under the
provenance dataset **`fenyves_derived_synapse_sign`** (which mints no connections of its own). For
`pre → post`: the presynaptic fast transmitter (`decode_fast_transmitters` of the neuron-graph NT
code, refined by Wang 2024 hermaphrodite calls) is matched to the ionotropic receptor genes `post`
expresses — **cation-only → excitatory, anion-only → inhibitory, mixed → complex (unset), none →
unpredicted (unset)**. Every in-scope connection gets a `predicted_sign` evidence note recording the
outcome, so complex/unpredicted stay queryable though `sign` is unset. Dauer is excluded (its
stage-specific remodeling makes adult CeNGEN least applicable); other stages share the standard
adult-expression caveat.

**Comparison to Fenyves (decomposed in increment 4) — read the counting carefully.** On Cook 2019
hermaphrodite the in-scope chemical connections are 704 excitatory / 465 inhibitory / 1,842 complex /
1,283 unpredicted. Three pitfalls before comparing to Fenyves:
- **Coverage is not the main divergence.** Fenyves' "73% predicted" counts *complex as predicted*;
  counted that way ours is ≈70%. The ≈27% definite-sign-only fraction is a different quantity.
  (Fenyves counts synapses; we count connections — reconciled in increment 4.)
- **Muscle is an extension, not a gap.** Fenyves signed a **neuron→neuron** network (20,589 synapses,
  297 neurons); muscle is not scored. The ~1,000 body-wall-muscle NMJs that fall to `unpredicted`
  here (post absent from neurons-only CeNGEN) appear only because we sign all hermaphrodite chemical
  synapses, not because Fenyves covered them. A muscle receptor-expression substrate (Altun /
  WormAtlas) would extend *beyond* Fenyves.
- **The real divergence is complex rate + E:I, from expression density.** Fenyves 2020 (per-synapse)
  ≈3.5:1 E:I, 16.6% complex; ours ≈1.5:1 with far more complex, because raw threshold-4 calls the
  anion channels broadly (avr-14 68/128 classes, unc-38 92/128, lgc-47 45/128, unc-49 38/128) and
  finds 58/62 genes expressed, vs Fenyves' curated binary matrix (WormBase + WormAtlas + CeNGEN,
  42/62). Fenyves themselves report adding CeNGEN lowered their E:I — this is a stronger dose.

There is no independent "~4:1" prior; ~4:1 is Fenyves' measured value, and any per-connection split
depends on how complex is resolved. Fenyves' own group (Szilágyi et al., *BMC Biol* 2026;24:167, DOI
10.1186/s12915-026-02641-4) resolves complex by the anion/cation receptor-count ratio (ties →
unknown), reporting 3,638 connections → 1,555 / 553 / 1,530. CIRCE keeps the stored chemical `sign`
**raw** (complex unset); a receptor-ratio-resolved view is an increment-4 derived comparison. The
engine is unit-tested in isolation (`tests/test_synapse_sign.py`); the quantitative reconciliation
against both references is the increment-4 validation analysis.
