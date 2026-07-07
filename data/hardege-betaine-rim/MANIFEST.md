# Pinned reference — Hardege et al., betaine neurotransmission (RIM)

Evidence for **betaine** as a neurotransmitter of the **RIM** interneurons — the basis for the
`RIM` betaine assignment in [`../curation/neurotransmitter_curation.csv`]
(../curation/neurotransmitter_curation.csv) and for adding `b` (betaine) to the transmitter
vocabulary.

## Source

- **Paper:** Hardege I, Morud J, Yu J, Wilson TS, Schroeder FC, Schafer WR. *Neuronally produced
  betaine acts via a ligand-gated ion channel to control behavioral states.* PNAS 119(48), 2022.
- **DOI:** 10.1073/pnas.2201783119 · **PMID:** 36413500 · **PMCID:** PMC9860315
- **Verified (UTC):** 2026-07-07 (abstract via Europe PMC; the PNAS site is Cloudflare-gated to
  headless fetch). Abstract retained in `abstract.txt`.

## Why this qualifies as a transmitter identity (not just uptake)

Meets the release-mechanism bar we apply to transmitter calls (cf.
[`../wang-neurotransmitter-atlas/`](../wang-neurotransmitter-atlas/MANIFEST.md) and the ASI
uptake caveat):

- **Synthesis** — betaine is *produced* in the RIM interneurons.
- **Vesicular packaging** — loaded into synaptic vesicles by the vesicular monoamine
  transporter **CAT-1 / VMAT** expressed in RIM.
- **Receptor + function** — acts on a betaine-gated chloride channel **LGC-41**; betaine-synthesis
  mutants mis-control the local↔global foraging switch, rescued by restoring betaine specifically
  to RIM.

Contrast: ASI's `betaine (uptake)` (Wang 2024) is uptake-only and remains **unassigned**; RIM's
betaine is synthesized and released, so it is applied.

## Applied

`RIML, RIMR`: `lt → blt` (glutamate + tyramine + **betaine**). This is the first use of the `b`
code; see the betaine-vocabulary note in `../curation/MANIFEST.md`.
