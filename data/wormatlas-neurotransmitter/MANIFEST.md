# Pinned reference — WormAtlas neurotransmitter table

A crawled snapshot of WormAtlas's compiled C. elegans neurotransmitter table, with a
structured parse. This is a **cited neurotransmitter source** — captured to help address the
`nt` provenance gap noted in [`../neuron-graph/MANIFEST.md`](../neuron-graph/MANIFEST.md)
(NemaNode's `neurotransmitter` field is inherited with no cited per-neuron source).

## Crawl

- **URL:** https://www.wormatlas.org/neurotransmitterstable.htm
- **Crawled (UTC):** 2026-07-02T22:21:00Z
- **Method:** headless Chrome DOM dump (`google-chrome --headless --dump-dom`). `curl`/
  `requests` (even with `certifi`) reject the site's TLS chain — it serves an incomplete
  chain that strict verifiers won't complete but browsers do — so the snapshot is the
  rendered DOM Chrome retrieved, not raw server bytes.
- **Title:** "Neurotransmitters Table"
- **`neurotransmitterstable.html` sha256:**
  `c9337b7366bf31330de085a7c44765a018cfd10c5ade06b3010ffefbefcf9fbf`

## Files

| File | What |
|------|------|
| `neurotransmitterstable.html` | the crawled snapshot (source of truth) |
| `parse_neurotransmitters.py` | reproducible parser (`uv run --with 'beautifulsoup4,lxml' python parse_neurotransmitters.py`) |
| `neurotransmitters.csv` / `.json` | parsed evidence rows (126) |

## Parse

The page is one large table split into neurotransmitter sections (acetylcholine, dopamine,
tyramine/octopamine, serotonin, GABA, glutamate), each an evidence table with columns
**Description / Gene / Detection method / Localization / References**. The parser emits one
row per evidence entry, tagged with its neurotransmitter (sections delimited by the
"Summary List of X" rows). Counts: acetylcholine 25, dopamine 26, tyramine/octopamine 11,
serotonin 31, GABA 22, glutamate 11 = **126 rows**.

**Caveats — read the `localization` column as WormAtlas wrote it:**
- Neuron assignments live in `localization`; the parser does **not** tokenize it into
  individual neurons. It contains footnote superscripts (as trailing numbers, e.g. `ALA 44`),
  male-specific `(m)` and hermaphrodite `(h)` markers, ranges (`DD1-6`, `VD1-13`), qualifiers
  (`(weak)`, `[ ... ]`), and non-neuronal cells (muscle, intestine, hypodermis).
- Each row is one line of *evidence* (an assay/reporter with a reference), not a definitive
  assignment; a neuron's identity is supported by the union of rows across methods/labs.
- Deriving a clean neuron→neurotransmitter mapping from this (normalizing names, resolving
  markers, reconciling references) is a downstream curation step, intentionally not done here.

## Attribution

© WormAtlas. This snapshot is retained for research reference; cite WormAtlas (and the
primary references listed in the `references` column) when using these assignments.
