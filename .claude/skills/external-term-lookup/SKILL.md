---
name: external-term-lookup
description: Use this to confirm WBBT (and other external ontology) terms and term IDs. NEVER guess ontology terms or WBbt IDs — always confirm against the pinned ontology or a live lookup. Adapted for celegans-connectome-kg (WBBT grounding) from the go-ontology skill of the same name.
---

The connectome KG grounds every cell in the **C. elegans Gross Anatomy Ontology (WBBT)**.
Anatomy URIs (`WBbt:NNNNNNN`) must be **confirmed**, never guessed — a wrong WBbt ID is a
silent data error. Use this skill when resolving the match work-list
(`outputs/match_worklist.csv`), spot-checking an auto-match, or adding any WBbt ID by hand.

## 1. Search the pinned WBBT snapshot (preferred — reproducible, fast)

The pinned ontology is `data/wbbt/wbbt.json` (OBO-graph JSON, the same artifact the matcher
indexes; provenance in `data/wbbt/MANIFEST.md`). Searching it directly is reproducible
because it is the exact version this build is pinned to.

Substring search across labels + synonyms (pure stdlib — works without any extra install):

```bash
uv run python -c "
import json
q = 'body wall muscle'.casefold()   # <-- your query
g = json.load(open('data/wbbt/wbbt.json'))['graphs'][0]
for n in g['nodes']:
    if 'WBbt_' not in n.get('id','') or n.get('type') != 'CLASS': continue
    if n.get('meta',{}).get('deprecated'): continue
    names = [n.get('lbl','')] + [s['val'] for s in n.get('meta',{}).get('synonyms',[])]
    if any(q in (x or '').casefold() for x in names):
        print('WBbt:'+n['id'].split('WBbt_')[1], '|', n.get('lbl'), '|',
              [s['val'] for s in n.get('meta',{}).get('synonyms',[])])
"
```

Exact name lookup via the project's own matcher index (needs the `match` module on `main`):

```bash
uv run python -c "
from celegans_connectome_kg.match.wbbt import WBBTIndex
idx = WBBTIndex.from_obograph('data/wbbt/wbbt.json')
for h in idx.lookup('AVAL'):           # <-- candidate name
    print(h.curie, h.kind, idx.terms[h.curie].label)
"
```

The matcher distinguishes hit *kinds*: `label`/`exact` are confident; `related`/`broad`/
`narrow` are weak (the reason a term lands in the work-list as ambiguous). Confirm a weak
hit by reading the term's definition before adopting its ID.

## 2. Live lookup via OLS — for terms not in the pinned snapshot

If a term may exist upstream but is missing locally (e.g. the pin predates a newly added
term), confirm against the live ontology. The lightweight option is the OLS REST API:

```bash
curl -s 'https://www.ebi.ac.uk/ols4/api/search?q=pharynx&ontology=wbbt' \
  | python3 -c "import json,sys; d=json.load(sys.stdin); [print(x['obo_id'],'|',x['label']) for x in d['response']['docs']]"
```

OAK (`runoak`) also has an `ols:wbbt` adapter (`uv run --with oaklib runoak -i ols:wbbt
search <q>`). Note: OAK's `obograph:` adapter is **very slow** on the full `wbbt.json`
(loads the whole graph per call) — prefer the §1 search for the local file; reach for OAK
only when you need relationship/tree queries, and convert to semsql first if so.

Treat any live result as a lead only: anything used in the build must exist in the
**pinned** `data/wbbt/wbbt.json`. If the right term is genuinely missing from the pin, that
is a signal to re-pin a newer WBBT release (update `data/wbbt/MANIFEST.md`), not to
hand-enter an ID the snapshot does not contain.
