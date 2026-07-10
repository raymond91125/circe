"""Aggregate the Cook et al. 2020 pharyngeal connectome from the combined synapse list (SI 3).

Source: Cook, Crouse, Yakovlev, Nguyen, Hall, Emmons (2020) J Comp Neurol 528:2767-2784;
PMC7601127. Supplemental Data 3 (`SI3_combined_synapse_list.csv`), the per-synapse list exported
from the Elegance reconstruction databases (SI 1 `n2w`, SI 2 `jsa`; SI 3 also carries an N2T
series). Each row is one (combined) synapse: pre -> post1..post4 (polyadic), a type
(chemical/electrical), a `series`, and a `sections` EM-serial-section count.

We aggregate to a weighted edge list, weight = **total EM serial sections** summed over every
synapse connecting a pair (each polyadic pre->post_i pairing counts the synapse's `sections`),
across all series. This is the SAME weight definition as Cook 2019 (synapse number x size in EM
sections), so the pharynx is directly comparable to the 2019 pharyngeal subset. We deliberately do
NOT reproduce SI 4's per-series average (its published CSV export is internally inconsistent).

Non-cell endpoints (raw object ids `obj...`, `unk`) are dropped; lowercase name variants are left
for `cook_name_aliases.csv` to reconcile, including `g1vl`/`g1vr` (non-standard names for the g1
gland ventral processes) which are compiled into the `g1` gland cell.

Usage: python aggregate_edges.py  ->  writes edges.csv (source,target,weight,type)
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent
DROP = re.compile(r"^(obj\d+|unk\d*)$", re.I)
# g1vl/g1vr are non-standard names for g1 gland ventral processes; kept and compiled to the
# `g1` gland cell via cook_name_aliases.csv (g1vl/g1vr -> g1). Only raw object/unk ids are dropped.
EXCLUDE: set[str] = set()


def aggregate() -> list[tuple[str, str, int, str]]:
    edges: dict[tuple[str, str, str], float] = defaultdict(float)
    with (HERE / "SI3_combined_synapse_list.csv").open() as fh:
        for r in csv.DictReader(fh):
            pre = (r["pre"] or "").strip()
            if not pre or DROP.match(pre) or pre in EXCLUDE:
                continue
            typ = (r["type"] or "").strip()
            sections = float(r["sections"] or 0)
            for i in (1, 2, 3, 4):
                post = (r[f"post{i}"] or "").strip()
                if not post or DROP.match(post) or post in EXCLUDE:
                    continue
                if typ == "electrical":
                    a, b = sorted((pre, post))  # undirected, canonical orientation
                    edges[(a, b, "electrical")] += sections
                else:
                    edges[(pre, post, "chemical")] += sections
    return [(a, b, int(w), t) for (a, b, t), w in edges.items()]


def main() -> None:
    edges = sorted(aggregate())
    out = HERE / "edges.csv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "target", "weight", "type"])
        w.writerows(edges)
    chem = sum(1 for e in edges if e[3] == "chemical")
    elec = sum(1 for e in edges if e[3] == "electrical")
    print(f"wrote {out}: {len(edges)} edges ({chem} chemical, {elec} electrical)")


if __name__ == "__main__":
    main()
