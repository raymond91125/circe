#!/usr/bin/env python
"""Reconcile neuron neurotransmitter assignments across three sources.

Compares, per cell class:
  - neuron-graph `neurons.json`  (the `nt` field carried in the KG; provenance-less, see
    data/neuron-graph/MANIFEST.md)
  - WormAtlas neurotransmitter table  (curated "By Class" summaries; data/wormatlas-neurotransmitter/)
  - Wang et al. 2024 eLife atlas  (Table S2 consolidated calls; data/wang-neurotransmitter-atlas/)

Writes neurotransmitter_reconciliation.csv and prints a discrepancy summary.

Run:  uv run --with 'beautifulsoup4,lxml,openpyxl' python analysis/reconcile_neurotransmitters.py

CAVEATS — this is an automated cross-check, not ground truth:
  * WormAtlas "By Class" mixes hermaphrodite + male phenotypes and encodes sex / developmental
    stage / evidence-strength in `(m)`/`(h)` markers and footnotes. We drop `(m)` male and
    non-neuronal entries, but do NOT parse footnotes — so stage/sex switches can mislead. Known
    example: AIM is glutamatergic/serotonergic in hermaphrodites but the male post-L3 switch to
    ACh makes AIM appear in the cholinergic summary (WormAtlas footnote 4); it is NOT a real
    hermaphrodite discrepancy.
  * neuron-graph classes for ventral-cord motor neurons are `DAn/DDn/…`; matched by stripping
    the trailing `n`. Ranges (`DD1-6`) are matched at the class-prefix level.
"""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

from bs4 import BeautifulSoup

from celegans_connectome_kg.ingest.neuron_graph import read_cells

ROOT = Path(__file__).resolve().parents[1]
NEURONS = ROOT / "data" / "neuron-graph" / "neurons.json"
WORMATLAS_HTML = ROOT / "data" / "wormatlas-neurotransmitter" / "neurotransmitterstable.html"
WANG_CSV = ROOT / "data" / "wang-neurotransmitter-atlas" / "assignments.csv"
OUT = Path(__file__).resolve().parent / "neurotransmitter_reconciliation.csv"

# canonical single-letter codes
NAME = {
    "a": "ACh",
    "d": "dopamine",
    "g": "GABA",
    "l": "glutamate",
    "o": "octopamine",
    "s": "serotonin",
    "t": "tyramine",
    "b": "betaine",
}


def norm_class(cls: str) -> str:
    """neuron-graph ventral-cord motor classes are DAn/DDn/…; strip the placeholder n."""
    return cls[:-1] if re.fullmatch(r"[A-Z]{2,3}n", cls) else cls


def load_neurons_json() -> tuple[dict[str, set[str]], dict[str, list]]:
    cells = read_cells(NEURONS)
    by_class: dict[str, list] = defaultdict(list)
    for c in cells:
        if c.cell_class:
            by_class[c.cell_class].append(c)
    codes = {
        cls: {ch for c in cc for ch in (c.neurotransmitter or "") if ch in NAME}
        for cls, cc in by_class.items()
    }
    return codes, by_class


def is_neuron(cc: list) -> bool:
    return any(
        (c.nemanode_type or "") != "b" and any(x in (c.nemanode_type or "") for x in "simn")
        for c in cc
    )


def load_wormatlas() -> dict[str, set[str]]:
    """WormAtlas curated 'By Class' summaries -> class-prefix -> nt codes (hermaphrodite)."""
    nt_of = {
        "Cholinergic": "a",
        "Dopaminergic": "d",
        "Tyraminergic": "t",
        "Octopaminergic": "o",
        "Serotonergic": "s",
        "GABAergic": "g",
        "Glutamatergic": "l",
    }
    text = re.sub(
        r"\s+", " ", BeautifulSoup(WORMATLAS_HTML.read_text(errors="replace"), "lxml").get_text(" ")
    )
    entry = re.compile(r"([A-Z][A-Za-z]*\d*(?:-\d+)?(?:/[A-Z])?)\s*(\([^)]*\))?")
    out: dict[str, set[str]] = defaultdict(set)
    for block in re.split(r"SUMMARY\s*-\s*", text)[1:]:
        m = re.match(r"([A-Za-z]+)\s+Neurons", block)
        if not m or m.group(1) not in nt_of:
            continue
        code = nt_of[m.group(1)]
        by_class = re.search(r"By Class[^:]*:\s*(.*?)\s*By Body Region", block)
        if not by_class:
            continue
        for name, paren in entry.findall(by_class.group(1)):
            if paren and re.search(r"\bm\b", paren):  # male-specific
                continue
            if name.lower().startswith(("uv", "gonadal")):  # non-neuronal
                continue
            core = name.split("/")[0]
            rng = re.match(r"([A-Z]+)\d+-\d+$", core)
            out[rng.group(1) if rng else re.match(r"[A-Z][A-Za-z0-9]*", core).group(0)].add(code)
    out["RIC"].add("o")  # octopaminergic block has no "By Class" line; hermaphrodite = RIC
    return out


def load_wang(by_class: dict[str, list]) -> dict[str, set[str]]:
    code = {
        "ach": "a",
        "glu": "l",
        "gaba": "g",
        "da": "d",
        "5-ht": "s",
        "5-htp": "s",  # serotonin synthesis precursor (e.g. MI)
        "octopamine": "o",
        "tyramine": "t",
        "betaine": "b",
    }
    cell_codes: dict[str, set[str]] = {}
    for row in csv.DictReader(WANG_CSV.open()):
        s = row["neurotransmitters"].lower()
        cell_codes[row["neuron"]] = {
            v for k, v in code.items() if re.search(r"\b" + re.escape(k) + r"\b", s)
        }
    out: dict[str, set[str]] = defaultdict(set)
    for cls, cc in by_class.items():
        for c in cc:
            out[cls] |= cell_codes.get(c.name, set())
    return out


def fmt(codes: set[str]) -> str:
    return ",".join(NAME[c] for c in sorted(codes)) or "-"


def main() -> None:
    nj, by_class = load_neurons_json()
    wa = load_wormatlas()
    wang = load_wang(by_class)

    rows = []
    for cls, cc in sorted(by_class.items()):
        if not is_neuron(cc):
            continue
        n = nj.get(cls, set())
        w = wa.get(cls, set()) | wa.get(norm_class(cls), set())
        g = wang.get(cls, set())
        agree = (n == w or not w) and (n == g or not g)
        rows.append(
            {
                "class": cls,
                "neurons_json": fmt(n),
                "wormatlas": fmt(w),
                "wang": fmt(g),
                "agree": "yes" if agree else "no",
            }
        )

    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["class", "neurons_json", "wormatlas", "wang", "agree"]
        )
        writer.writeheader()
        writer.writerows(rows)

    disc = [r for r in rows if r["agree"] == "no"]
    print(f"classes compared: {len(rows)}  |  disagreements: {len(disc)}")
    print(f"{'class':6} {'neurons.json':20} {'wormatlas':20} {'wang':20}")
    for r in disc:
        print(f"{r['class']:6} {r['neurons_json']:20} {r['wormatlas']:20} {r['wang']:20}")
    print(f"\nwrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
