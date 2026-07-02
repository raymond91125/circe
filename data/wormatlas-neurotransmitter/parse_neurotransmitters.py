#!/usr/bin/env python
"""Parse the pinned WormAtlas neurotransmitter table into structured rows.

Source: neurotransmitterstable.html (crawled snapshot — see MANIFEST.md).
The page is one large table organized into neurotransmitter sections, each an evidence
table with columns Description / Gene / Detection method / Localization / References. The
neuron assignments live in the Localization column. This extracts those evidence rows
faithfully, tagged by neurotransmitter (sections are delimited by the "Summary List of X"
rows). It deliberately does NOT tokenize Localization into individual neurons — that text
carries footnotes, (m) male-specific and (h) hermaphrodite markers, ranges (DD1-6), and
non-neuronal cells, and is best read as WormAtlas wrote it.

Run:  uv run --with 'beautifulsoup4,lxml' python parse_neurotransmitters.py
Outputs: neurotransmitters.csv, neurotransmitters.json (next to this script).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from bs4 import BeautifulSoup

HERE = Path(__file__).resolve().parent
HTML = HERE / "neurotransmitterstable.html"

# "Summary List of X (Neurons)" rows delimit the sections → canonical neurotransmitter.
SUMMARY_TO_NT = [
    ("ACh", "acetylcholine"),
    ("Dopaminergic", "dopamine"),
    ("TA & OA", "tyramine/octopamine"),
    ("5HT", "serotonin"),
    ("GABA", "GABA"),
    ("Glu", "glutamate"),
]
HEADER_CELLS = {"DESCRIPTION", "GENE NAME", "DETECTION METHOD", "LOCALIZATION", "REFERENCES"}
FIELDS = [
    "neurotransmitter",
    "description",
    "gene",
    "detection_method",
    "localization",
    "references",
]


def _section(text: str) -> str | None:
    if "Summary List" not in text:
        return None
    for key, nt in SUMMARY_TO_NT:
        if key.lower() in text.lower():
            return nt
    return None


def parse(html_path: Path = HTML) -> list[dict]:
    soup = BeautifulSoup(html_path.read_text(encoding="utf-8", errors="replace"), "lxml")
    # The main data table is the one whose first section is ACETYLCHOLINE.
    table = max(soup.find_all("table"), key=lambda t: len(t.find_all("tr")))
    current: str | None = None
    records: list[dict] = []
    for row in table.find_all("tr"):
        cells = [c.get_text(" ", strip=True) for c in row.find_all(["td", "th"])]
        if not cells:
            continue
        section = _section(" ".join(cells))
        if section:
            current = section
            continue
        if len(cells) == 5 and cells[0].upper() not in HEADER_CELLS and current:
            desc, gene, method, loc, refs = cells
            records.append(
                {
                    "neurotransmitter": current,
                    "description": desc,
                    "gene": gene,
                    "detection_method": method,
                    "localization": loc,
                    "references": refs,
                }
            )
    return records


def main() -> None:
    records = parse()
    with (HERE / "neurotransmitters.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(records)
    (HERE / "neurotransmitters.json").write_text(json.dumps(records, indent=1))
    print(f"parsed {len(records)} evidence rows -> neurotransmitters.csv / .json")


if __name__ == "__main__":
    main()
