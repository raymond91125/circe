"""Parse the Wang et al. 2024 (eLife 95402) male neurotransmitter data into per-sex assignments.

Produces ``sex_neurotransmitters.csv`` (cell, sex, neurotransmitter, confidence, note), the build
input for the reified NeurotransmitterAssignment records. Two sources:

- **Supp File 3** (``Supp3_male_expression.xlsx``): the male reporter atlas. The consolidated
  "Neurotransmitter(s)" column (col 21) is mapped mechanically to the neuron-graph code scheme
  (a=ACh, l=Glu, g=GABA, d=dopamine, s=serotonin, u=unknown). "(new)" = newly identified; a
  leading "*" or a "?" marks a tentative call -> confidence "putative", else "reported". These are
  the male-specific neurons.

- **Supp File 4** (``Supp4_dimorphic.xlsx``): sex-shared neurons whose transmitter identity differs
  by sex. The male call is derived (curated) from the documented reporter change applied to the
  hermaphrodite base call (from the built KG); recorded below with a note. AVG is excluded: its
  difference is unc-17 expression level only (cholinergic in both sexes), not an identity change.

Run:  python data/wang-neurotransmitter-atlas/parse_male_atlas.py
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

import openpyxl

#: Supp3 groups serially-repeated neurons as "DX1/2", "EF3/4"; expand to the member cell names.
_GROUP = re.compile(r"^([A-Za-z]+)(\d+)/(\d+)$")


def _expand(name: str) -> list[str]:
    m = _GROUP.match(name)
    if m:
        prefix, a, b = m.group(1), m.group(2), m.group(3)
        return [f"{prefix}{a}", f"{prefix}{b}"]
    return [name]

HERE = Path(__file__).resolve().parent
OUT = HERE / "sex_neurotransmitters.csv"

# Curated sexually-dimorphic identity differences (Supp File 4), one entry per class. herm = the
# hermaphrodite code already in the KG (neuron-graph + reconciliation); male = herm modified by the
# Supp4 reporter change. Members expand the class to its KG cells.
DIMORPHIC = [
    # class, members,          herm,  male,   note
    ("ADF", ["ADFL", "ADFR"], "as", "ags", "male gains GABA (unc-47+) [Wang 2024 Supp4]"),
    ("AS10", ["AS10"], "a", "ag", "male gains GABA (unc-47+) [Wang 2024 Supp4]"),
    ("AS11", ["AS11"], "a", "ag", "male gains GABA (unc-47+) [Wang 2024 Supp4]"),
    ("PHC", ["PHCL", "PHCR"], "bl", "bgl", "male gains GABA (unc-47+); eat-4 stronger [Supp4]"),
    ("PDB", ["PDB"], "a", "ag", "male gains GABA (unc-47+) [Wang 2024 Supp4]"),
    ("PVN", ["PVNL", "PVNR"], "abl", "abgl", "male gains GABA (unc-47+); eat-4 stronger [Supp4]"),
    ("AIM", ["AIML", "AIMR"], "ls", "as", "Glu->ACh switch in male (eat-4 down, unc-17 up) [Supp4]"),
    ("PVW", ["PVWL", "PVWR"], "u", "s", "male gains serotonin (anti-5-HT+) [Wang 2024 Supp4]"),
]


def _code(call: str) -> tuple[str, str]:
    """(neurotransmitter code, confidence) from a Supp3 consolidated call."""
    c = call.strip()
    conf = "putative" if c.startswith("*") or "?" in c else "reported"
    low = c.lstrip("*").strip().lower()
    if low.startswith("ach"):
        code = "a"
    elif low.startswith("glu"):
        code = "l"
    elif low.startswith("gaba"):
        code = "g"
    elif low.startswith("da"):
        code = "d"
    elif low.startswith("5-ht"):
        code = "s"
    else:  # "unknown (orphan…)", "…unknown monoamine?" -> unknown
        code = "u"
    return code, conf


def main() -> None:
    rows: list[dict] = []

    # --- Supp3: male-specific neurons ---
    wb = openpyxl.load_workbook(HERE / "Supp3_male_expression.xlsx", read_only=True, data_only=True)
    ws = wb["Supp File 3"]
    for r in list(ws.iter_rows(values_only=True))[5:]:
        neuron, call = r[2], r[21]
        if not (neuron and call and str(neuron).strip()):
            continue
        code, conf = _code(str(call))
        for cell in _expand(str(neuron).strip()):
            rows.append(
                {
                    "cell": cell,
                    "sex": "male",
                    "neurotransmitter": code,
                    "confidence": conf,
                    "note": f"Wang 2024 Supp3 call: {str(call).strip()}",
                }
            )
    wb.close()

    # --- Supp4: sex-shared dimorphic neurons (both sexes recorded) ---
    for _cls, members, herm, male, note in DIMORPHIC:
        for cell in members:
            rows.append(
                {"cell": cell, "sex": "hermaphrodite", "neurotransmitter": herm, "confidence": "reported", "note": note}
            )
            rows.append(
                {"cell": cell, "sex": "male", "neurotransmitter": male, "confidence": "reported", "note": note}
            )

    # De-duplicate (cell, sex): the grouped serial neurons (DX1/2, EF3/4) are listed under more
    # than one section of Supp3, so their expansions repeat. The call is identical -> keep one.
    seen: dict[tuple[str, str], dict] = {}
    for r in rows:
        key = (r["cell"], r["sex"])
        if key in seen:
            if seen[key]["neurotransmitter"] != r["neurotransmitter"]:
                raise ValueError(f"conflicting NT for {key}: {seen[key]} vs {r}")
            continue
        seen[key] = r
    rows = list(seen.values())

    rows.sort(key=lambda x: (x["cell"], x["sex"]))
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cell", "sex", "neurotransmitter", "confidence", "note"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} assignments)")


if __name__ == "__main__":
    main()
