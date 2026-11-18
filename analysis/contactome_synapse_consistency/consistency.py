"""Contactome QC: do CIRCE's synapses fall within the membrane contactome, same-volume?

A chemical synapse or gap junction is physically impossible without a membrane apposition, so WITHIN a
single EM reconstruction every synaptic pair must also be a contact pair (synapses are a subset of
contacts). CIRCE now holds three contactomes, each scored on the SAME EM volume as a synaptic dataset
it already has, which lets us run that identity as a consistency check in BOTH directions:

  * synapse WITHOUT contact  -> a cross-scoring inconsistency / candidate error (should be ~empty).
  * contact WITHOUT synapse  -> a "potential synapse" (cells that touch but do not synapse; Peters'
                                rule / Stepanyants-Chklovskii). Reported as a count only, for context.

The same-volume pairings (synaptic dataset <-> contactome dataset, same animal/reconstruction):
  white_1986_jsh  <-> brittin_2021_l4_contactome      (L4, JSH series)
  white_1986_n2u  <-> brittin_2021_adult_contactome   (adult, N2U series)
  yim_2024_dauer  <-> yim_2024_dauer_contactome       (dauer; Yim scored both layers)

Prior art (reconciled per ../../CLAUDE.md): Brittin et al. 2021 (Nature 591:105, PMID 33627874,
PMC11648602) state the premise ("the space of possible synaptic and gap junction contact is restricted
by the existence of a physical membrane contact") but do NOT report this QC. Their validation is
contact-vs-contact (their algorithm vs White's contact neighborhoods for AIAR/AIBR/AQR -- where contact
data flagged a White *contact* mis-score), and the only synapse fraction they give is on *reproducible*
M4 contacts ("77% and 85% of synaptic and gap junction contacts occur at M4 contacts"), NOT the
all-contact fraction computed here. So this is the reverse (synapse->contact), all-contact, cross-
scoring consistency direction, which they do not report.

Deterministic (sorted; no RNG; independent of PYTHONHASHSEED). Reads the built graph's JSON projection
(run ``uv run cckg build`` first); the committed results.json captures the numbers. Run from a safe cwd:

    uv run python -I analysis/contactome_synapse_consistency/consistency.py
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CONN = REPO / "outputs" / "connectome.json"
OUT = Path(__file__).parent / "results.json"

#: (synaptic dataset, contactome dataset, human label) -- same EM volume on each row.
PAIRINGS = [
    ("white_1986_jsh", "brittin_2021_l4_contactome", "L4 (White JSH / Brittin L4)"),
    ("white_1986_n2u", "brittin_2021_adult_contactome", "adult (White N2U / Brittin adult)"),
    ("yim_2024_dauer", "yim_2024_dauer_contactome", "dauer (Yim synapses / Yim contactome)"),
]
HEAVY = (
    5  # weight >= HEAVY sections: a heavy miss is a real contact gap / naming mismatch to inspect
)


def _strip(s: str) -> str:
    return str(s).split("/")[-1]


def load_graph() -> dict:
    if not CONN.exists():
        raise SystemExit(f"{CONN} not found -- run `uv run cckg build` first.")
    return json.load(open(CONN))


def index(graph: dict):
    """-> (contact_pairs[ds], contact_cells[ds], syn_weight[ds][pair], syn_types[ds][pair])."""
    contact_pairs: dict[str, set] = defaultdict(set)
    contact_cells: dict[str, set] = defaultdict(set)
    syn_weight: dict[str, dict] = defaultdict(lambda: defaultdict(float))
    syn_types: dict[str, dict] = defaultdict(lambda: defaultdict(set))
    for c in graph["connections"]:
        ds, ct = _strip(c["dataset"]), _strip(c["connection_type"])
        a, b = _strip(c["pre"]), _strip(c["post"])
        if a == b:
            continue
        pair = frozenset((a, b))
        if ct == "adjacency":
            contact_pairs[ds].add(pair)
            contact_cells[ds].update((a, b))
        elif ct in ("chemical", "gap_junction"):
            syn_weight[ds][pair] += float(c["weight"])
            syn_types[ds][pair].add("gap_junction" if ct == "gap_junction" else "chemical")
    return contact_pairs, contact_cells, syn_weight, syn_types


def analyse(graph: dict) -> dict:
    contact_pairs, contact_cells, syn_weight, syn_types = index(graph)
    cases = []
    for sds, cds, label in PAIRINGS:
        cpairs, cells = contact_pairs[cds], contact_cells[cds]
        # Restrict synapse pairs to cells the contactome actually reconstructed (both endpoints), so a
        # miss means "touch not scored", never "cell absent from the contactome's volume". NOTE this is
        # necessary but NOT sufficient for White<->Brittin: the contactomes are nerve-ring-only but the
        # White synapse scoring extends beyond it, so a miss can also be an out-of-volume synapse. That
        # scope confound makes the White<->Brittin fractions LOWER BOUNDS; which cause applies to any one
        # miss is not resolvable here. The dauer case has no scope confound (Yim scored synapses AND
        # contacts from the same reconstruction), so it is the one clean subset estimate.
        spairs = {p for p in syn_weight[sds] if len(p) == 2 and set(p) <= cells}
        # Context only: is the miss pair scored as a contact in ANOTHER contactome (different animal AND
        # stage)? This does NOT establish the cause of the miss -- contact in another stage does not
        # transfer to this volume (neurites remodel), and "apposed elsewhere" is equally consistent with
        # a contact-scoring gap as with an out-of-volume synapse. It only shows the cells are not
        # categorically non-adjacent. The cause of an individual miss is not resolvable here (CIRCE has
        # no per-synapse coordinates); only the dauer same-scope control gives a clean subset estimate.
        elsewhere = {p for other, ps in contact_pairs.items() if other != cds for p in ps}
        misses = sorted(
            (
                {
                    "pair": sorted(p),
                    "weight": round(syn_weight[sds][p], 3),
                    "types": sorted(syn_types[sds][p]),
                    "contact_in_other_contactome": p in elsewhere,
                }
                for p in spairs - cpairs
            ),
            key=lambda m: (-m["weight"], m["pair"]),
        )
        at_contact = len(spairs) - len(misses)
        # contact-without-synapse ("potential synapses"): contact pairs whose cells both have synapses
        # in the paired dataset but no synapse between them. Cells restricted likewise, for symmetry.
        syn_cells = {c for p in syn_weight[sds] for c in p}
        potential = sum(
            1 for p in cpairs if p not in syn_weight[sds] and set(p) <= (cells & syn_cells)
        )
        cases.append(
            {
                "label": label,
                "synapse_dataset": sds,
                "contactome_dataset": cds,
                "synapse_pairs_in_contactome_cells": len(spairs),
                "contact_pairs": len(cpairs),
                "synapses_at_a_contact": at_contact,
                "fraction_at_contact": round(at_contact / len(spairs), 4) if spairs else None,
                "synapses_without_contact": len(misses),
                "misses_weight_1": sum(1 for m in misses if m["weight"] <= 1),
                "misses_contact_in_other_contactome": sum(
                    1 for m in misses if m["contact_in_other_contactome"]
                ),
                "misses_heavy": [m for m in misses if m["weight"] >= HEAVY],
                "all_misses": misses,
                "potential_synapses_contact_without_synapse": potential,
            }
        )
    return {
        "description": "Same-volume synapse<->contact consistency (synapses must occur at a contact).",
        "heavy_threshold_sections": HEAVY,
        "cases": cases,
    }


def main() -> None:
    results = analyse(load_graph())
    OUT.write_text(json.dumps(results, indent=2) + "\n")
    for c in results["cases"]:
        print(f"\n=== {c['label']} ===")
        print(
            f"  synapses at a contact: {c['synapses_at_a_contact']}/"
            f"{c['synapse_pairs_in_contactome_cells']} = {c['fraction_at_contact']:.1%}"
        )
        print(
            f"  synapse WITHOUT contact: {c['synapses_without_contact']} "
            f"(weight=1: {c['misses_weight_1']}; heavy >= {HEAVY}: {len(c['misses_heavy'])}; "
            f"apposed in another contactome: {c['misses_contact_in_other_contactome']})"
        )
        for m in c["misses_heavy"]:
            tag = " [apposed in another stage]" if m["contact_in_other_contactome"] else ""
            print(
                f"      {m['pair'][0]:8} {m['pair'][1]:8} weight={m['weight']:>5g} "
                f"{','.join(m['types'])}{tag}"
            )
        print(
            f"  potential synapses (contact without synapse): "
            f"{c['potential_synapses_contact_without_synapse']}"
        )
    print(f"\nwrote {OUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
