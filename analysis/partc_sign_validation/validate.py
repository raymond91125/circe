"""Validate CIRCE's Part C predicted chemical-synapse sign (increment 4).

Three checks over the built connectome (``outputs/connectome.json``), all deterministic (sorted,
no RNG), primary target ``cook_2019_hermaphrodite`` (= WormWiring, Fenyves' own connectome):

  1. HEADLINE STATS vs the two verified Fenyves references, counted the comparable way.
     - RAW: the sign as stored in the KG (complex left unset).
     - RECEPTOR-RATIO-RESOLVED: complex relabelled excitatory/inhibitory by the majority of
       matched cation vs anion receptor genes (ties -> unknown) -- replicating Fenyves' own group
       (Szilagyi et al., BMC Biol 2026;24:167). This is an ANALYSIS view only; the KG stays raw.

  2. BIOLOGICAL REALITY CHECK vs measured functional sign (Randi 2023, Part B). For (pre,post)
     pairs that have BOTH a predicted chemical sign and a measured functional sign, how often do
     they agree, and what is the confusion? Functional sign is a NET signal-propagation effect
     (possibly polysynaptic / extrasynaptic), not a monosynaptic ground truth -- so this bounds,
     not proves, correctness. Reported honestly.

  3. NT-IDENTITY SANITY + the inhibitory-overcall diagnostic: predicted sign by presynaptic
     transmitter, and which anion receptor genes drive the functional-overlap disagreements.

Reads only committed inputs + the build output; writes results.json next to this file. Run from the
repo root after ``uv run cckg build``:

    uv run python -I analysis/partc_sign_validation/validate.py
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

from celegans_connectome_kg.ingest.synapse_sign import (
    decode_fast_transmitters,
    read_receptor_catalog,
    read_receptor_expression,
)

REPO = Path(__file__).resolve().parents[2]
CONNECTOME = REPO / "outputs" / "connectome.json"
DATA = REPO / "data" / "cengen-ionotropic-receptors"
WANG = REPO / "data" / "wang-neurotransmitter-atlas" / "sex_neurotransmitters.csv"
OUT = Path(__file__).parent / "results.json"
PRIMARY = "cook_2019_hermaphrodite"
SIGN_DS = "fenyves_derived_synapse_sign"

# Verified reference numbers (see README for sources/quotes).
FENYVES_2020 = {  # per SYNAPSE, neuron->neuron, 20,589 synapses
    "excitatory": 9070,
    "inhibitory": 2580,
    "complex": 3413,
    "unpredicted": 5526,
    "E:I": 3.52,
    "predicted_pct": 73.2,
}
SZILAGYI_2026 = {  # per CONNECTION, 3,638; complex resolved by receptor ratio (ties -> unknown)
    "excitatory": 1555,
    "inhibitory": 553,
    "unknown": 1530,
    "E:I": 2.81,
}


def _strip(s: str) -> str:
    return str(s).split("/")[-1]


def _ratio(e: int, i: int) -> float | None:
    return round(e / i, 2) if i else None


def main() -> None:
    conn = json.loads(CONNECTOME.read_text())
    catalog = read_receptor_catalog(DATA / "ionotropic_receptor_genes.csv")
    expr = read_receptor_expression(DATA / "ionotropic_receptor_expression.csv")

    # Presynaptic NT identity exactly as the build uses it: neuron-graph (nt-curated) + Wang herm.
    nt_by_cell: dict[str, str] = {
        c["name"]: str(c.get("neurotransmitter") or "") for c in conn["cells"]
    }
    for line in WANG.read_text().splitlines()[1:]:
        cell, sex, nt, *_ = line.split(",")
        if sex == "hermaphrodite":
            nt_by_cell[cell] = nt

    # ----- gather, per connection, the stored (raw) predicted outcome + a recomputed ion tally -----
    def ion_counts(pre: str, post: str) -> Counter:
        """Matched cation/anion receptor genes for pre's transmitter(s) that post expresses."""
        tx = decode_fast_transmitters(nt_by_cell.get(pre))
        c: Counter = Counter()
        for wb in expr.get(post, ()):
            r = catalog.get(wb)
            if r and not r.excluded and r.transmitter in tx:
                c[r.ion] += 1
        return c

    def resolve(cation: int, anion: int) -> str:
        if cation and not anion:
            return "excitatory"
        if anion and not cation:
            return "inhibitory"
        if not cation and not anion:
            return "unpredicted"
        if cation > anion:
            return "excitatory"
        if anion > cation:
            return "inhibitory"
        return "unknown"  # exact tie among mixed receptors

    # per (pre,post) on the primary connectome: raw outcome, raw stored sign, resolved label, ions
    raw_outcome: dict[tuple[str, str], str] = {}
    raw_sign: dict[tuple[str, str], str | None] = {}
    resolved: dict[tuple[str, str], str] = {}
    ions: dict[tuple[str, str], Counter] = {}
    mismatch_selfcheck = 0
    for x in conn["connections"]:
        if x["connection_type"] != "chemical" or _strip(x["dataset"]) != PRIMARY:
            continue
        outcome = next(
            (e["value"] for e in x.get("evidence", []) if _strip(e["dataset"]) == SIGN_DS), None
        )
        if outcome is None:
            continue
        pre, post = _strip(x["pre"]), _strip(x["post"])
        raw_outcome[(pre, post)] = outcome
        raw_sign[(pre, post)] = x.get("sign")
        ic = ion_counts(pre, post)
        ions[(pre, post)] = ic
        resolved[(pre, post)] = resolve(ic["cation"], ic["anion"])
        # self-consistency: recomputed raw label must match the KG's stored outcome
        recomputed = (
            "excitatory"
            if ic["cation"] and not ic["anion"]
            else "inhibitory"
            if ic["anion"] and not ic["cation"]
            else "complex"
            if ic["cation"] and ic["anion"]
            else "unpredicted"
        )
        if recomputed != outcome:
            mismatch_selfcheck += 1

    n = len(raw_outcome)
    raw_tally = Counter(raw_outcome.values())
    resolved_tally = Counter(resolved.values())

    def coverage_block(tally: Counter, categories_definite=("excitatory", "inhibitory")) -> dict:
        definite = sum(tally[c] for c in categories_definite)
        predicted_incl_complex = definite + tally.get("complex", 0)
        return {
            "counts": dict(tally),
            "n_in_scope": n,
            "definite_pct": round(100 * definite / n, 1) if n else None,
            "predicted_incl_complex_pct": round(100 * predicted_incl_complex / n, 1) if n else None,
            "E:I": _ratio(tally["excitatory"], tally["inhibitory"]),
        }

    # ----- biological reality check vs measured functional sign (Randi 2023) -----
    func_sign: dict[tuple[str, str], str] = {}
    for x in conn["connections"]:
        if x["connection_type"] == "functional" and x.get("sign") in ("excitatory", "inhibitory"):
            func_sign[(_strip(x["pre"]), _strip(x["post"]))] = x["sign"]

    def concordance(pred_label: dict[tuple[str, str], str]) -> dict:
        pairs = sorted(
            p for p, s in pred_label.items() if s in ("excitatory", "inhibitory") and p in func_sign
        )
        confusion = Counter((pred_label[p], func_sign[p]) for p in pairs)
        agree = sum(v for (a, b), v in confusion.items() if a == b)
        return {
            "overlap": len(pairs),
            "agree": agree,
            "agreement_pct": round(100 * agree / len(pairs), 1) if pairs else None,
            "confusion_pred_x_measured": {f"{a}|{b}": v for (a, b), v in sorted(confusion.items())},
        }

    # diagnostic: for inhibitory-predicted but measured-excitatory pairs, which anion genes fire
    wb_symbol = {wb: r.symbol for wb, r in catalog.items()}
    overcall_driver: Counter = Counter()
    for p in sorted(raw_sign):
        if raw_sign[p] == "inhibitory" and func_sign.get(p) == "excitatory":
            pre, post = p
            tx = decode_fast_transmitters(nt_by_cell.get(pre))
            for wb in expr.get(post, ()):
                r = catalog.get(wb)
                if r and not r.excluded and r.ion == "anion" and r.transmitter in tx:
                    overcall_driver[wb_symbol[wb]] += 1

    # ----- NT-identity sanity: predicted raw outcome by presynaptic transmitter -----
    nt_breakdown: dict[str, Counter] = defaultdict(Counter)
    for (pre, post), outcome in raw_outcome.items():
        for tx in sorted(decode_fast_transmitters(nt_by_cell.get(pre))):
            nt_breakdown[tx][outcome] += 1

    results = {
        "primary_connectome": PRIMARY,
        "self_check_recomputed_vs_stored_mismatches": mismatch_selfcheck,
        "headline": {
            "raw": coverage_block(raw_tally),
            "receptor_ratio_resolved": coverage_block(resolved_tally),
            "reference_fenyves_2020_per_synapse": FENYVES_2020,
            "reference_szilagyi_2026_per_connection": SZILAGYI_2026,
        },
        "reality_check_vs_measured_functional": {
            "note": (
                "Functional sign (Randi 2023) is a net signal-propagation effect, not a "
                "monosynaptic ground truth; this bounds, not proves, correctness."
            ),
            "n_measured_functional_signed_pairs": len(func_sign),
            "raw": concordance(raw_sign),
            "receptor_ratio_resolved": concordance(resolved),
            "inhibitory_overcall_anion_drivers": dict(overcall_driver.most_common()),
        },
        "nt_identity_sanity": {tx: dict(c) for tx, c in sorted(nt_breakdown.items())},
    }
    OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    # ----- console summary -----
    print(f"primary: {PRIMARY}   in-scope chemical connections: {n}")
    print(f"self-check (recomputed==stored): {n - mismatch_selfcheck}/{n} OK")
    r, rr = results["headline"]["raw"], results["headline"]["receptor_ratio_resolved"]
    print(f"\nRAW      : {r['counts']}")
    print(
        f"           definite={r['definite_pct']}%  incl-complex={r['predicted_incl_complex_pct']}%  E:I={r['E:I']}"
    )
    print(f"RESOLVED : {rr['counts']}")
    print(f"           definite={rr['definite_pct']}%  E:I={rr['E:I']}")
    print(
        f"ref Fenyves2020 (synapses): E:I={FENYVES_2020['E:I']} predicted={FENYVES_2020['predicted_pct']}%"
    )
    print(f"ref Szilagyi2026 (conns)  : E:I={SZILAGYI_2026['E:I']}")
    rc = results["reality_check_vs_measured_functional"]
    print("\nreality check vs measured functional sign:")
    print(
        f"  RAW      overlap={rc['raw']['overlap']} agree={rc['raw']['agreement_pct']}%  {rc['raw']['confusion_pred_x_measured']}"
    )
    print(
        f"  RESOLVED overlap={rc['receptor_ratio_resolved']['overlap']} agree={rc['receptor_ratio_resolved']['agreement_pct']}%  {rc['receptor_ratio_resolved']['confusion_pred_x_measured']}"
    )
    print(f"  inhibitory-overcall anion drivers: {rc['inhibitory_overcall_anion_drivers']}")
    print("\nNT-identity sanity (raw outcome by presynaptic transmitter):")
    for tx, c in results["nt_identity_sanity"].items():
        print(f"  {tx:14} {c}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
