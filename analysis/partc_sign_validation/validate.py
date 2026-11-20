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

  4. EXTERNAL VALIDATION vs Fenyves 2020 S1 Data (optional; needs the file). Two parts:
     - LOGIC CHECK: run our sign rule on Fenyves' OWN receptor-match booleans (their S1 Data sheet
       5) and confirm we reproduce their per-connection label -- isolates the rule from the data.
     - DATA-VINTAGE DECOMPOSITION: compare our per-connection predictions to theirs on the shared
       connections, and attribute each disagreement to an NT-source difference vs an expression
       difference (ours = CeNGEN threshold-4; theirs = a curated WormBase/WormAtlas/CeNGEN matrix).
     S1 Data (xlsx, ~1 MB, linkgroup.hu; a validation reference, NOT vendored) is given via
     ``$FENYVES_S1`` or argv[1]; its sha256 + the numbers are captured in results.json.

Reads committed inputs + the build output (+ optional S1 Data); writes results.json next to this
file. Run from the repo root after ``uv run cckg build``:

    uv run python -I analysis/partc_sign_validation/validate.py                       # checks 1-3
    FENYVES_S1=/path/to/S1_Data.xlsx uv run python -I analysis/partc_sign_validation/validate.py
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
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
FENYVES_S1_URL = "https://linkgroup.hu/docs/S1_Data.xlsx"  # NOT vendored; validation reference only
# Fenyves per-connection label (S1 Data sheet 5 col 16) -> our outcome vocabulary.
_FEN_MAP = {"+": "excitatory", "-": "inhibitory", "complex": "complex", "no pred": "unpredicted"}
_FEN_NT = {"Glu": "glutamate", "ACh": "acetylcholine", "GABA": "GABA"}


def _strip(s: str) -> str:
    return str(s).split("/")[-1]


def _ratio(e: int, i: int) -> float | None:
    return round(e / i, 2) if i else None


def _rule_from_booleans(gp, gn, ap, an, bp, bn) -> str:
    """Our sign rule expressed over Fenyves' six receptor-match booleans (+/- x Glu/ACh/GABA),
    returning their vocabulary (+/-/complex/no pred). Used only for the logic check."""
    cation, anion = (gp or ap or bp), (gn or an or bn)
    if cation and not anion:
        return "+"
    if anion and not cation:
        return "-"
    if cation and anion:
        return "complex"
    return "no pred"


def fenyves_s1_comparison(
    path: Path,
    raw_outcome: dict[tuple[str, str], str],
    nt_by_cell: dict[str, str],
) -> dict:
    """Logic check + data-vintage decomposition vs Fenyves 2020 S1 Data (sheets 1 + 5)."""
    import openpyxl

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    # their NT per neuron (sheet 1): dominant + alternative -> fast-transmitter set
    their_fast: dict[str, frozenset[str]] = {}
    for r in wb["1. NT expr"].iter_rows(min_row=2, values_only=True):
        if r[0]:
            their_fast[r[0]] = frozenset(_FEN_NT[x] for x in (r[1], r[2]) if x in _FEN_NT)
    # sheet 5: per chemical connection, their label + the six receptor-match booleans
    fen_pred: dict[tuple[str, str], str] = {}
    logic_ok = logic_n = 0
    for r in wb["5. Sign prediction"].iter_rows(min_row=3, values_only=True):
        if r[5] != "chemical" or not r[0] or not r[3]:
            continue
        fen_pred[(r[0], r[3])] = r[16]
        booleans = [str(r[j]) == "True" for j in range(6, 12)]
        logic_n += 1
        logic_ok += _rule_from_booleans(*booleans) == r[16]

    # output comparison on connections both have AND we make a sign call for (jointly in-scope)
    common = sorted(p for p in raw_outcome if p in fen_pred)
    confusion = Counter((raw_outcome[p], _FEN_MAP[fen_pred[p]]) for p in common)
    agree = sum(v for (a, b), v in confusion.items() if a == b)
    both_definite = [
        p
        for p in common
        if raw_outcome[p] in ("excitatory", "inhibitory")
        and _FEN_MAP[fen_pred[p]] in ("excitatory", "inhibitory")
    ]
    bd_agree = sum(1 for p in both_definite if raw_outcome[p] == _FEN_MAP[fen_pred[p]])
    # attribute disagreements: NT(source) differs ours-vs-theirs, else expression-driven
    disagree = [p for p in common if raw_outcome[p] != _FEN_MAP[fen_pred[p]]]
    nt_driven = sum(
        1
        for p in disagree
        if decode_fast_transmitters(nt_by_cell.get(p[0])) != their_fast.get(p[0], frozenset())
    )
    return {
        "source_url": FENYVES_S1_URL,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "logic_check": {
            "note": "our rule on Fenyves' own receptor-match booleans (S1 sheet 5) vs their label",
            "reproduced": logic_ok,
            "total": logic_n,
            "pct": round(100 * logic_ok / logic_n, 2) if logic_n else None,
        },
        "fenyves_per_connection_label_dist": dict(Counter(fen_pred.values())),
        "output_comparison_raw": {
            "jointly_in_scope": len(common),
            "agreement_pct": round(100 * agree / len(common), 1) if common else None,
            "both_make_definite_call": len(both_definite),
            "agree_on_direction_when_both_definite_pct": (
                round(100 * bd_agree / len(both_definite), 1) if both_definite else None
            ),
            "confusion_ours_x_fenyves": {f"{a}|{b}": v for (a, b), v in sorted(confusion.items())},
        },
        "divergence_attribution": {
            "disagreements": len(disagree),
            "nt_source_driven": nt_driven,
            "expression_driven": len(disagree) - nt_driven,
        },
    }


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

    # ----- optional external validation vs Fenyves S1 Data -----
    s1_path = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("FENYVES_S1")
    if s1_path and Path(s1_path).exists():
        results["external_fenyves_s1"] = fenyves_s1_comparison(
            Path(s1_path), raw_outcome, nt_by_cell
        )

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
    if "external_fenyves_s1" in results:
        ex = results["external_fenyves_s1"]
        lc, oc, da = ex["logic_check"], ex["output_comparison_raw"], ex["divergence_attribution"]
        print("\nexternal validation vs Fenyves S1 Data:")
        print(
            f"  LOGIC CHECK (our rule on their data): {lc['reproduced']}/{lc['total']} = {lc['pct']}%"
        )
        print(
            f"  output agreement (raw): {oc['agreement_pct']}% over {oc['jointly_in_scope']} shared"
        )
        print(
            f"  agree on DIRECTION when both make a definite call: "
            f"{oc['agree_on_direction_when_both_definite_pct']}% of {oc['both_make_definite_call']}"
        )
        print(
            f"  disagreements: {da['disagreements']}  expression-driven={da['expression_driven']}  "
            f"NT-driven={da['nt_source_driven']}"
        )
    else:
        print("\n(external Fenyves S1 check skipped -- set FENYVES_S1=/path/to/S1_Data.xlsx)")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
