"""Threshold sensitivity of the Part C predicted chemical-synapse sign (increment 4 follow-on).

The biological reality check (``validate.py`` / README Result 2-3) found the one real limitation
of the sign layer: **glutamatergic synapses are over-called inhibitory**, driven by the broadly
expressed glutamate-gated chloride channels ``avr-14`` and ``glc-4``. The layer uses CeNGEN's
**threshold 4** matrix -- CeNGEN's *most stringent* (fewest genes called expressed). This script
asks the obvious follow-up: **is the over-call a tunable artifact that a different CeNGEN threshold
tames, or a genuine data-level limitation?**

CeNGEN publishes four binary matrices (1 liberal -> 4 stringent; thresholds are nested:
expr(t4) subset expr(t3) subset expr(t2) subset expr(t1)). Threshold 4 is already the floor, so the
only directions available are *less* stringent. We reclassify the same in-scope Cook-2019
hermaphrodite chemical connections at each threshold (swapping ONLY the postsynaptic receptor
expression; the NT identity, pre/post pairs and receptor polarity table are threshold-independent)
and measure, per threshold: the GluCl breadth, the glutamate exc/inh/complex split, the raw and
receptor-ratio-resolved E:I, and the functional-sign agreement vs Randi 2023.

The four CeNGEN matrices are NOT vendored (~10-14 MB each). Threshold 4 is read from the committed
projection (``ionotropic_receptor_expression.csv``); thresholds 1-3 are projected here from the raw
matrices, whose directory is given via ``$CENGEN_MATRICES_DIR`` or argv[1] (expects files named
``021821_liberal_threshold1.csv`` / ``021821_medium_threshold2.csv`` /
``021821_conservative_threshold3.csv``; sha256s in ``threshold_sensitivity_results.json``).

Deterministic (sorted, no RNG). Run from the repo root after ``uv run cckg build``:

    CENGEN_MATRICES_DIR=/path/to/matrices uv run python -I \
        analysis/partc_sign_validation/threshold_sensitivity.py
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

from celegans_connectome_kg.ingest.synapse_sign import (
    classify_sign,
    decode_fast_transmitters,
    read_receptor_catalog,
    read_receptor_expression,
)

REPO = Path(__file__).resolve().parents[2]
CONNECTOME = REPO / "outputs" / "connectome.json"
DATA = REPO / "data" / "cengen-ionotropic-receptors"
NEURONS = DATA.parent / "neuron-graph" / "neurons.json"
WANG = REPO / "data" / "wang-neurotransmitter-atlas" / "sex_neurotransmitters.csv"
OUT = Path(__file__).parent / "threshold_sensitivity_results.json"
PRIMARY = "cook_2019_hermaphrodite"
SIGN_DS = "fenyves_derived_synapse_sign"

# the six glutamate-gated chloride channels (anion glutamate receptors); avr-14/glc-4 are the
# over-call drivers the reality check flagged.
GLUCL = ["avr-14", "avr-15", "glc-1", "glc-2", "glc-3", "glc-4"]

MATRIX_FILES = {
    1: "021821_liberal_threshold1.csv",
    2: "021821_medium_threshold2.csv",
    3: "021821_conservative_threshold3.csv",
}

# Import the committed projection's class->cell map + ">0" test so the 19 aggregate class mappings
# never drift between the shipped build and this analysis.
_spec = importlib.util.spec_from_file_location("_cengen_build", DATA / "build_expression.py")
_bx = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_bx)
SPECIAL, _pos = _bx.SPECIAL, _bx._pos


def _strip(s: str) -> str:
    return str(s).split("/")[-1]


def _ratio(e: int, i: int) -> float | None:
    return round(e / i, 2) if i else None


def project(matrix_path: Path, wanted_bare: set[str]) -> dict[str, set[str]]:
    """Project a CeNGEN class matrix onto CIRCE cells -> {cell: {"WB:WBGene..."}}.

    Mirrors ``build_expression.py`` exactly (same SPECIAL map, same ">0" test, same conservative
    union collapses), so threshold 4 run through this reproduces the committed expression CSV.
    """
    cells = json.loads(NEURONS.read_text())
    names = {c["name"] for c in cells}
    by_class: dict[str, list[str]] = defaultdict(list)
    for c in cells:
        by_class[c["classes"]].append(c["name"])

    def resolve(cl: str) -> list[str]:
        if cl in SPECIAL:
            members = SPECIAL[cl]
        elif cl in by_class:
            members = by_class[cl]
        elif cl in names:
            members = [cl]
        else:
            return []
        return [m for m in members if m in names]

    expr_by_wb: dict[str, set[str]] = {}
    with open(matrix_path, newline="") as f:
        rd = csv.reader(f)
        classes = [h.strip().strip('"') for h in next(rd)[3:]]
        for row in rd:
            wb = row[2].strip().strip('"')
            if wb in wanted_bare:
                expr_by_wb[wb] = {classes[i] for i, v in enumerate(row[3:]) if _pos(v)}

    out: dict[str, set[str]] = defaultdict(set)
    for wb, cls in expr_by_wb.items():
        for cl in cls:
            for cell in resolve(cl):
                out[cell].add(f"WB:{wb}")
    return out


def main() -> None:
    mdir = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("CENGEN_MATRICES_DIR")
    if not mdir:
        sys.exit(
            "set CENGEN_MATRICES_DIR (or argv[1]) to the dir holding the threshold 1-3 matrices"
        )
    mdir = Path(mdir)

    conn = json.loads(CONNECTOME.read_text())
    catalog = read_receptor_catalog(DATA / "ionotropic_receptor_genes.csv")
    wb_symbol = {wb: r.symbol for wb, r in catalog.items()}
    wanted_bare = {wb.removeprefix("WB:") for wb in catalog}
    glucl_wb = {r.wbgene for r in catalog.values() if r.symbol in GLUCL}

    # NT identity exactly as the build/validate uses it: neuron-graph (nt-curated) + Wang herm.
    nt_by_cell: dict[str, str] = {
        c["name"]: str(c.get("neurotransmitter") or "") for c in conn["cells"]
    }
    for line in WANG.read_text().splitlines()[1:]:
        cell, sex, nt, *_ = line.split(",")
        if sex == "hermaphrodite":
            nt_by_cell[cell] = nt

    # In-scope pairs (threshold-independent): chemical, primary connectome, pre has a fast NT.
    in_scope: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for x in conn["connections"]:
        if x["connection_type"] != "chemical" or _strip(x["dataset"]) != PRIMARY:
            continue
        pre, post = _strip(x["pre"]), _strip(x["post"])
        if (pre, post) in seen:
            continue
        if decode_fast_transmitters(nt_by_cell.get(pre)):
            seen.add((pre, post))
            in_scope.append((pre, post))
    in_scope.sort()

    # measured functional sign (Randi 2023) for the reality check.
    func_sign: dict[tuple[str, str], str] = {}
    for x in conn["connections"]:
        if x["connection_type"] == "functional" and x.get("sign") in ("excitatory", "inhibitory"):
            func_sign[(_strip(x["pre"]), _strip(x["post"]))] = x["sign"]

    # expression per threshold: t4 = committed CSV (canonical shipped build); t1-3 = projected here.
    expr_by_t: dict[int, dict[str, set[str]]] = {}
    sha_by_t: dict[int, str] = {}
    expr_by_t[4] = read_receptor_expression(DATA / "ionotropic_receptor_expression.csv")
    for t, fname in MATRIX_FILES.items():
        p = mdir / fname
        if not p.exists():
            sys.exit(f"missing matrix: {p}")
        sha_by_t[t] = hashlib.sha256(p.read_bytes()).hexdigest()
        expr_by_t[t] = project(p, wanted_bare)

    def ion_counts(expr: dict[str, set[str]], pre: str, post: str) -> Counter:
        tx = decode_fast_transmitters(nt_by_cell.get(pre))
        c: Counter = Counter()
        for wb in expr.get(post, ()):
            r = catalog.get(wb)
            if r and not r.excluded and r.transmitter in tx:
                c[r.ion] += 1
        return c

    def resolve_ratio(cation: int, anion: int) -> str:
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
        return "unknown"

    def analyse(t: int) -> dict:
        expr = expr_by_t[t]
        raw_outcome: dict[tuple[str, str], str] = {}
        raw_sign: dict[tuple[str, str], str | None] = {}
        resolved: dict[tuple[str, str], str] = {}
        glu_breakdown: Counter = Counter()
        for pre, post in in_scope:
            tx = decode_fast_transmitters(nt_by_cell.get(pre))
            sign, outcome = classify_sign(tx, expr.get(post, set()), catalog)
            raw_outcome[(pre, post)] = outcome
            raw_sign[(pre, post)] = sign
            ic = ion_counts(expr, pre, post)
            resolved[(pre, post)] = resolve_ratio(ic["cation"], ic["anion"])
            if "glutamate" in tx:
                glu_breakdown[outcome] += 1

        raw = Counter(raw_outcome.values())
        res = Counter(resolved.values())
        n = len(in_scope)

        # GluCl breadth: cells (and CeNGEN classes via distinct expression) calling each anion GluCl.
        glucl_cells = {g: 0 for g in GLUCL}
        for cell, genes in expr.items():
            for wb in genes & glucl_wb:
                glucl_cells[wb_symbol[wb]] += 1

        # reality check vs measured functional sign.
        def concordance(pred: dict[tuple[str, str], str]) -> dict:
            pairs = [
                p for p, s in pred.items() if s in ("excitatory", "inhibitory") and p in func_sign
            ]
            conf = Counter((pred[p], func_sign[p]) for p in pairs)
            agree = sum(v for (a, b), v in conf.items() if a == b)
            return {
                "overlap": len(pairs),
                "agreement_pct": round(100 * agree / len(pairs), 1) if pairs else None,
                "confusion": {f"{a}|{b}": v for (a, b), v in sorted(conf.items())},
            }

        # which anion genes drive inhibitory-predicted / measured-excitatory (the over-call).
        overcall: Counter = Counter()
        for p, s in raw_sign.items():
            if s == "inhibitory" and func_sign.get(p) == "excitatory":
                tx = decode_fast_transmitters(nt_by_cell.get(p[0]))
                for wb in expr.get(p[1], ()):
                    r = catalog.get(wb)
                    if r and not r.excluded and r.ion == "anion" and r.transmitter in tx:
                        overcall[wb_symbol[wb]] += 1

        definite = raw["excitatory"] + raw["inhibitory"]
        return {
            "raw_counts": dict(raw),
            "definite_pct": round(100 * definite / n, 1),
            "complex_pct": round(100 * raw["complex"] / n, 1),
            "raw_E:I": _ratio(raw["excitatory"], raw["inhibitory"]),
            "resolved_counts": dict(res),
            "resolved_E:I": _ratio(res["excitatory"], res["inhibitory"]),
            "glutamate_outcome": dict(glu_breakdown),
            "glutamate_E:I": _ratio(glu_breakdown["excitatory"], glu_breakdown["inhibitory"]),
            "glucl_cells_expressing": glucl_cells,
            "n_cells_with_any_receptor": len(expr),
            "reality_raw": concordance(raw_sign),
            "reality_resolved": concordance(resolved),
            "inhibitory_overcall_drivers": dict(overcall.most_common()),
            "_raw_outcome": raw_outcome,  # internal, popped before serialise
        }

    per_t = {t: analyse(t) for t in (1, 2, 3, 4)}

    # ----- self-consistency guards -----
    # (a) threshold 4 here must reproduce the committed validate.py raw counts.
    t4 = per_t[4]["raw_counts"]
    expected_t4 = {"complex": 1842, "excitatory": 704, "inhibitory": 465, "unpredicted": 1283}
    t4_ok = t4 == expected_t4

    # (b) nesting: expr(t4) subset expr(t3) subset expr(t2) subset expr(t1) at the (cell,gene) level.
    def flat(expr: dict[str, set[str]]) -> set[tuple[str, str]]:
        return {(c, g) for c, gs in expr.items() for g in gs}

    sets = {t: flat(expr_by_t[t]) for t in (1, 2, 3, 4)}
    nesting_ok = sets[4] <= sets[3] <= sets[2] <= sets[1]

    for t in (1, 2, 3, 4):
        per_t[t].pop("_raw_outcome")  # internal scratch, not serialised

    results = {
        "primary_connectome": PRIMARY,
        "n_in_scope_chemical_connections": len(in_scope),
        "matrix_sha256": {f"threshold{t}": sha_by_t[t] for t in (1, 2, 3)},
        "self_checks": {
            "threshold4_reproduces_committed_raw_counts": t4_ok,
            "expression_nested_t4_subset_t3_subset_t2_subset_t1": nesting_ok,
            "expression_edge_counts": {f"threshold{t}": len(sets[t]) for t in (1, 2, 3, 4)},
        },
        "per_threshold": {f"threshold{t}": per_t[t] for t in (1, 2, 3, 4)},
    }
    OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    # ----- console summary -----
    print(f"in-scope chemical connections: {len(in_scope)}")
    print(f"self-check t4==committed: {t4_ok}   nesting t4<=t3<=t2<=t1: {nesting_ok}")
    print(f"expression edges by threshold: {results['self_checks']['expression_edge_counts']}")
    print(
        "\nthr | raw exc/inh/cpx/unp            | def%  cpx%  rawE:I resE:I | glu exc/inh/cpx  gluE:I"
    )
    for t in (1, 2, 3, 4):
        d = per_t[t]
        rc = d["raw_counts"]
        g = d["glutamate_outcome"]
        tag = " (stringent, shipped)" if t == 4 else ""
        print(
            f" {t}  | {rc.get('excitatory', 0):4}/{rc.get('inhibitory', 0):4}/"
            f"{rc.get('complex', 0):4}/{rc.get('unpredicted', 0):4} "
            f"| {d['definite_pct']:4} {d['complex_pct']:4} {str(d['raw_E:I']):>5} {str(d['resolved_E:I']):>5}"
            f"  | {g.get('excitatory', 0):3}/{g.get('inhibitory', 0):3}/{g.get('complex', 0):4}"
            f"  {str(d['glutamate_E:I']):>5}{tag}"
        )
    print("\nGluCl breadth (cells calling each anion glutamate receptor):")
    for t in (1, 2, 3, 4):
        print(f"  t{t}: {per_t[t]['glucl_cells_expressing']}")
    print("\nreality check vs measured functional sign (raw / resolved agreement):")
    for t in (1, 2, 3, 4):
        rr, rs = per_t[t]["reality_raw"], per_t[t]["reality_resolved"]
        print(
            f"  t{t}: raw {rr['agreement_pct']}% (n={rr['overlap']})  "
            f"resolved {rs['agreement_pct']}% (n={rs['overlap']})"
        )
    print("\ninhibitory-overcall anion drivers (predicted inh / measured exc):")
    for t in (1, 2, 3, 4):
        print(f"  t{t}: {per_t[t]['inhibitory_overcall_drivers']}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
