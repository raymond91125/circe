"""Per-synapse reconciliation of the Part C predicted chemical-synapse sign (increment 4 follow-on).

``validate.py`` / README Results 1 and 4 compare CIRCE's sign tally to Fenyves 2020 **per
connection**, while Fenyves' headline is **per synapse** (20,589 synapses). This script closes that
counting-unit gap directly: it weights each Cook-2019-hermaphrodite chemical connection by its
synapse count (the connectome ``weight`` = Cook's aggregated synapse number) and re-tallies the sign
predictions per synapse, restricted to the **neuron->neuron** scope Fenyves actually scored (muscle
NMJs are CIRCE's documented extension beyond him, excluded here for like-for-like).

Question: is the per-connection E:I gap vs Fenyves (ours 1.51 raw vs his 3.52) a counting-unit
artifact? Answer (see README Result 6): no. Weighting by synapse count moves raw E:I the *wrong*
way (1.51 -> 1.33) because the strongest edges are disproportionately complex; the gap is the
complex-handling convention + expression density, exactly as Results 1/4 found. Once complex is
receptor-ratio-resolved, per-synapse E:I = 3.13, beside Fenyves' 3.52, in his unit and scope.

Out-of-scope presynapses (no fast transmitter, e.g. purely monoaminergic) are folded into
``unpredicted`` here, matching how Fenyves' neuron->neuron denominator labels them "no pred".

Deterministic. Run from the repo root after ``uv run cckg build``:

    uv run python -I analysis/partc_sign_validation/per_synapse_reconciliation.py
"""

from __future__ import annotations

import json
from collections import Counter
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
WANG = REPO / "data" / "wang-neurotransmitter-atlas" / "sex_neurotransmitters.csv"
OUT = Path(__file__).parent / "per_synapse_reconciliation_results.json"
PRIMARY = "cook_2019_hermaphrodite"

FENYVES_2020_PER_SYNAPSE = {  # neuron->neuron, 20,589 synapses
    "excitatory": 9070,
    "inhibitory": 2580,
    "complex": 3413,
    "unpredicted": 5526,
    "total": 20589,
    "connections": 3638,
    "E:I": 3.52,
    "predicted_incl_complex_pct": 73.2,
}


def _strip(s: str) -> str:
    return str(s).split("/")[-1]


def _ratio(e: float, i: float) -> float | None:
    return round(e / i, 2) if i else None


def main() -> None:
    conn = json.loads(CONNECTOME.read_text())
    catalog = read_receptor_catalog(DATA / "ionotropic_receptor_genes.csv")
    expr = read_receptor_expression(DATA / "ionotropic_receptor_expression.csv")
    cell_type = {c["name"]: c["cell_type"] for c in conn["cells"]}

    # NT identity exactly as the build/validate uses it: neuron-graph (nt-curated) + Wang herm.
    nt_by_cell: dict[str, str] = {
        c["name"]: str(c.get("neurotransmitter") or "") for c in conn["cells"]
    }
    for line in WANG.read_text().splitlines()[1:]:
        cell, sex, nt, *_ = line.split(",")
        if sex == "hermaphrodite":
            nt_by_cell[cell] = nt

    def ion_counts(pre: str, post: str) -> Counter:
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
        return "unknown"

    def tally(neuron_to_neuron_only: bool) -> dict:
        per_conn: Counter = Counter()
        per_syn: Counter = Counter()
        resolved_syn: Counter = Counter()
        for x in conn["connections"]:
            if x["connection_type"] != "chemical" or _strip(x["dataset"]) != PRIMARY:
                continue
            pre, post = _strip(x["pre"]), _strip(x["post"])
            if neuron_to_neuron_only and not (
                cell_type.get(pre) == "neuron" and cell_type.get(post) == "neuron"
            ):
                continue
            w = x.get("weight", 0) or 0
            tx = decode_fast_transmitters(nt_by_cell.get(pre))
            # out-of-scope (no fast transmitter) folded to unpredicted, matching Fenyves "no pred".
            outcome = (
                "unpredicted" if not tx else classify_sign(tx, expr.get(post, set()), catalog)[1]
            )
            per_conn[outcome] += 1
            per_syn[outcome] += w
            ic = ion_counts(pre, post)
            resolved_syn[resolve(ic["cation"], ic["anion"])] += w

        n_conn = sum(per_conn.values())
        n_syn = sum(per_syn.values())

        def block(t: Counter, total: float) -> dict:
            definite = t["excitatory"] + t["inhibitory"]
            incl_complex = definite + t.get("complex", 0)
            return {
                "counts": {k: round(v, 1) if isinstance(v, float) else v for k, v in t.items()},
                "E:I": _ratio(t["excitatory"], t["inhibitory"]),
                "definite_pct": round(100 * definite / total, 1) if total else None,
                "predicted_incl_complex_pct": round(100 * incl_complex / total, 1)
                if total
                else None,
            }

        return {
            "n_connections": n_conn,
            "n_synapses": round(n_syn),
            "per_connection": block(per_conn, n_conn),
            "per_synapse_raw": block(per_syn, n_syn),
            "per_synapse_receptor_ratio_resolved": block(resolved_syn, n_syn),
        }

    all_scope = tally(neuron_to_neuron_only=False)
    n2n = tally(neuron_to_neuron_only=True)

    # self-check: exc/inh/complex per-connection are scope-dependent only via unpredicted; the
    # definite+complex counts must equal validate.py's committed raw (704/465/1842) in both scopes.
    committed = {"excitatory": 704, "inhibitory": 465, "complex": 1842}
    selfcheck = all(
        all_scope["per_connection"]["counts"].get(k) == v for k, v in committed.items()
    ) and all(n2n["per_connection"]["counts"].get(k) == v for k, v in committed.items())

    results = {
        "primary_connectome": PRIMARY,
        "self_check_definite_and_complex_match_committed_raw": selfcheck,
        "scope_all_chemical_incl_muscle": all_scope,
        "scope_neuron_to_neuron": n2n,
        "reference_fenyves_2020_per_synapse": FENYVES_2020_PER_SYNAPSE,
        "note": (
            "neuron->neuron is the like-for-like scope vs Fenyves (he did not score muscle). "
            "Weighting per synapse does NOT close the raw E:I gap (1.51 per-connection -> 1.33 "
            "per-synapse, the wrong direction, because the strongest edges are disproportionately "
            "complex); the gap is complex-handling + expression density. Receptor-ratio-resolved "
            "per-synapse E:I lands beside Fenyves' 3.52."
        ),
    }
    OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    # ----- console summary -----
    f = FENYVES_2020_PER_SYNAPSE
    print(f"primary: {PRIMARY}")
    print(f"self-check (exc/inh/complex match committed raw): {selfcheck}")
    for name, s in (("ALL (incl muscle)", all_scope), ("NEURON->NEURON", n2n)):
        print(f"\n=== {name}: {s['n_connections']} connections / {s['n_synapses']} synapses ===")
        pc, ps, pr = (
            s["per_connection"],
            s["per_synapse_raw"],
            s["per_synapse_receptor_ratio_resolved"],
        )
        print(f"  per-connection   : {pc['counts']}  E:I={pc['E:I']}  def={pc['definite_pct']}%")
        print(f"  per-synapse raw  : {ps['counts']}  E:I={ps['E:I']}  def={ps['definite_pct']}%")
        print(f"  per-synapse resol: {pr['counts']}  E:I={pr['E:I']}  def={pr['definite_pct']}%")
    print(
        f"\nFenyves 2020 (neuron->neuron, per synapse): "
        f"exc={f['excitatory']} inh={f['inhibitory']} complex={f['complex']} unp={f['unpredicted']} "
        f"total={f['total']}  E:I={f['E:I']}  predicted(incl complex)={f['predicted_incl_complex_pct']}%"
    )
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
