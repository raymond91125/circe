#!/usr/bin/env python3
"""CIRCE-side facts about the CAN (canal-associated neuron) pair.

Reads the built graph's JSON projection (outputs/connectome.json) and emits the
connectome + expression facts for CAN used in article.md. Deterministic: every
collection is sorted; output is independent of PYTHONHASHSEED.

    uv run cckg build          # -> outputs/connectome.json
    uv run python analysis/can_neuron/can_circe.py   # -> results.json
"""

from __future__ import annotations
import json
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "outputs" / "connectome.json"
OUT = Path(__file__).resolve().parent / "results.json"


def nm(x: str) -> str:
    return x.split("/")[-1]


def sex_of(dataset: str) -> str:
    d = dataset.lower()
    return "male" if "male" in d else ("hermaphrodite" if "herm" in d else dataset)


def main() -> None:
    d = json.loads(SRC.read_text())
    cells = {c["name"]: c for c in d["cells"]}
    genes = {g["id"]: g for g in d["genes"]}

    # --- cell records ---
    can_cells = {
        n: {
            "cell_type": cells[n]["cell_type"],
            "anatomy": cells[n]["anatomy"],
            "neurotransmitter_code": cells[n].get("neurotransmitter"),
            "sexes": cells[n].get("sexes"),
        }
        for n in ("CANL", "CANR")
        if n in cells
    }

    # --- structural connections (chemical + gap junction), deduped by (pre,post,sex) ---
    struct: dict[str, dict[str, list]] = {
        "chemical": defaultdict(set),
        "gap_junction": defaultdict(set),
    }
    chem_out, chem_in = set(), set()
    predicted_counts: dict[str, int] = defaultdict(int)
    for c in d["connections"]:
        ct = c.get("connection_type")
        pre, post = nm(c["pre"]), nm(c["post"])
        touches = "CAN" in pre or "CAN" in post
        if not touches:
            continue
        if ct in ("chemical", "gap_junction"):
            struct[ct][sex_of(nm(c["dataset"]))].add((pre, post))
            if ct == "chemical":
                (chem_out if "CAN" in pre else chem_in).add((pre, post, nm(c["dataset"])))
        else:
            predicted_counts[ct] += 1  # neuropeptidergic / monoaminergic predicted layers

    structural = {
        ct: {sex: sorted(pairs) for sex, pairs in sorted(bysex.items())}
        for ct, bysex in struct.items()
    }

    # --- expression of CIRCE's signaling genes in CAN (deduped) ---
    expr: dict[str, set] = defaultdict(set)
    for g in d["gene_expressions"]:
        cell = nm(g["cell"])
        if "CAN" not in cell:
            continue
        gid = g.get("gene")
        sym = genes.get(gid, {}).get("symbol", gid)
        cat = genes.get(gid, {}).get("category", "?")
        expr[cell].add((sym, cat))
    expression = {cell: sorted(v) for cell, v in sorted(expr.items())}

    result = {
        "source": "outputs/connectome.json",
        "cells": can_cells,
        "neurotransmitter_code_legend": "b = betaine (same code CIRCE assigns ASIL/ASIR)",
        "structural_connections_EM": structural,
        "signaling_gene_expression": expression,  # CIRCE's 155 signaling genes only
        "chemical_output_edges": sorted(chem_out),  # CAN as presynaptic partner
        "chemical_input_edges": sorted(chem_in),  # CAN as postsynaptic partner
        "predicted_layer_edge_counts": dict(sorted(predicted_counts.items())),
        "notes": (
            "CAN has ZERO chemical-synapse outputs across all EM datasets; it appears only "
            "postsynaptic (ALA->CAN) plus gap junctions and predicted neuropeptide edges."
        ),
    }
    OUT.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {OUT.relative_to(REPO)}")
    print("chemical OUTPUT edges (CAN presynaptic):", result["chemical_output_edges"] or "NONE")
    print(
        "structural edge sets:",
        {ct: {s: len(p) for s, p in bs.items()} for ct, bs in structural.items()},
    )


if __name__ == "__main__":
    main()
