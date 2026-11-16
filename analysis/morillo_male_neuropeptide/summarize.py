"""Summary statistics for the derived male long-range neuropeptide connectome (Phase C).

Deterministic (sorted, no RNG; independent of PYTHONHASHSEED). Reads the built graph
(outputs/connectome.json) and the source CSVs, writes results.json next to this file, and prints a
short report. Run after `uv run cckg build`:

    python analysis/morillo_male_neuropeptide/summarize.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CONN = REPO / "outputs" / "connectome.json"
MECH = REPO / "data" / "ripoll-2023-neuropeptide" / "mechanistic"
MALE = REPO / "data" / "morillo-2025-male"
DS = "cckg:dataset/morillo_2025_male_neuropeptide_lr"


def _strip(s: str) -> str:
    return str(s).split("/")[-1]


def main() -> None:
    d = json.load(open(CONN))
    cells = {c["id"]: c for c in d["cells"]}
    male_only = {c["id"] for c in d["cells"] if [str(s) for s in c.get("sexes", [])] == ["male"]}

    edges = [c for c in d["connections"] if c["dataset"] == DS]
    endpoints = {c["pre"] for c in edges} | {c["post"] for c in edges}
    weights = [int(c["weight"]) for c in edges]
    ms_edges = [c for c in edges if c["pre"] in male_only or c["post"] in male_only]
    ms_neurons = sorted(cells[n]["name"] for n in (endpoints & male_only))

    # Pair coverage: how many of the 92 Beets pairs have BOTH ligand and receptor in the male atlas.
    mlig = {r["symbol"] for r in csv.DictReader(open(MALE / "morillo_genes.csv"))}
    mrec = {r["symbol"] for r in csv.DictReader(open(MALE / "morillo_receptor_genes.csv"))}
    pairs = list(csv.DictReader(open(MECH / "npp_gpcr_pairs.csv")))
    usable = [p for p in pairs if p["ligand"] in mlig and p["gpcr"] in mrec]

    # ins-39 expression-layer validation (Haque 2024/2025 ground truth: male-biased, ciliated sensory)
    def expr_cells(dataset_suffix: str, symbol_wb: str) -> list[str]:
        # cell id is "cckg:cell/NAME" (split "/"); gene id is "WB:WBGene..." (split ":").
        return sorted(
            _strip(e["cell"])
            for e in d["gene_expressions"]
            if _strip(e["dataset"]) == dataset_suffix and str(e["gene"]).split(":")[-1] == symbol_wb
        )

    ins39_wb = "WBGene00017668"
    ins39_male = expr_cells("morillo_2025_male", ins39_wb)
    ins39_herm = expr_cells("morillo_2025_hermaphrodite", ins39_wb)

    out = {
        "dataset": _strip(DS),
        "edges": len(edges),
        "nodes": len(endpoints),
        "weight_min": min(weights),
        "weight_max": max(weights),
        "weight_total": sum(weights),
        "edges_touching_male_specific_neuron": len(ms_edges),
        "male_specific_neurons_in_network": ms_neurons,
        "pairs_total": len(pairs),
        "pairs_usable_in_male": len(usable),
        "ligands_in_male_atlas": sum(p["ligand"] in mlig for p in pairs),
        "receptors_in_male_atlas": sum(p["gpcr"] in mrec for p in pairs),
        "ins39_male_cells": ins39_male,
        "ins39_herm_cells": ins39_herm,
    }
    (Path(__file__).parent / "results.json").write_text(json.dumps(out, indent=2) + "\n")
    for k, v in out.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
