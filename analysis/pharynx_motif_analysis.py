"""Network-motif analysis of the Cook 2020 pharyngeal connectome (from the CIRCE KG projection).

Chemical synapses = directed graph; gap junctions = undirected graph. We compute the dyad and
triad census on the chemical network and test each triad's abundance against 1000 degree-preserving
random networks (Maslov-Sneppen directed edge swaps), reporting Z-scores and the triad significance
profile (Milo et al. 2004). Gap junctions get an undirected treatment (triangles/clustering).
"""

import json
import math
import random
from collections import Counter
from pathlib import Path

import networkx as nx

random.seed(20260714)  # deterministic run (fixed seed for reproducibility)

# Repo-relative: run `cckg export` first so outputs/neuron-graph-pharynx/ exists.
# Run me with:  uv run --with networkx python analysis/pharynx_motif_analysis.py
ROOT = Path(__file__).resolve().parents[1]
CONN = ROOT / "outputs" / "neuron-graph-pharynx" / "connections.json"
CELLS = ROOT / "outputs" / "neuron-graph-pharynx" / "cells.json"

conns = json.load(open(CONN))
cells = json.load(open(CELLS))


# --- Node classification (neuron / muscle / marginal / gland / other) from name conventions ------
def kind(name: str) -> str:
    n = name
    if n[:2] in ("pm",):
        return "muscle"
    if n[:2] in ("mc",):
        return "marginal"
    if n[:1] == "g" and n[1:2].isdigit():
        return "gland"
    return "neuron"


node_kind = {c["name"]: kind(c["name"]) for c in cells}
print("=== Cook 2020 pharyngeal connectome ===")
print("cells:", len(cells), "|", dict(Counter(node_kind.values())))

# --- Build graphs (binary; drop self-loops for triad analysis) ----------------------------------
D = nx.DiGraph()  # chemical, directed
U = nx.Graph()  # gap junctions, undirected
D.add_nodes_from(node_kind)
U.add_nodes_from(node_kind)
selfloops = 0
for c in conns:
    pre, post, t = c["pre"], c["post"], c["type"]
    if pre == post:
        selfloops += 1
        continue
    if t == "chemical":
        D.add_edge(pre, post)
    elif t == "electrical":
        U.add_edge(pre, post)

print(f"chemical: {D.number_of_edges()} directed edges (self-loops dropped: {selfloops})")
print(f"gap junctions: {U.number_of_edges()} undirected edges")

n = D.number_of_nodes()
m = D.number_of_edges()
print("\n--- Chemical directed network ---")
print(f"nodes={n}  edges={m}  density={m / (n * (n - 1)):.4f}  mean out/in-degree={m / n:.2f}")

# Dyad census + reciprocity
mutual = sum(1 for u, v in D.edges() if D.has_edge(v, u)) // 2
asym = m - 2 * mutual
null = n * (n - 1) // 2 - mutual - asym
print(f"dyads: mutual={mutual}  asymmetric={asym}  null={null}  reciprocity={2 * mutual / m:.3f}")

# --- Triad census (16 Holland-Leinhardt classes) ------------------------------------------------
real = nx.triadic_census(D)
TRIADS = [
    "003",
    "012",
    "102",
    "021D",
    "021U",
    "021C",
    "111D",
    "111U",
    "030T",
    "030C",
    "201",
    "120D",
    "120U",
    "120C",
    "210",
    "300",
]
# Connected-triad classes only (>=2 edges) are meaningful for motif significance.
CONNECTED = [
    "021D",
    "021U",
    "021C",
    "111D",
    "111U",
    "030T",
    "030C",
    "201",
    "120D",
    "120U",
    "120C",
    "210",
    "300",
]
NAMES = {
    "021D": "divergent (a<-b->c)",
    "021U": "convergent (a->b<-c)",
    "021C": "chain (a->b->c)",
    "111D": "111D",
    "111U": "111U",
    "030T": "feedforward loop (FFL)",
    "030C": "3-cycle (feedback)",
    "201": "201",
    "120D": "120D",
    "120U": "120U",
    "120C": "120C",
    "210": "210",
    "300": "fully mutual (clique)",
}


# --- Degree-preserving null model: directed edge swaps (preserve in+out degree of every node) ----
def randomize(g: nx.DiGraph, swaps_per_edge: int = 10) -> nx.DiGraph:
    r = g.copy()
    # sorted() so the choice pool order is independent of PYTHONHASHSEED -> fully reproducible.
    edges = sorted(r.edges())
    target = swaps_per_edge * len(edges)
    done = 0
    attempts = 0
    while done < target and attempts < target * 20:
        attempts += 1
        (a, b), (c, d) = random.choice(edges), random.choice(edges)
        if len({a, b, c, d}) < 4:
            continue
        if r.has_edge(a, d) or r.has_edge(c, b):
            continue
        r.remove_edge(a, b)
        r.remove_edge(c, d)
        r.add_edge(a, d)
        r.add_edge(c, b)
        # Maintain the pool in place (deterministic given the seed) — re-sorting each swap is O(n log n).
        edges.remove((a, b))
        edges.remove((c, d))
        edges.append((a, d))
        edges.append((c, b))
        done += 1
    return r


N_NULL = 1000
null_counts = {t: [] for t in TRIADS}
for i in range(N_NULL):
    rc = nx.triadic_census(randomize(D))
    for t in TRIADS:
        null_counts[t].append(rc[t])

print(f"\n--- Triad significance profile ({N_NULL} degree-preserving nulls) ---")
print(f"{'triad':6} {'name':28} {'real':>6} {'null_mean':>10} {'null_sd':>8} {'Z':>8}")
zscores = {}
for t in CONNECTED:
    vals = null_counts[t]
    mean = sum(vals) / len(vals)
    sd = (sum((x - mean) ** 2 for x in vals) / len(vals)) ** 0.5
    z = (real[t] - mean) / sd if sd > 0 else 0.0
    zscores[t] = z
    print(f"{t:6} {NAMES.get(t, t):28} {real[t]:>6} {mean:>10.1f} {sd:>8.2f} {z:>+8.2f}")

# Triad significance profile (normalized Z, Milo 2004)
norm = math.sqrt(sum(z * z for z in zscores.values())) or 1.0
tsp = {t: zscores[t] / norm for t in CONNECTED}

print("\n--- Motif verdict (|Z|>=2 significant) ---")
over = sorted([t for t in CONNECTED if zscores[t] >= 2], key=lambda t: -zscores[t])
under = sorted([t for t in CONNECTED if zscores[t] <= -2], key=lambda t: zscores[t])
print("OVER-represented:", [(t, NAMES.get(t, t), round(zscores[t], 1)) for t in over])
print("UNDER-represented:", [(t, NAMES.get(t, t), round(zscores[t], 1)) for t in under])

# --- Hubs (chemical) ----------------------------------------------------------------------------
print("\n--- Top hubs (chemical out / in degree) ---")
outdeg = sorted(D.out_degree(), key=lambda x: -x[1])[:8]
indeg = sorted(D.in_degree(), key=lambda x: -x[1])[:8]
print("out:", [(u, d, node_kind[u]) for u, d in outdeg])
print("in :", [(u, d, node_kind[u]) for u, d in indeg])

# --- Gap junction undirected network ------------------------------------------------------------
print("\n--- Gap junction undirected network ---")
un = U.number_of_nodes()
um = U.number_of_edges()
tri = sum(nx.triangles(U).values()) // 3
print(
    f"nodes-with-gj={sum(1 for _, d in U.degree() if d > 0)}  edges={um}  triangles={tri}"
    f"  transitivity={nx.transitivity(U):.3f}"
)
gj_deg = sorted(U.degree(), key=lambda x: -x[1])[:8]
print("top gj degree:", [(u, d, node_kind[u]) for u, d in gj_deg])

# --- Neuron-only subnetwork (isolate inter-neuronal motifs from pure-sink muscles) --------------
neurons = [c for c in node_kind if node_kind[c] == "neuron"]
Dn = D.subgraph(neurons).copy()
Dn.remove_edges_from(nx.selfloop_edges(Dn))
nn, mn = Dn.number_of_nodes(), Dn.number_of_edges()
mut_n = sum(1 for u, v in Dn.edges() if Dn.has_edge(v, u)) // 2
realn = nx.triadic_census(Dn)
nulln = {t: [] for t in TRIADS}
for _ in range(N_NULL):
    rc = nx.triadic_census(randomize(Dn))
    for t in TRIADS:
        nulln[t].append(rc[t])
zn = {}
for t in CONNECTED:
    vals = nulln[t]
    mean = sum(vals) / len(vals)
    sd = (sum((x - mean) ** 2 for x in vals) / len(vals)) ** 0.5
    zn[t] = (realn[t] - mean) / sd if sd > 0 else 0.0
print(
    f"\n=== NEURON-ONLY subnetwork ({nn} neurons, {mn} chemical edges, "
    f"reciprocity={2 * mut_n / mn:.3f}) ==="
)
print("OVER :", [(t, NAMES.get(t, t), round(zn[t], 1)) for t in CONNECTED if zn[t] >= 2])
print("UNDER:", [(t, NAMES.get(t, t), round(zn[t], 1)) for t in CONNECTED if zn[t] <= -2])

# --- Persist results for the report -------------------------------------------------------------
out = {
    "neuron_only": {
        "nodes": nn,
        "edges": mn,
        "reciprocity": round(2 * mut_n / mn, 3),
        "z": {t: round(zn[t], 2) for t in CONNECTED},
        "real": {t: realn[t] for t in TRIADS},
    },
    "nodes": n,
    "chem_edges": m,
    "gj_edges": um,
    "node_kinds": dict(Counter(node_kind.values())),
    "density": m / (n * (n - 1)),
    "reciprocity": 2 * mutual / m,
    "dyads": {"mutual": mutual, "asym": asym, "null": null},
    "triads_real": {t: real[t] for t in TRIADS},
    "z": {t: round(zscores[t], 2) for t in CONNECTED},
    "tsp": {t: round(tsp[t], 3) for t in CONNECTED},
    "over": over,
    "under": under,
    "hubs_out": [(u, d) for u, d in outdeg],
    "hubs_in": [(u, d) for u, d in indeg],
    "gj_triangles": tri,
    "gj_transitivity": round(nx.transitivity(U), 3),
    "names": NAMES,
}
json.dump(out, open(ROOT / "analysis" / "pharynx_motif_results.json", "w"), indent=2)
print("\nwrote analysis/pharynx_motif_results.json")
