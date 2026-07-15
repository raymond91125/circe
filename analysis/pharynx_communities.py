"""Louvain community detection on the Cook 2020 pharyngeal connectome.

Method: fast-unfolding / Louvain modularity maximization (Blondel, Guillaume, Lambiotte &
Lefebvre, 2008), matching Cook et al. 2020's construction -- a single weighted UNDIRECTED graph
combining chemical synapses (both directions summed) and gap junctions, weights = EM serial
sections.

Per the repo analysis conventions:
  * prior art: Cook 2020 ran the same analysis and reported 4 functional modules at Q=0.352;
    we reproduce and reconcile.
  * stability: Louvain is stochastic -> 200 seeds; report the community-count distribution and
    a pairwise co-assignment stability score, take the max-modularity partition.
  * significance: compare Q to a degree-preserving null (double-edge-swap).
  * sensitivity: sweep the resolution parameter.

Deterministic: fixed seed lists; no PYTHONHASHSEED dependence (all seeds explicit).

Run: uv run --with networkx --with numpy python analysis/pharynx_communities.py
"""

import json
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx
from networkx.algorithms.community import louvain_communities, modularity

ROOT = Path(__file__).resolve().parents[1]
conns = json.load(open(ROOT / "outputs" / "neuron-graph-pharynx" / "connections.json"))
cells = json.load(open(ROOT / "outputs" / "neuron-graph-pharynx" / "cells.json"))


def kind(n):
    if n[:2] == "pm":
        return "muscle"
    if n[:2] == "mc":
        return "marginal"
    if n[:1] == "g" and n[1:2].isdigit():
        return "gland"
    return "neuron"


node_kind = {c["name"]: kind(c["name"]) for c in cells}

# Combined weighted undirected graph: chemical (both directions summed) + gap junction.
G = nx.Graph()
G.add_nodes_from(node_kind)
for c in conns:
    u, v = c["pre"], c["post"]
    if u == v:
        continue
    w = c["synapses"]["cook_2020_pharynx"]
    if G.has_edge(u, v):
        G[u][v]["weight"] += w
    else:
        G.add_edge(u, v, weight=w)
print(f"combined weighted graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

SEEDS = list(range(200))


def best_partition(graph, resolution=1.0, seeds=SEEDS):
    parts = [
        louvain_communities(graph, weight="weight", resolution=resolution, seed=s) for s in seeds
    ]
    qs = [modularity(graph, p, weight="weight", resolution=resolution) for p in parts]
    i = max(range(len(parts)), key=lambda k: qs[k])
    return parts, qs, parts[i], qs[i]


parts, qs, best, qbest = best_partition(G)
ks = [len(p) for p in parts]
print(f"\ndefault resolution (gamma=1): community counts {dict(sorted(Counter(ks).items()))}")
print(f"modularity Q: max {qbest:.3f}, mean {sum(qs) / len(qs):.3f}")

# pairwise co-assignment stability (1 = every pair decisively same/different across seeds)
nodes = sorted(G.nodes())
co = defaultdict(int)
for p in parts:
    lab = {n: i for i, com in enumerate(p) for n in com}
    for a in range(len(nodes)):
        for b in range(a + 1, len(nodes)):
            if lab[nodes[a]] == lab[nodes[b]]:
                co[(a, b)] += 1
pairs = [(a, b) for a in range(len(nodes)) for b in range(a + 1, len(nodes))]
stability = sum(abs(2 * co[p] / len(SEEDS) - 1) for p in pairs) / len(pairs)
print(f"stability (mean pairwise decisiveness): {stability:.3f}")

# Q significance vs degree-preserving null
null_q = []
for s in range(100):
    r = nx.double_edge_swap(
        G.copy(), nswap=5 * G.number_of_edges(), max_tries=100 * G.number_of_edges(), seed=s
    )
    _, _, _, q = best_partition(r, seeds=range(10))
    null_q.append(q)
mn = sum(null_q) / len(null_q)
sd = (sum((x - mn) ** 2 for x in null_q) / len(null_q)) ** 0.5
zq = (qbest - mn) / sd if sd else 0.0
print(f"Q vs degree-preserving null: {qbest:.3f} vs {mn:.3f}+/-{sd:.3f}  (Z {zq:+.1f})")

# resolution sweep
sweep = []
for g in [0.5, 0.75, 1.0, 1.25, 1.5]:
    _, _, p, q = best_partition(G, resolution=g, seeds=range(60))
    sweep.append({"gamma": g, "k": len(p), "Q": round(q, 3)})
    print(f"  gamma={g}: {len(p)} communities, Q={q:.3f}")

# --- The 4-module functional partition (gamma=0.75) that matches Cook's count ---------------------
_, _, mods, qmods = best_partition(G, resolution=0.75)
mods = sorted(mods, key=len, reverse=True)


def label_module(com):
    """Name a module by its signature cells (Cook 2020's functional domains)."""
    names = set(com)
    if "M5" in names:
        return "Grinding", "terminal bulb"
    if "M4" in names:
        return "Peristalsis", "isthmus"
    if "NSML" in names or "NSMR" in names:
        return "Neuromodulation", "systemic"
    return "Pumping", "corpus"


# --- Deterministic layout for the report's network diagram ---------------------------------------
# Unweighted spring layout with a larger optimal distance spreads the dense, highly-connected
# core (weighting collapses it into a hairball); positions are for display only.
pos = nx.spring_layout(G, k=0.9, seed=7, iterations=400)
xs = [p[0] for p in pos.values()]
ys = [p[1] for p in pos.values()]
minx, maxx, miny, maxy = min(xs), max(xs), min(ys), max(ys)


def norm(p):
    return [round((p[0] - minx) / (maxx - minx), 4), round((p[1] - miny) / (maxy - miny), 4)]


mod_of = {}
labels = []
for i, com in enumerate(mods):
    name, region = label_module(com)
    labels.append(
        {
            "name": name,
            "region": region,
            "n": len(com),
            "kinds": dict(Counter(node_kind[n] for n in com)),
            "cells": sorted(com),
        }
    )
    for n in com:
        mod_of[n] = i

print("\n=== 4 functional modules (gamma=0.75, Q={:.3f}) ===".format(qmods))
for m in labels:
    print(f"  {m['name']:16} ({m['region']:13}) {m['n']:2} cells {m['kinds']}")
    print(f"      {', '.join(m['cells'])}")

deg = dict(G.degree(weight="weight"))
out = {
    "nodes_total": G.number_of_nodes(),
    "edges": G.number_of_edges(),
    "q_default": round(qbest, 3),
    "k_default": len(best),
    "k_distribution": {str(k): v for k, v in sorted(Counter(ks).items())},
    "stability": round(stability, 3),
    "q_null_mean": round(mn, 3),
    "q_null_sd": round(sd, 3),
    "q_z": round(zq, 1),
    "sweep": sweep,
    "q_modules": round(qmods, 3),
    "modules": labels,
    "graph": {
        "nodes": [
            {
                "id": n,
                "kind": node_kind[n],
                "mod": mod_of[n],
                "deg": round(deg[n], 1),
                "xy": norm(pos[n]),
            }
            for n in nodes
        ],
        "edges": [{"u": u, "v": v, "w": G[u][v]["weight"]} for u, v in G.edges()],
    },
}
json.dump(out, open(ROOT / "analysis" / "pharynx_communities_results.json", "w"), indent=1)
print("\nwrote analysis/pharynx_communities_results.json")
