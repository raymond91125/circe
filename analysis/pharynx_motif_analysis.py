"""Network-motif analysis of the Cook 2020 pharyngeal connectome (from the CIRCE KG projection).

Motif detection is a two-level question, and the null model must match the level:

  * Dyadic (2-node): is the network more RECIPROCAL than chance? Tested against a null that
    preserves only in/out-degree (so reciprocity is free to vary).
  * Triadic (3-node): are there 3-node building blocks beyond what the dyad structure explains?
    Tested against a null that ALSO preserves the number of mutual (bidirectional) dyads
    (Milo et al. 2002; the standard for directed motifs, and what Cook et al. 2020 Fig 6 used).

Controlling for reciprocity at the triad level is essential: otherwise a network with genuinely
high reciprocity shows every mutual-containing triad (300, 201, 210, ...) as trivially "enriched"
-- a shadow of the dyadic reciprocity, not a triadic motif. This script reports both levels.
"""

import json
import random
from collections import Counter
from pathlib import Path

import networkx as nx

random.seed(20260714)  # deterministic run (edge pools are sorted below -> hash-seed independent)

ROOT = Path(__file__).resolve().parents[1]
conns = json.load(open(ROOT / "outputs" / "neuron-graph-pharynx" / "connections.json"))
cells = json.load(open(ROOT / "outputs" / "neuron-graph-pharynx" / "cells.json"))


# --- Node classification (neuron / muscle / marginal / gland) from name conventions -------------
def kind(name: str) -> str:
    if name[:2] == "pm":
        return "muscle"
    if name[:2] == "mc":
        return "marginal"
    if name[:1] == "g" and name[1:2].isdigit():
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
Wt = {}  # (pre, post) -> EM serial-section weight, for the worked example
selfloops = 0
for c in conns:
    pre, post, t = c["pre"], c["post"], c["type"]
    if pre == post:
        selfloops += 1
        continue
    if t == "chemical":
        D.add_edge(pre, post)
        Wt[(pre, post)] = c["synapses"]["cook_2020_pharynx"]
    elif t == "electrical":
        U.add_edge(pre, post)

TRIADS = ["003","012","102","021D","021U","021C","111D","111U",
          "030T","030C","201","120D","120U","120C","210","300"]  # fmt: skip
CONNECTED = ["021D","021U","021C","111D","111U","030T","030C","201","120D","120U","120C","210","300"]  # fmt: skip
NAMES = {
    "021D": "divergent fan-out",
    "021U": "convergent fan-in",
    "021C": "directed chain",
    "111D": "mutual pair, one incoming",
    "111U": "mutual pair, one outgoing",
    "030T": "feed-forward loop",
    "030C": "directed 3-cycle",
    "201": "mutual chain",
    "120D": "shared source, reciprocal pair",
    "120U": "reciprocal pair, shared target",
    "120C": "mutual pair in a cycle",
    "210": "near-clique",
    "300": "fully-mutual triangle",
}


def zscore(real_count, null_counts):
    mean = sum(null_counts) / len(null_counts)
    sd = (sum((v - mean) ** 2 for v in null_counts) / len(null_counts)) ** 0.5
    return ((real_count - mean) / sd if sd else 0.0), mean


def split_dyads(g):
    """(mutual undirected pairs, single directed edges) of a DiGraph."""
    mutual = sorted({tuple(sorted((u, v))) for u, v in g.edges() if g.has_edge(v, u)})
    single = sorted((u, v) for u, v in g.edges() if not g.has_edge(v, u))
    return mutual, single


# --- Null 1: preserve in/out-degree only (reciprocity free to vary) -> DYADIC test ---------------
def randomize_degree(g, k=10):
    r = g.copy()
    edges = sorted(r.edges())
    done = attempts = 0
    while done < k * len(edges) and attempts < 200 * len(edges):
        attempts += 1
        (a, b), (c, d) = random.choice(edges), random.choice(edges)
        if len({a, b, c, d}) < 4 or r.has_edge(a, d) or r.has_edge(c, b):
            continue
        r.remove_edge(a, b)
        r.remove_edge(c, d)
        r.add_edge(a, d)
        r.add_edge(c, b)
        edges.remove((a, b))
        edges.remove((c, d))
        edges.append((a, d))
        edges.append((c, b))
        done += 1
    return r


# --- Null 2: also preserve the mutual-dyad count (Milo/Cook) -> TRIADIC test ---------------------
def randomize_reciprocity(mutual, single, nodes, k=10):
    """Swap singles among singles and mutuals among mutuals -> preserves single-in, single-out
    and mutual degree per node (so reciprocity is held fixed)."""
    S, M = list(single), list(mutual)
    Sset, Mset = set(S), set(M)

    def occ(u, v):
        return (u, v) in Sset or tuple(sorted((u, v))) in Mset

    for _ in range(k * len(S)):
        i, j = random.randrange(len(S)), random.randrange(len(S))
        (a, b), (c, d) = S[i], S[j]
        if len({a, b, c, d}) < 4:
            continue
        if occ(a, d) or occ(c, b) or occ(d, a) or occ(b, c):
            continue
        Sset.discard((a, b))
        Sset.discard((c, d))
        Sset.add((a, d))
        Sset.add((c, b))
        S[i], S[j] = (a, d), (c, b)
    for _ in range(k * len(M)):
        i, j = random.randrange(len(M)), random.randrange(len(M))
        (a, b), (c, d) = M[i], M[j]
        if len({a, b, c, d}) < 4:
            continue
        n1, n2 = tuple(sorted((a, d))), tuple(sorted((c, b)))
        if n1 in Mset or n2 in Mset or any(x in Sset for x in [(a, d), (d, a), (c, b), (b, c)]):
            continue
        Mset.discard((a, b))
        Mset.discard((c, d))
        Mset.add(n1)
        Mset.add(n2)
        M[i], M[j] = n1, n2
    r = nx.DiGraph()
    r.add_nodes_from(nodes)
    for u, v in Sset:
        r.add_edge(u, v)
    for a, b in Mset:
        r.add_edge(a, b)
        r.add_edge(b, a)
    return r


N_NULL = 1000


def analyse(g, label):
    m = g.number_of_edges()
    mutual, single = split_dyads(g)
    real = nx.triadic_census(g)
    nodes = sorted(g.nodes())
    recip = 2 * len(mutual) / m

    # dyadic: mutual-dyad count vs degree-only null
    deg_nulls = [randomize_degree(g) for _ in range(N_NULL)]
    mut_deg = [len(split_dyads(rr)[0]) for rr in deg_nulls]
    zmut, mmut = zscore(len(mutual), mut_deg)

    # triadic: census vs reciprocity-preserving null
    rec_nulls = [
        nx.triadic_census(randomize_reciprocity(mutual, single, nodes)) for _ in range(N_NULL)
    ]
    z = {}
    for t in CONNECTED:
        z[t], _ = zscore(real[t], [rc[t] for rc in rec_nulls])

    print(f"\n=== {label}: {g.number_of_nodes()} nodes, {m} chemical edges ===")
    print(
        f"reciprocity {recip:.3f}  |  mutual dyads {len(mutual)} observed vs "
        f"{mmut:.1f} degree-matched (Z {zmut:+.1f})  <- DYADIC finding"
    )
    print("triad Z (reciprocity-preserving null), most enriched first:")
    for t in sorted(CONNECTED, key=lambda t: -z[t]):
        tag = "OVER" if z[t] >= 2 else ("under" if z[t] <= -2 else "")
        print(f"  {t:5} {NAMES[t]:30} obs {real[t]:4}  Z {z[t]:+5.1f}  {tag}")
    over = [t for t in CONNECTED if z[t] >= 2]
    print(f"  -> {len(over)} over-represented triplet(s): {over}")
    return {
        "nodes": g.number_of_nodes(),
        "edges": m,
        "reciprocity": round(recip, 3),
        "mut_obs": len(mutual),
        "mut_deg_null": round(mmut, 1),
        "z_mut": round(zmut, 1),
        "real": {t: real[t] for t in TRIADS},
        "z": {t: round(z[t], 1) for t in CONNECTED},
        "over": over,
    }


whole = analyse(D, "WHOLE NETWORK")
Dn = D.subgraph([c for c in node_kind if node_kind[c] == "neuron"]).copy()
Dn.remove_edges_from(nx.selfloop_edges(Dn))
neuron = analyse(Dn, "NEURON-ONLY")

# --- Hubs + gap junctions -----------------------------------------------------------------------
outdeg = sorted(D.out_degree(), key=lambda x: -x[1])[:8]
indeg = sorted(D.in_degree(), key=lambda x: -x[1])[:8]
gj_deg = sorted(U.degree(), key=lambda x: -x[1])[:8]
gj_tri = sum(nx.triangles(U).values()) // 3
print(f"\nhubs out: {[(u, d) for u, d in outdeg]}")
print(f"hubs in : {[(u, d) for u, d in indeg]}")
print(
    f"gap junctions: {U.number_of_edges()} edges, {gj_tri} triangles, "
    f"transitivity {nx.transitivity(U):.3f}"
)

# --- Worked example: a fully-mutual (300) triangle (abundant, but NOT a motif once reciprocity
#     is controlled -- it is a shadow of the pairwise reciprocity) -------------------------------
mutual, _ = split_dyads(D)
UM = nx.Graph()
UM.add_edges_from(mutual)
example = None
for a, b in mutual:
    common = sorted(set(UM[a]) & set(UM[b]))
    if common:
        x, y, z2 = a, b, common[0]
        example = [(x, y), (y, x), (y, z2), (z2, y), (x, z2), (z2, x)]
        break
print("\nexample fully-mutual triangle:")
for u, v in example:
    print(f"  {u} -> {v}  weight {Wt[(u, v)]}")

out = {
    "whole": whole,
    "neuron": neuron,
    "self_loops_dropped": selfloops,
    "node_kinds": dict(Counter(node_kind.values())),
    "gj_edges": U.number_of_edges(),
    "gj_triangles": gj_tri,
    "gj_transitivity": round(nx.transitivity(U), 3),
    "hubs_out": [[u, d] for u, d in outdeg],
    "hubs_in": [[u, d] for u, d in indeg],
    "gj_hubs": [[u, d] for u, d in gj_deg],
    "example_300": [[u, v, Wt[(u, v)]] for u, v in example],
    "names": NAMES,
}
json.dump(out, open(ROOT / "analysis" / "pharynx_motif_results.json", "w"), indent=2)
print("\nwrote analysis/pharynx_motif_results.json")
