"""Indirect-connectivity paths between two cells in the CIRCE knowledge graph.

Answers the two questions the connectome graph is built for:
  * are cell A and cell B connected *indirectly* (through other cells)?
  * what are the paths -- the shortest one, and all simple paths up to a hop limit?

SPARQL answers *reachability* (``queries/01_reachability.rq``) and short, fixed-length paths
(``queries/02_two_hop_paths.rq``), but a property path cannot return an arbitrary-length path, the
shortest path, or a weight-thresholded path, and it cannot constrain the connection type mid-path.
This tool does: it pulls the requested connection type(s)' edge set from the graph via pyoxigraph
(the store ``cckg verify`` uses) and runs the path algorithms in networkx.

``--type`` selects which edges to traverse -- one type, a comma-separated list, or the alias
``synaptic`` (= ``chemical,gap_junction``, i.e. the full wiring graph: chemical + electrical
synapses). With several types in one graph a path may mix hop types (e.g. a chemical hop then a gap
hop); each hop is labelled with its type. ``chemical`` / ``neuropeptidergic`` / ``monoaminergic``
are directed (pre->post); ``gap_junction`` / ``functional`` are symmetric and loaded both ways.

Traverse deliberately: the predicted neuropeptide layer alone has ~125k edges, so throwing every
type into one graph connects almost everything to everything. ``synaptic`` is the meaningful
"are they wired together" question; single types answer channel-specific ones.

Determinism: edges are aggregated across datasets (and across types) by MAX weight, the graph is
built in sorted edge order, and every emitted path list is sorted -- so output is independent of
PYTHONHASHSEED and of networkx iteration order.

Run (networkx is an ad-hoc dependency, as in ``pharynx_communities.py``):
  uv run --with networkx python analysis/connectivity_paths/paths.py AWCL AVBL --type synaptic
  uv run --with networkx python analysis/connectivity_paths/paths.py AWCL AVBL --type chemical
  uv run --with networkx python analysis/connectivity_paths/paths.py AWCL AVBL \\
      --type synaptic --min-weight 'chemical=5,gap_junction=2' --max-hops 3 --out results.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx
from pyoxigraph import RdfFormat, Store

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TTL = REPO / "outputs" / "connectome.ttl"

#: Connection types that are undirected (an edge implies both directions).
SYMMETRIC = {"gap_junction", "functional"}

#: Every connection_type present in the graph.
VALID_TYPES = {"chemical", "gap_junction", "neuropeptidergic", "monoaminergic", "functional"}

#: Convenience names that expand to a set of connection types.
TYPE_ALIASES = {"synaptic": ("chemical", "gap_junction")}

#: Edges of the requested type(s), weight aggregated across datasets (MAX), per (pre, post, type).
EDGE_QUERY = """
PREFIX cckg: <https://wormbase.org/resources/connectome/>
PREFIX dcterms: <http://purl.org/dc/terms/>
SELECT ?pre ?post ?type (MAX(?w) AS ?weight) WHERE {
  VALUES ?type { %s }
  ?c a cckg:Connection ; cckg:connection_type ?type ;
     cckg:pre ?preC ; cckg:post ?postC ; cckg:weight ?w .
  ?preC dcterms:title ?pre .
  ?postC dcterms:title ?post .
}
GROUP BY ?pre ?post ?type
ORDER BY ?pre ?post ?type
"""


def parse_types(arg: str) -> list[str]:
    """Expand a ``--type`` argument (single / comma list / alias) to an ordered, de-duped list."""
    expanded: list[str] = []
    for tok in (t.strip() for t in arg.split(",")):
        expanded.extend(TYPE_ALIASES.get(tok, (tok,)) if tok else ())
    seen: set[str] = set()
    out: list[str] = []
    for t in expanded:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out


def parse_min_weights(arg: str, types: list[str]) -> dict[str, float]:
    """A bare number is the floor for every type; 'type=val,...' sets per-type floors (others 0)."""
    arg = arg.strip()
    if "=" not in arg:
        try:
            return dict.fromkeys(types, float(arg or 0))
        except ValueError:
            raise SystemExit(f"--min-weight: {arg!r} is not a number") from None
    out = dict.fromkeys(types, 0.0)
    for tok in filter(None, (t.strip() for t in arg.split(","))):
        key, sep, val = tok.partition("=")
        key = key.strip()
        if not sep or key not in types:
            raise SystemExit(f"--min-weight: {tok!r} must be TYPE=NUMBER with TYPE among {types}")
        try:
            out[key] = float(val)
        except ValueError:
            raise SystemExit(f"--min-weight: {val!r} is not a number") from None
    return out


def _add_edge(g: nx.DiGraph, u: str, v: str, ctype: str, w: float) -> None:
    """Add/merge one directed edge, keeping MAX weight and the set of contributing types."""
    if g.has_edge(u, v):
        d = g[u][v]
        d["weight"] = max(d["weight"], w)
        d["types"].add(ctype)
    else:
        g.add_edge(u, v, weight=w, types={ctype})


def load_graph(store: Store, types: list[str], min_weights: dict[str, float]) -> nx.DiGraph:
    """Build one directed graph over all requested types (symmetric types get both directions).

    Each edge is filtered by the floor for *its own* connection type (``min_weights``), so different
    types can carry different strength thresholds.
    """
    g = nx.DiGraph()
    values = " ".join(json.dumps(t) for t in types)
    for s in store.query(EDGE_QUERY % values):
        ctype, pre, post = s["type"].value, s["pre"].value, s["post"].value
        w = float(s["weight"].value)
        if w < min_weights[ctype]:
            continue
        _add_edge(g, pre, post, ctype, w)
        if ctype in SYMMETRIC:
            _add_edge(g, post, pre, ctype, w)
    return g


def describe(g: nx.DiGraph, path: list[str]) -> dict[str, object]:
    """Annotate a node path with per-hop weight and type(s)."""
    hops = [g[path[i]][path[i + 1]] for i in range(len(path) - 1)]
    weights = [h["weight"] for h in hops]
    return {
        "path": path,
        "hops": len(path) - 1,
        "weights": weights,
        "hop_types": [sorted(h["types"]) for h in hops],
        "min_weight": min(weights) if weights else None,
    }


def analyze(g: nx.DiGraph, a: str, b: str, max_hops: int) -> dict[str, object]:
    """Reachability, shortest path(s), and all simple paths up to ``max_hops`` edges."""
    out: dict[str, object] = {"direct": g.has_edge(a, b), "reachable": nx.has_path(g, a, b)}
    if not out["reachable"]:
        return out
    shortest = sorted(nx.all_shortest_paths(g, a, b))
    out["shortest_hops"] = len(shortest[0]) - 1
    out["shortest_path"] = describe(g, shortest[0])
    out["n_shortest_paths"] = len(shortest)
    simple = sorted(nx.all_simple_paths(g, a, b, cutoff=max_hops), key=lambda p: (len(p), p))
    out["paths_within_max_hops"] = [describe(g, p) for p in simple]
    out["n_paths_within_max_hops"] = len(simple)
    return out


def render(p: dict[str, object], show_types: bool) -> str:
    """One-line rendering of a described path, with hop weights (and types when mixed)."""
    path, weights, hop_types = p["path"], p["weights"], p["hop_types"]
    segs = []
    for i in range(p["hops"]):
        label = f"{weights[i]:g}"
        if show_types:
            label += ":" + "/".join(hop_types[i])
        segs.append(f"{path[i]} --{label}-->")
    return " ".join(segs) + f" {path[-1]}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("source", help="source cell name, e.g. AWCL")
    ap.add_argument("target", help="target cell name, e.g. AVBL")
    ap.add_argument(
        "--type",
        default="chemical",
        help="type(s) to traverse: a connection_type, a comma list, or 'synaptic' "
        "(=chemical,gap_junction). Default: chemical",
    )
    ap.add_argument("--max-hops", type=int, default=3, help="max edges in enumerated simple paths")
    ap.add_argument(
        "--min-weight",
        default="0",
        help="edge strength floor: one number for every type, or per-type like "
        "'chemical=5,gap_junction=2'",
    )
    ap.add_argument("--out", type=Path, help="write the result as JSON to this path")
    args = ap.parse_args()

    types = parse_types(args.type)
    unknown = [t for t in types if t not in VALID_TYPES]
    if not types or unknown:
        raise SystemExit(
            f"unknown connection type(s): {unknown or args.type!r}; "
            f"valid: {sorted(VALID_TYPES)} (or alias 'synaptic')"
        )
    min_weights = parse_min_weights(args.min_weight, types)

    if not TTL.exists():
        raise SystemExit(f"missing {TTL} -- run `uv run cckg export` first")
    store = Store()
    store.load(TTL.read_bytes(), format=RdfFormat.TURTLE)
    g = load_graph(store, types, min_weights)

    a, b = args.source, args.target
    label = args.type if len(types) == 1 else f"{args.type} [{'+'.join(types)}]"
    for n in (a, b):
        if n not in g:
            raise SystemExit(
                f"{n!r} has no {label} edges (min_weight={args.min_weight}) in the graph"
            )

    result = {
        "source": a,
        "target": b,
        "types_requested": args.type,
        "types": types,
        "min_weights": min_weights,
        "max_hops": args.max_hops,
        **analyze(g, a, b, args.max_hops),
    }
    if args.out:
        args.out.write_text(json.dumps(result, indent=2) + "\n")

    show_types = len(types) > 1
    print(f"{a} -> {b}  via {label} (min_weight={args.min_weight})")
    print(f"  direct edge: {'yes' if result['direct'] else 'no'}")
    print(f"  reachable (indirectly): {'YES' if result['reachable'] else 'no'}")
    if not result["reachable"]:
        return
    print(
        f"  shortest: {result['shortest_hops']} hops  "
        f"({render(result['shortest_path'], show_types)})  "
        f"[{result['n_shortest_paths']} shortest path(s)]"
    )
    print(f"  simple paths <= {args.max_hops} hops: {result['n_paths_within_max_hops']}")
    if result["shortest_hops"] > args.max_hops:
        print(
            f"  (none within {args.max_hops} hops because the shortest is "
            f"{result['shortest_hops']} -- they ARE connected; raise --max-hops to "
            f"{result['shortest_hops']} or more to list the routes)"
        )
    for p in result["paths_within_max_hops"][:15]:
        print(f"    {render(p, show_types)}")


if __name__ == "__main__":
    main()
