"""Rewrite a Turtle file into a canonical, deterministic form (in place).

``linkml.generators.owlgen`` emits the OWL/TBox as Turtle whose blank-node labels and
statement order vary run-to-run (the OWL restrictions are anonymous nodes, and neither the
blank-node identifiers nor their serialization order are stable — not even under a fixed
``PYTHONHASHSEED``). That makes the committed ``connectome.owl.ttl`` churn by hundreds of
lines on every regeneration, drowning real schema changes in noise.

This canonicalizes the graph with rdflib's isomorphic-canonicalization (content-hashed,
stable blank-node labels; see ``rdflib.compare.to_canonical_graph``) and re-serializes,
preserving the source's namespace prefixes. The result is byte-identical across runs for
identical schema content, so the artifact only changes when the schema actually does.

Usage: python scripts/canonicalize_ttl.py <file.ttl> [<file.ttl> ...]
"""

from __future__ import annotations

import sys
from pathlib import Path

import rdflib
from rdflib.compare import to_canonical_graph


def canonicalize(path: Path) -> None:
    source = rdflib.Graph()
    source.parse(str(path), format="turtle")

    # Content-hashed canonical blank-node labels -> deterministic across runs.
    canonical = to_canonical_graph(source)

    # to_canonical_graph returns a read-only aggregate without prefix bindings; copy the
    # triples into a fresh graph carrying the source's namespaces so the output stays readable.
    out = rdflib.Graph()
    for prefix, namespace in source.namespaces():
        out.bind(prefix, namespace)
    for triple in canonical:
        out.add(triple)

    path.write_text(out.serialize(format="turtle"))


def main(argv: list[str]) -> int:
    if not argv:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    for arg in argv:
        canonicalize(Path(arg))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
