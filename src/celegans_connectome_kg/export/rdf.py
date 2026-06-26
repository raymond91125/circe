"""Phase 3 export: serialize the LinkML Connectome.

RDF/OWL (Turtle) is the primary serialization (triplestore / SPARQL / research). The
LinkML-native JSON written by the build stage round-trips back into the data model here, so
``build`` and ``export`` are independently runnable: build → ``connectome.json`` → export →
``connectome.ttl``.
"""

from __future__ import annotations

from pathlib import Path

from linkml_runtime.dumpers import json_dumper, rdflib_dumper
from linkml_runtime.loaders import json_loader

from celegans_connectome_kg.build.datamodel import datamodel, schemaview


def write_json(connectome: object, path: Path) -> None:
    """Write the LinkML-native JSON (the build artifact / export input)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json_dumper.dumps(connectome))


def load_json(path: Path) -> object:
    """Load a Connectome back from LinkML-native JSON."""
    return json_loader.load(str(path), target_class=datamodel().Connectome)


def to_turtle(connectome: object) -> str:
    """Serialize a Connectome to RDF/Turtle."""
    return rdflib_dumper.dumps(connectome, schemaview())


def write_turtle(connectome: object, path: Path) -> None:
    """Write the Connectome as RDF/Turtle."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(to_turtle(connectome))
