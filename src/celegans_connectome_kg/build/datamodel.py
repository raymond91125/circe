"""Access to the LinkML data model and SchemaView.

The schema (``schema/connectome.yaml``) is the single source of truth; the Python classes
are compiled from it at runtime (cached per process) rather than committed, so they can
never drift from the schema.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from types import ModuleType

from linkml.generators.pythongen import PythonGenerator
from linkml_runtime import SchemaView

SCHEMA_PATH = Path(__file__).resolve().parents[1] / "schema" / "connectome.yaml"


@lru_cache(maxsize=1)
def datamodel() -> ModuleType:
    """Compile and return the schema's Python module (Connectome, Cell, …)."""
    return PythonGenerator(str(SCHEMA_PATH)).compile_module()


@lru_cache(maxsize=1)
def schemaview() -> SchemaView:
    """Return a SchemaView over the schema (used by the RDF dumper)."""
    return SchemaView(str(SCHEMA_PATH))
