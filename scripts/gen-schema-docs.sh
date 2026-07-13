#!/bin/sh
#
# Generate the published schema artifacts from the LinkML source of truth
# (src/celegans_connectome_kg/schema/connectome.yaml):
#
#   docs/schema/connectome.owl.ttl      OWL (schema/TBox) for triplestores & reasoners
#   docs/schema/connectome.schema.json  JSON Schema for validating data
#   docs/schema/*.md                     browsable schema docs (MkDocs source; gitignored)
#
# The OWL + JSON Schema are committed (small, downloadable). The Markdown is
# regenerated here and by CI, then published as a site via MkDocs + GitHub Pages
# (see .github/workflows/schema-docs.yml). Build the site locally with:
#   uv run --with mkdocs-material mkdocs serve
#
# Usage: sh scripts/gen-schema-docs.sh
set -e

root=$(cd "$(dirname "$0")/.." && pwd)
schema="$root/src/celegans_connectome_kg/schema/connectome.yaml"
out="$root/docs/schema"
mkdir -p "$out"

echo "==> OWL              -> docs/schema/connectome.owl.ttl"
uv run python -m linkml.generators.owlgen "$schema" > "$out/connectome.owl.ttl"

echo "==> JSON Schema      -> docs/schema/connectome.schema.json"
uv run python -m linkml.generators.jsonschemagen "$schema" > "$out/connectome.schema.json"

echo "==> Markdown docs    -> docs/schema/*.md"
uv run python -m linkml.generators.docgen -d "$out" "$schema"

echo "==> Done. Committed artifacts: connectome.owl.ttl, connectome.schema.json."
echo "    Preview the docs site: uv run --with mkdocs-material mkdocs serve"
