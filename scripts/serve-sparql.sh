#!/bin/sh
# Serve a read-only SPARQL endpoint over the CIRCE knowledge graph.
# Loads outputs/connectome.ttl (run `uv run cckg build && uv run cckg export` first
# if it's missing) into an in-memory Oxigraph store and serves a query console.
# Usage: sh scripts/serve-sparql.sh [--port 7878] [--ttl PATH] [--host 0.0.0.0]
exec uv run python "$(dirname "$0")/serve_sparql.py" "$@"
