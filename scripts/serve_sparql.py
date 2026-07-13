#!/usr/bin/env python3
"""Read-only SPARQL endpoint over the CIRCE knowledge graph.

Loads an RDF/Turtle build (default: outputs/connectome.ttl) into an in-memory
Oxigraph store and serves:

  GET  /                      a minimal query console (textarea + sample queries)
  GET|POST /query             SPARQL query -> SPARQL-JSON (default), CSV, TSV, or
                              Turtle (for CONSTRUCT/DESCRIBE); ?format= or Accept
  GET  /queries               JSON list of shipped sample-query names
  GET  /queries/<name>        raw .rq text of a sample query

Read-only: SPARQL Update is refused. Common CIRCE prefixes are auto-prepended, so
queries can use cckg:, WBbt:, dcterms:, rdf:, rdfs:, xsd: without declaring them.

Usage:
  uv run python scripts/serve_sparql.py [--ttl PATH] [--host H] [--port P]
"""

from __future__ import annotations

import argparse
import io
import json
import re
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

import pyoxigraph as ox

ROOT = Path(__file__).resolve().parents[1]
QUERY_DIR = ROOT / "src" / "celegans_connectome_kg" / "verify" / "queries"

PREFIXES = """\
PREFIX cckg: <https://wormbase.org/resources/connectome/>
PREFIX WBbt: <http://purl.obolibrary.org/obo/WBbt_>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
"""

_UPDATE = re.compile(r"\b(INSERT|DELETE|LOAD|CLEAR|DROP|CREATE|ADD|MOVE|COPY)\b", re.I)

_RESULT_FORMATS = {
    "json": (ox.QueryResultsFormat.JSON, "application/sparql-results+json"),
    "csv": (ox.QueryResultsFormat.CSV, "text/csv"),
    "tsv": (ox.QueryResultsFormat.TSV, "text/tab-separated-values"),
    "xml": (ox.QueryResultsFormat.XML, "application/sparql-results+xml"),
}

PAGE = """<!doctype html><meta charset=utf-8>
<title>CIRCE SPARQL</title>
<style>
 body{font:14px/1.5 system-ui,sans-serif;margin:24px;max-width:1000px}
 h1{font-size:20px} textarea{width:100%;height:200px;font-family:monospace;font-size:13px}
 button,select{font-size:14px;padding:4px 8px} pre{background:#f5f7f8;padding:12px;overflow:auto;max-height:50vh}
 .samples button{margin:2px;background:#eef;border:1px solid #ccd;border-radius:4px;cursor:pointer}
 a{color:#2b6cb0}
</style>
<h1>CIRCE knowledge graph &mdash; SPARQL</h1>
<p>Read-only endpoint. Common prefixes (<code>cckg:</code>, <code>WBbt:</code>, &hellip;) are auto-added.
Schema docs: <a href="https://raymond91125.github.io/circe/" target=_blank>raymond91125.github.io/circe</a>.</p>
<div class=samples id=samples></div>
<textarea id=q>SELECT ?type (COUNT(*) AS ?n) WHERE {
  ?c a cckg:Connection ; cckg:connection_type ?type .
} GROUP BY ?type ORDER BY DESC(?n)</textarea>
<p><button onclick=run()>Run</button>
 <select id=fmt><option value=json>JSON</option><option value=csv>CSV</option></select></p>
<pre id=out>Results appear here.</pre>
<script>
fetch('queries').then(r=>r.json()).then(names=>{
 const d=document.getElementById('samples');
 names.forEach(n=>{const b=document.createElement('button');b.textContent=n;
  b.onclick=()=>fetch('queries/'+n).then(r=>r.text()).then(t=>document.getElementById('q').value=t);
  d.appendChild(b);});});
function run(){
 const q=document.getElementById('q').value, f=document.getElementById('fmt').value;
 document.getElementById('out').textContent='Running…';
 fetch('query?format='+f,{method:'POST',headers:{'Content-Type':'application/sparql-query'},body:q})
  .then(async r=>{const t=await r.text();document.getElementById('out').textContent=r.ok?t:('Error '+r.status+'\\n'+t);})
  .catch(e=>document.getElementById('out').textContent=e);}
</script>
"""


class Handler(BaseHTTPRequestHandler):
    store: ox.Store = None  # set in main

    def _send(self, code, body, ctype="text/plain; charset=utf-8"):
        if isinstance(body, str):
            body = body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self):
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def do_GET(self):
        u = urlparse(self.path)
        if u.path == "/":
            return self._send(200, PAGE, "text/html; charset=utf-8")
        if u.path == "/queries":
            names = sorted(p.stem for p in QUERY_DIR.glob("*.rq"))
            return self._send(200, json.dumps(names), "application/json")
        if u.path.startswith("/queries/"):
            f = QUERY_DIR / (u.path[len("/queries/") :] + ".rq")
            if not f.is_file() or f.parent != QUERY_DIR:
                return self._send(404, "no such query")
            return self._send(200, f.read_text())
        if u.path == "/query":
            q = parse_qs(u.query).get("query", [None])[0]
            fmt = parse_qs(u.query).get("format", ["json"])[0]
            return self._run(q, fmt)
        self._send(404, "not found")

    def do_POST(self):
        u = urlparse(self.path)
        if u.path != "/query":
            return self._send(404, "not found")
        n = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(n).decode()
        ctype = self.headers.get("Content-Type", "")
        q = body if "sparql-query" in ctype else parse_qs(body).get("query", [None])[0]
        fmt = parse_qs(u.query).get("format", ["json"])[0]
        self._run(q, fmt)

    def _run(self, query, fmt):
        if not query or not query.strip():
            return self._send(400, "empty query")
        # Friendly early rejection of updates. Strip line comments first so update keywords
        # inside comments (e.g. "add") don't trip it; read-only is guaranteed regardless,
        # since store.query() below raises on any update operation.
        if _UPDATE.search(re.sub(r"(?m)#.*$", "", query)):
            return self._send(403, "read-only endpoint: SPARQL Update is not allowed")
        try:
            res = self.store.query(PREFIXES + query)
        except Exception as e:  # noqa: BLE001 — surface parse/eval errors to the caller
            return self._send(400, f"query error: {e}")
        buf = io.BytesIO()
        if isinstance(res, ox.QueryTriples):  # CONSTRUCT / DESCRIBE
            res.serialize(output=buf, format=ox.RdfFormat.TURTLE)
            return self._send(200, buf.getvalue(), "text/turtle")
        qfmt, ctype = _RESULT_FORMATS.get(fmt.lower(), _RESULT_FORMATS["json"])
        res.serialize(output=buf, format=qfmt)
        self._send(200, buf.getvalue(), ctype + "; charset=utf-8")

    def log_message(self, *_):  # quiet
        pass


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ttl", default=str(ROOT / "outputs" / "connectome.ttl"))
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=7878)
    a = ap.parse_args()

    print(f"loading {a.ttl} ...")
    store = ox.Store()
    store.load(path=a.ttl, format=ox.RdfFormat.TURTLE)
    print(f"loaded {len(store):,} triples")
    Handler.store = store
    srv = ThreadingHTTPServer((a.host, a.port), Handler)
    print(f"CIRCE SPARQL endpoint: http://{a.host}:{a.port}/  (Ctrl-C to stop)")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
