"""Match cell names to WBBT terms and emit the match report + curation work-list.

Bucketing (per docs/PLANNING.md):
  - matched    — exactly one term via a strong kind (label / exact synonym).
  - ambiguous  — several strong terms, OR only weak (related/broad/narrow) synonym hits.
  - unmatched  — no hit at all.

The ambiguous + unmatched tail is written to a flat work-list for a human to resolve. No
cell-type (neuron/muscle) classification happens here — that is a build-stage concern; the
work-list carries the raw neuron-graph fields as curation context.
"""

from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from celegans_connectome_kg.ingest.neuron_graph import CellRecord
from celegans_connectome_kg.match.wbbt import STRONG_KINDS, Hit, WBBTIndex

MatchStatus = Literal["matched", "curated", "ambiguous", "unmatched"]


@dataclass(frozen=True)
class CellMatch:
    """The result of matching one cell name against the WBBT index."""

    cell_name: str
    status: MatchStatus
    wbbt_id: str | None
    match_kind: str | None
    reason: str
    candidates: list[Hit] = field(default_factory=list)


def match_cell(name: str, index: WBBTIndex) -> CellMatch:
    """Match a single cell name, bucketing into matched / ambiguous / unmatched."""
    hits = index.lookup(name)
    if not hits:
        return CellMatch(name, "unmatched", None, None, "no lexical hit in WBBT")

    strong = {h.curie for h in hits if h.kind in STRONG_KINDS}
    if len(strong) == 1:
        curie = next(iter(strong))
        kind = next(h.kind for h in hits if h.curie == curie and h.kind in STRONG_KINDS)
        return CellMatch(name, "matched", curie, kind, f"unique {kind} match", hits)
    if len(strong) > 1:
        return CellMatch(name, "ambiguous", None, None, f"{len(strong)} strong matches", hits)
    # Only weak synonym hits (related/broad/narrow) → low confidence.
    weak_curies = {h.curie for h in hits}
    kinds = ",".join(sorted({h.kind for h in hits}))
    reason = f"only weak synonym ({kinds})"
    if len(weak_curies) > 1:
        reason = f"{len(weak_curies)} weak-synonym candidates ({kinds})"
    return CellMatch(name, "ambiguous", None, None, reason, hits)


def match_cells(
    cells: list[CellRecord],
    index: WBBTIndex,
    curation: dict[str, str] | None = None,
) -> list[CellMatch]:
    """Match every cell (sorted by name). Curated cells take precedence over the lexical result."""
    curation = curation or {}
    results = []
    for cell in sorted(cells, key=lambda c: c.name):
        if cell.name in curation:
            results.append(
                CellMatch(
                    cell.name,
                    "curated",
                    curation[cell.name],
                    "curated",
                    "manual curation (data/curation)",
                    index.lookup(cell.name),
                )
            )
        else:
            results.append(match_cell(cell.name, index))
    return results


def summarize(matches: list[CellMatch]) -> Counter[str]:
    """Count matches by status."""
    return Counter(m.status for m in matches)


def _format_candidates(hits: list[Hit]) -> str:
    """``WBbt:0004742(related);WBbt:0004013(label)`` — stable, human-readable."""
    return ";".join(f"{h.curie}({h.kind})" for h in hits)


def write_report_csv(matches: list[CellMatch], path: Path) -> None:
    """Write the full match report: one row per cell, every status."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["cell_name", "status", "wbbt_id", "match_kind", "reason", "candidates"])
        for m in matches:
            writer.writerow(
                [
                    m.cell_name,
                    m.status,
                    m.wbbt_id or "",
                    m.match_kind or "",
                    m.reason,
                    _format_candidates(m.candidates),
                ]
            )


def write_worklist_csv(matches: list[CellMatch], cells: list[CellRecord], path: Path) -> None:
    """Write the curation work-list: only ambiguous + unmatched, with curation context.

    ``resolved_wbbt_id`` is left blank for a human to fill in.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    by_name = {c.name: c for c in cells}
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "cell_name",
                "status",
                "reason",
                "candidates",
                "nemanode_type",
                "cell_class",
                "resolved_wbbt_id",
            ]
        )
        for m in matches:
            if m.status in ("matched", "curated"):
                continue
            cell = by_name.get(m.cell_name)
            writer.writerow(
                [
                    m.cell_name,
                    m.status,
                    m.reason,
                    _format_candidates(m.candidates),
                    cell.nemanode_type if cell else "",
                    cell.cell_class if cell else "",
                    "",
                ]
            )
