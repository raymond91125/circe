"""Ingest the Cook et al. 2019 whole-animal connectomes (both sexes) from SI 5.

Parses the labeled adjacency-matrix sheets of the pinned
``data/cook-2019-connectome/SI5_connectome_adjacency_matrices_corrected_2020.xlsx``
(see that directory's MANIFEST) into normalized connection + cell records, one dataset per sex
(``cook_2019_hermaphrodite``, ``cook_2019_male``).

Sheets used per sex: ``{sex} chemical`` (directed: row = pre, col = post) and
``{sex} gap jn symmetric`` (undirected — emitted once per unordered pair to match the KG's
single-record-per-gap-junction convention; the mirror entry is skipped).

IMPORTANT — weight semantics: Cook weights are **EM serial-section counts** (synapse number ×
size), NOT raw synapse counts like the neuron-graph datasets (White 1986, Witvliet). They are
carried through verbatim; downstream must keep them tagged to their dataset/source so weights are
never silently compared across sources.

Sheet layout: each matrix has region-grouping headers in the outer rows/columns; the cell-name
header row and header column are found by scanning for the row/column with the most cell-like
labels (robust to the outer region rows/cols and to blank separators).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import openpyxl

from celegans_connectome_kg.ingest.neuron_graph import ConnectionRecord

#: Pinned SI 5 workbook (relative to the repo root).
DEFAULT_COOK_XLSX = (
    Path(__file__).resolve().parents[3]
    / "data"
    / "cook-2019-connectome"
    / "SI5_connectome_adjacency_matrices_corrected_2020.xlsx"
)

#: sex -> dataset id.
DATASET_IDS = {
    "hermaphrodite": "cook_2019_hermaphrodite",
    "male": "cook_2019_male",
}

_CHEMICAL = "chemical"
_GAP_JUNCTION = "gap_junction"
# (connection_type, sheet-name template, symmetric?)
_SHEETS = (
    (_CHEMICAL, "{sex} chemical", False),
    (_GAP_JUNCTION, "{sex} gap jn symmetric", True),
)

# A cell/organ label: starts with a letter, short, no whitespace. Used only to *locate* the
# header row/column; extraction then takes every non-empty label in that row/column.
_CELL_LABEL = re.compile(r"^[A-Za-z][A-Za-z0-9_.\-]{0,9}$")


@dataclass(frozen=True)
class CookDataset:
    """One Cook 2019 connectome dataset (one sex)."""

    id: str
    sex: str
    name: str
    description: str


@dataclass(frozen=True)
class CookData:
    """The ingested Cook 2019 bundle: sex-tagged datasets, connections, and cell sets."""

    datasets: list[CookDataset]
    connections: list[ConnectionRecord]
    cells_by_sex: dict[str, set[str]]


def _looks_like_cell(value: object) -> bool:
    return isinstance(value, str) and bool(_CELL_LABEL.match(value.strip()))


def _find_headers(grid: list[list]) -> tuple[int, int]:
    """Return (header_row, header_col) as the row/column with the most cell-like labels.

    The outer region-grouping rows/columns hold only a handful of labels, so the true
    cell-name header row and column win by count.
    """
    n_search = 8
    ncol = max((len(r) for r in grid), default=0)
    header_row = max(
        range(min(n_search, len(grid))),
        key=lambda i: sum(_looks_like_cell(v) for v in grid[i]),
    )
    header_col = max(
        range(min(n_search, ncol)),
        key=lambda j: sum(_looks_like_cell(r[j]) for r in grid if j < len(r)),
    )
    return header_row, header_col


def _parse_matrix(ws, connection_type: str, symmetric: bool) -> tuple[set[str], list[tuple]]:
    """Parse one adjacency sheet -> (cell names, [(pre, post, weight), ...]).

    For a symmetric (gap-junction) matrix each pair is present twice; we keep one edge per
    unordered pair (alphabetically sorted orientation) and drop the mirror.
    """
    grid = [list(row) for row in ws.iter_rows(values_only=True)]
    if not grid:
        return set(), []
    hr, hc = _find_headers(grid)

    colnames = {
        c: grid[hr][c].strip()
        for c in range(len(grid[hr]))
        if c != hc and isinstance(grid[hr][c], str) and grid[hr][c].strip()
    }
    rownames = {
        r: grid[r][hc].strip()
        for r in range(len(grid))
        if r != hr and hc < len(grid[r]) and isinstance(grid[r][hc], str) and grid[r][hc].strip()
    }

    cells = set(colnames.values()) | set(rownames.values())
    edges: list[tuple] = []
    seen: set[frozenset] = set()
    for r, pre in rownames.items():
        row = grid[r]
        for c, post in colnames.items():
            value = row[c] if c < len(row) else None
            if not isinstance(value, (int, float)) or value == 0:
                continue
            if symmetric:
                key = frozenset((pre, post))
                if key in seen:
                    continue
                seen.add(key)
                lo, hi = sorted((pre, post))  # canonical orientation (fresh vars, not pre/post)
                edges.append((lo, hi, float(value)))
            else:
                edges.append((pre, post, float(value)))
    return cells, edges


def read_cook(xlsx_path: Path = DEFAULT_COOK_XLSX) -> CookData:
    """Parse the SI 5 workbook into sex-tagged datasets, connections, and per-sex cell sets."""
    wb = openpyxl.load_workbook(str(xlsx_path), read_only=True, data_only=True)
    try:
        datasets: list[CookDataset] = []
        connections: list[ConnectionRecord] = []
        cells_by_sex: dict[str, set[str]] = {}
        for sex, dataset_id in DATASET_IDS.items():
            sex_cells: set[str] = set()
            for connection_type, template, symmetric in _SHEETS:
                cells, edges = _parse_matrix(
                    wb[template.format(sex=sex)], connection_type, symmetric
                )
                sex_cells |= cells
                for pre, post, weight in edges:
                    connections.append(
                        ConnectionRecord(
                            dataset_id=dataset_id,
                            pre=pre,
                            post=post,
                            connection_type=connection_type,
                            weight=weight,
                            syn=(),
                            ids=None,
                            pre_tid=None,
                            post_tid=None,
                        )
                    )
            cells_by_sex[sex] = sex_cells
            datasets.append(
                CookDataset(
                    id=dataset_id,
                    sex=sex,
                    name=f"Cook et al. 2019 ({sex})",
                    description=(
                        f"Whole-animal {sex} connectome, Cook et al. 2019 (Nature "
                        "571:63-71); SI 5 adjacency matrices, corrected July 2020. "
                        "Weights are EM serial-section counts (synapse number x size)."
                    ),
                )
            )
        return CookData(datasets=datasets, connections=connections, cells_by_sex=cells_by_sex)
    finally:
        wb.close()
