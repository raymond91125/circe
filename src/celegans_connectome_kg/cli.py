"""Command-line entry point for the pipeline.

Stages are stubs in Phase 0; each is filled in over the roadmap (see docs/PLANNING.md).
"""

from collections import Counter
from pathlib import Path

import click

DEFAULT_NEURON_GRAPH_DIR = Path("data/neuron-graph")
DEFAULT_WBBT_JSON = Path("data/wbbt/wbbt.json")
DEFAULT_OUTPUT_DIR = Path("outputs")


@click.group()
@click.version_option()
def main() -> None:
    """cckg — C. elegans connectome knowledge graph pipeline."""


@main.command()
@click.option(
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_NEURON_GRAPH_DIR,
    show_default=True,
    help="Pinned neuron-graph snapshot directory.",
)
def ingest(data_dir: Path) -> None:
    """Read pinned neuron-graph files into normalized records. [Phase 2]"""
    from celegans_connectome_kg.ingest.neuron_graph import load_neuron_graph

    data = load_neuron_graph(data_dir)
    by_type = Counter(c.connection_type for c in data.connections)
    by_dataset = Counter(c.dataset_id for c in data.connections)
    click.echo(f"cells:       {len(data.cells)}")
    click.echo(f"datasets:    {len(data.datasets)}")
    click.echo(f"connections: {len(data.connections)}")
    for t in ("chemical", "gap_junction", "functional"):
        click.echo(f"  {t:13} {by_type.get(t, 0)}")
    click.echo(f"datasets with connections: {len(by_dataset)}")


@main.command()
@click.option(
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_NEURON_GRAPH_DIR,
    show_default=True,
    help="Pinned neuron-graph snapshot directory.",
)
@click.option(
    "--wbbt",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_WBBT_JSON,
    show_default=True,
    help="Pinned WBBT OBO-graph JSON.",
)
@click.option(
    "--out-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT_DIR,
    show_default=True,
    help="Where to write the match report and work-list.",
)
def match(data_dir: Path, wbbt: Path, out_dir: Path) -> None:
    """Resolve cell names to WBbt anatomy URIs; emit the match report. [Phase 2]"""
    from celegans_connectome_kg.ingest.neuron_graph import read_cells
    from celegans_connectome_kg.match.matcher import (
        match_cells,
        summarize,
        write_report_csv,
        write_worklist_csv,
    )
    from celegans_connectome_kg.match.wbbt import WBBTIndex

    index = WBBTIndex.from_obograph(wbbt)
    cells = read_cells(data_dir / "neurons.json")
    matches = match_cells(cells, index)
    counts = summarize(matches)

    report_path = out_dir / "match_report.csv"
    worklist_path = out_dir / "match_worklist.csv"
    write_report_csv(matches, report_path)
    write_worklist_csv(matches, cells, worklist_path)

    click.echo(f"WBBT terms indexed: {len(index.terms)}")
    click.echo(f"cells:     {len(cells)}")
    for status in ("matched", "ambiguous", "unmatched"):
        click.echo(f"  {status:10} {counts.get(status, 0)}")
    click.echo(f"report:   {report_path}")
    click.echo(
        f"worklist: {worklist_path} ({counts.get('ambiguous', 0) + counts.get('unmatched', 0)} rows)"
    )


@main.command()
@click.option(
    "--data-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=DEFAULT_NEURON_GRAPH_DIR,
    show_default=True,
    help="Pinned neuron-graph snapshot directory.",
)
@click.option(
    "--wbbt",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_WBBT_JSON,
    show_default=True,
    help="Pinned WBBT OBO-graph JSON.",
)
@click.option(
    "--out-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT_DIR,
    show_default=True,
    help="Where to write the assembled data.",
)
def build(data_dir: Path, wbbt: Path, out_dir: Path) -> None:
    """Assemble LinkML data (cells, connections, datasets, evidence). [Phase 3]"""
    from celegans_connectome_kg.build.assemble import assemble
    from celegans_connectome_kg.export.rdf import write_json

    connectome, stats = assemble(data_dir, wbbt)
    out_path = out_dir / "connectome.json"
    write_json(connectome, out_path)

    click.echo(f"cells:       {stats.cells} ({stats.cells_with_anatomy} with WBBT anatomy)")
    click.echo(f"datasets:    {stats.datasets}")
    click.echo(f"connections: {stats.connections}")
    for t in ("chemical", "gap_junction", "functional"):
        click.echo(f"  {t:13} {stats.connections_by_type.get(t, 0)}")
    if stats.unknown_connection_cells:
        click.echo(
            f"note: {stats.unknown_connection_cells} connection endpoints not in the cell list"
        )
    click.echo(f"wrote: {out_path}")


@main.command()
@click.option(
    "--in",
    "in_path",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT_DIR / "connectome.json",
    show_default=True,
    help="LinkML-native JSON produced by `cckg build`.",
)
@click.option(
    "--out-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=DEFAULT_OUTPUT_DIR,
    show_default=True,
    help="Where to write the serialized graph.",
)
def export(in_path: Path, out_dir: Path) -> None:
    """Serialize RDF/OWL and the neuron-graph JSON projection. [Phase 3]"""
    from celegans_connectome_kg.export.rdf import load_json, write_turtle

    connectome = load_json(in_path)
    ttl_path = out_dir / "connectome.ttl"
    write_turtle(connectome, ttl_path)
    click.echo(f"wrote: {ttl_path}")


if __name__ == "__main__":
    main()
