"""Command-line entry point for the pipeline.

Stages are stubs in Phase 0; each is filled in over the roadmap (see docs/PLANNING.md).
"""

import click


@click.group()
@click.version_option()
def main() -> None:
    """cckg — C. elegans connectome knowledge graph pipeline."""


@main.command()
def ingest() -> None:
    """Read pinned neuron-graph files into normalized records. [Phase 2]"""
    raise click.ClickException("not yet implemented (Phase 2)")


@main.command()
def match() -> None:
    """Resolve cell names to WBbt anatomy URIs; emit the match report. [Phase 2]"""
    raise click.ClickException("not yet implemented (Phase 2)")


@main.command()
def build() -> None:
    """Assemble LinkML data (cells, connections, datasets, evidence). [Phase 3]"""
    raise click.ClickException("not yet implemented (Phase 3)")


@main.command()
def export() -> None:
    """Serialize RDF/OWL and the neuron-graph JSON projection. [Phase 3]"""
    raise click.ClickException("not yet implemented (Phase 3)")


if __name__ == "__main__":
    main()
