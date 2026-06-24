# Contributing

Thanks for your interest in the C. elegans Connectome Knowledge Graph. This document
covers how to set up a dev environment, the conventions the project follows, and what CI
expects before a change can merge.

## Development setup

This project uses [uv](https://docs.astral.sh/uv/) for environment and dependency
management. The lockfile (`uv.lock`) is committed and authoritative.

```bash
uv sync --extra dev      # create .venv and install runtime + dev dependencies
uv run cckg --help       # pipeline CLI
```

## Project conventions

- **The LinkML schema is the single source of truth.** RDF/OWL and the neuron-graph-shaped
  JSON projection are *generated* from it — never hand-edit generated artifacts in
  `outputs/` (they are gitignored and rebuilt by the pipeline).
- **Source data is pinned, not vendored.** Anything under `data/` is a reproducibility
  snapshot with provenance and checksums in its `MANIFEST.md`. neuron-graph and WBBT remain
  external upstreams; don't fork or re-import them into this repo.
- See [`docs/PLANNING.md`](docs/PLANNING.md) for the design, roadmap, and which phase each
  package directory belongs to.

## Before you open a pull request

CI runs ruff and pytest on Python 3.10, 3.11, and 3.12. Run the same checks locally:

```bash
uv run ruff check .          # lint
uv run ruff format .         # apply formatting (CI runs --check)
uv run pytest                # tests
```

All three must pass. `ruff format --check` is enforced in CI, so format your code before
pushing.

## Pull request guidelines

- Branch off `main` and keep PRs focused on a single concern.
- Write a clear description of *what* changed and *why*.
- Add or update tests for behavior changes. New code in `ingest/`, `match/`, `build/`, and
  `export/` should be covered by tests under `tests/`.
- Keep commits scoped and their messages in the imperative mood (e.g. "Add WBBT matcher").
- Make sure CI is green before requesting review.

## Reporting issues

Open a GitHub issue with enough detail to reproduce: what you ran, what you expected, and
what happened. For data discrepancies, include the relevant `MANIFEST.md` provenance so the
source snapshot is unambiguous.

## License

By contributing, you agree that your contributions are licensed under the
[MIT License](LICENSE).
