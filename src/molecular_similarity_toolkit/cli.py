"""Scheduler-independent command line entry point."""

import json
from pathlib import Path

import typer

from .config import ChemistryConfig, FingerprintConfig, Source
from .index import MoleculeIndex, build_index, export_index, plan_build
from .inspection import status as inspect_status
from .inspection import verify as inspect_verify
from .operations import filter_pool, plan_filter, plan_search, run_task
from .operations import merge as merge_api
from .operations import search as search_api
from .results import Result

app = typer.Typer(
    help="Exact molecular similarity with immutable local/HPC work units.",
    no_args_is_help=True,
    pretty_exceptions_enable=False,
)
plan_app = typer.Typer(
    help="Prepare immutable work units for any scheduler.", no_args_is_help=True
)
app.add_typer(plan_app, name="plan")


def _sources(path: Path):
    if path.suffix.lower() == ".json":
        value = json.loads(path.read_text())
        entries = value.get("sources", value) if isinstance(value, dict) else value
        if isinstance(entries, dict):
            entries = [entries]
        result = []
        for entry in entries:
            entry = dict(entry)
            if entry.get("path"):
                entry["path"] = str((path.parent / entry["path"]).resolve())
            result.append(Source.from_dict(entry))
        return result
    if path.is_dir():
        return MoleculeIndex.open(path)
    return Source(path)


def _config(path, cls):
    return cls.from_dict(json.loads(path.read_text())) if path else cls()


def _build(
    sources,
    out,
    chemistry,
    fingerprint,
    shard_size,
    workers,
    memory_gb,
    scratch,
    planning,
):
    result = (plan_build if planning else build_index)(
        _sources(sources),
        out,
        chemistry=_config(chemistry, ChemistryConfig),
        fingerprint=_config(fingerprint, FingerprintConfig),
        shard_size=shard_size,
        workers=workers,
        memory_gb=memory_gb,
        scratch=scratch,
    )
    typer.echo(str(result if planning else result.path))


@app.command()
def build(
    sources: Path,
    out: Path,
    chemistry: Path | None = None,
    fingerprint: Path | None = None,
    shard_size: int = 1_000_000,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: Path | None = None,
):
    """Build/resume sources (a data path or JSON Source configuration)."""
    _build(
        sources,
        out,
        chemistry,
        fingerprint,
        shard_size,
        workers,
        memory_gb,
        scratch,
        False,
    )


@plan_app.command("build")
def build_plan(
    sources: Path,
    out: Path,
    chemistry: Path | None = None,
    fingerprint: Path | None = None,
    shard_size: int = 1_000_000,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: Path | None = None,
):
    """Stage/deduplicate once, then plan independent HDF5 shard builds."""
    _build(
        sources,
        out,
        chemistry,
        fingerprint,
        shard_size,
        workers,
        memory_gb,
        scratch,
        True,
    )


@app.command()
def search(
    index: Path,
    queries: Path,
    out: Path,
    k: int = 10,
    exclude_identity: bool = False,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: Path | None = None,
):
    """Search every shard and atomically merge exact top-k results."""
    result = search_api(
        index,
        _sources(queries),
        k,
        out,
        exclude_identity=exclude_identity,
        workers=workers,
        memory_gb=memory_gb,
        scratch=scratch,
    )
    typer.echo(json.dumps(result.summary, indent=2))


@plan_app.command("search")
def search_plan(
    index: Path,
    queries: Path,
    out: Path,
    k: int = 10,
    exclude_identity: bool = False,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: Path | None = None,
    query_chunk_size: int = 256,
):
    """Prepare query-chunk × pool-shard tasks."""
    typer.echo(
        str(
            plan_search(
                index,
                _sources(queries),
                k,
                out,
                exclude_identity=exclude_identity,
                workers=workers,
                memory_gb=memory_gb,
                scratch=scratch,
                query_chunk_size=query_chunk_size,
            )
        )
    )


@app.command("filter")
def filter_command(
    index: Path,
    protected: Path,
    out: Path,
    threshold: str = "0.8",
    comparison: str = ">",
    invalid_protected: str = "error",
    workers: int = 4,
    memory_gb: float = 4,
    scratch: Path | None = None,
):
    """Filter pool molecules against protected molecules using exact boundaries."""
    result = filter_pool(
        index,
        _sources(protected),
        threshold,
        comparison,
        out,
        invalid_protected=invalid_protected,
        workers=workers,
        memory_gb=memory_gb,
        scratch=scratch,
    )
    typer.echo(json.dumps(result.summary, indent=2))


@plan_app.command("filter")
def filter_plan(
    index: Path,
    protected: Path,
    out: Path,
    threshold: str = "0.8",
    comparison: str = ">",
    invalid_protected: str = "error",
    workers: int = 4,
    memory_gb: float = 4,
    scratch: Path | None = None,
):
    """Prepare independent pool-shard filtering tasks."""
    typer.echo(
        str(
            plan_filter(
                index,
                _sources(protected),
                out,
                threshold=threshold,
                comparison=comparison,
                invalid_protected=invalid_protected,
                workers=workers,
                memory_gb=memory_gb,
                scratch=scratch,
            )
        )
    )


@app.command()
def task(plan: Path, task_id: int, workers: int | None = None):
    """Run one numeric task ID; safe to retry missing/corrupt units."""
    typer.echo(str(run_task(plan, task_id, workers=workers)))


@app.command()
def merge(plan: Path):
    """Refuse incomplete/corrupt/incompatible task outputs, then publish results."""
    result = merge_api(plan)
    typer.echo(str(result.path))


@app.command()
def status(path: Path):
    """Report completion (use verify for checksummed integrity)."""
    typer.echo(json.dumps(inspect_status(path), indent=2))


@app.command()
def verify(path: Path):
    """Read and checksum every immutable output; exit nonzero on failure."""
    result = inspect_verify(path)
    typer.echo(json.dumps(result, indent=2))
    if not result["ok"]:
        raise typer.Exit(1)


@app.command("compact")
def compact(index: Path, out: Path):
    """Export a verified portable index without preparation/build artifacts."""
    typer.echo(str(export_index(index, out).path))


@app.command("export")
def export_command(run: Path, out: Path, table: str | None = None):
    """Stream a completed result table to CSV/TSV."""
    typer.echo(str(Result.open(run).export(out, table=table)))


if __name__ == "__main__":
    app()
