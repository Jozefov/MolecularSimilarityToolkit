"""Resumable streaming normalization, disk-backed deduplication and immutable shards."""

import itertools
import json
import math
import shutil
import time
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from loky import ProcessPoolExecutor

from .adapter import write_shard
from .chemistry import RDKIT_VERSION, prepare
from .config import ChemistryConfig, FingerprintConfig, Source, validate_configs
from .readers import selected_rows
from .resources import effective_workers
from .storage import (
    atomic_json,
    batches,
    check_document,
    check_units,
    connect,
    digest,
    publish_unit,
    read_json,
    receipt_valid,
    seal,
    signature,
    sql_path,
    write_rows,
)
from .versions import (
    FORMAT_VERSION,
    PREPARATION_VERSION,
    __version__,
    check_preparation,
)

STAGE_SCHEMA = pa.schema(
    [
        ("row_id", pa.string()),
        ("source", pa.string()),
        ("source_row", pa.int64()),
        ("original_id", pa.string()),
        ("original_smiles", pa.string()),
        ("split", pa.string()),
        ("metadata", pa.string()),
        ("canonical_smiles", pa.string()),
        ("molecule_id", pa.string()),
        ("fingerprint", pa.binary()),
        ("popcount", pa.int64()),
        ("error", pa.string()),
    ]
)
MOLECULE_SCHEMA = pa.schema(
    [
        ("molecule_id", pa.string()),
        ("canonical_smiles", pa.string()),
        ("fingerprint", pa.binary()),
        ("popcount", pa.int64()),
        ("local_id", pa.int64()),
    ]
)


def source_descriptors(sources):
    values = []
    for index, source in enumerate(sources):
        value = source.to_dict()
        value["name"] = source.name or (
            Path(source.path).name if source.path else f"source-{index}"
        )
        value["checksum"] = (
            digest(source.path) if source.path else signature(value["_records"])
        )
        values.append(value)
    names = [v["name"] for v in values]
    if len(set(names)) != len(names):
        raise ValueError("source names must be unique")
    return values


def _normalize_batch(items, context):
    import os
    import time

    started = time.process_time()
    (
        source_index,
        name,
        id_column,
        smiles_column,
        split_column,
        metadata_columns,
        chemistry,
        fingerprint,
    ) = context
    normalized = []
    for source_row, original in items:
        original_id = original.get(id_column or "id")
        row = {
            "row_id": signature([source_index, source_row]),
            "source": name,
            "source_row": source_row,
            "original_id": str(original_id)
            if original_id not in (None, "")
            else f"row-{source_row}",
            "original_smiles": str(original.get(smiles_column))
            if original.get(smiles_column) is not None
            else "",
            "split": None
            if split_column is None
            else str(original.get(split_column, "")),
            "metadata": json.dumps(
                {key: original.get(key) for key in metadata_columns},
                sort_keys=True,
                default=str,
            ),
        }
        try:
            if original.get("_reader_error"):
                raise ValueError(original["_reader_error"])
            canonical, mol_id, fp, count = prepare(
                row["original_smiles"], chemistry, fingerprint
            )
            row.update(
                canonical_smiles=canonical,
                molecule_id=mol_id,
                fingerprint=fp,
                popcount=count,
                error=None,
            )
        except (ValueError, RuntimeError, OverflowError) as error:
            row.update(
                canonical_smiles=None,
                molecule_id=None,
                fingerprint=None,
                popcount=None,
                error=str(error),
            )
        normalized.append(row)
    return normalized, os.getpid(), time.process_time() - started


def stage_sources(
    sources,
    root: Path,
    chemistry,
    fingerprint,
    *,
    chunk_size=8192,
    workers=4,
    memory_gb=4,
):
    """Reuse completed staging without replay; otherwise use bounded process batches."""
    started = time.perf_counter()
    expected_root = read_json(root / "build.json")["signature"]
    completed_path = root / "staging.json"
    if completed_path.exists():
        completed = read_json(completed_path)
        try:
            check_document(completed)
        except ValueError:
            completed = {"signature": expected_root, "required_units": []}
        if completed["signature"] != expected_root:
            raise ValueError("staging configuration changed")
        if completed["required_units"] and all(
            receipt_valid(root, u["id"], u["signature"])
            for u in completed["required_units"]
        ):
            return [
                root / "units" / u["id"] / "rows.parquet"
                for u in completed["required_units"]
            ]
    paths, units = [], []
    processed = 0
    worker_pids = set()
    cpu_seconds = 0.0
    process_workers = effective_workers(workers, memory_gb, per_worker_gb=0.35)
    executor = None
    try:
        for source_index, source in enumerate(sources):
            name = source.name or (
                Path(source.path).name if source.path else f"source-{source_index}"
            )
            context = (
                source_index,
                name,
                source.id_column,
                source.smiles_column,
                source.split_column,
                source.metadata_columns,
                chemistry,
                fingerprint,
            )
            pending, chunk = [], 0

            def flush(items, number, source_index=source_index, context=context):
                nonlocal executor, processed, cpu_seconds
                unit = f"stage-{source_index:04d}-{number:08d}"
                expected = signature({"root": expected_root, "unit": unit})

                def make(directory):
                    nonlocal executor, cpu_seconds
                    if process_workers > 1 and len(items) >= 1024:
                        if executor is None:
                            executor = ProcessPoolExecutor(
                                max_workers=process_workers,
                            )
                        blocks = [items[i : i + 256] for i in range(0, len(items), 256)]
                        results = executor.map(
                            _normalize_batch, blocks, itertools.repeat(context)
                        )
                    else:
                        results = [_normalize_batch(items, context)]
                    normalized = []
                    for rows, pid, cpu in results:
                        normalized.extend(rows)
                        worker_pids.add(pid)
                        cpu_seconds += cpu
                    write_rows(directory / "rows.parquet", normalized, STAGE_SCHEMA)

                publish_unit(root, unit, expected, make)
                paths.append(root / "units" / unit / "rows.parquet")
                units.append({"id": unit, "signature": expected})
                processed += len(items)
                atomic_json(
                    root / "progress.json",
                    {
                        "state": "staging",
                        "staged_chunks": len(paths),
                        "selected_rows": processed,
                        "source_index": source_index,
                        "process_workers": process_workers,
                        "worker_pids": sorted(worker_pids),
                        "preparation_cpu_seconds": cpu_seconds,
                        "preparation_wall_seconds": time.perf_counter() - started,
                    },
                )

            for source_row, original in selected_rows(source):
                pending.append((source_row, original))
                if len(pending) >= chunk_size:
                    flush(pending, chunk)
                    chunk += 1
                    pending = []
            if pending:
                flush(pending, chunk)
        if not paths:
            unit = "stage-empty"
            expected = signature({"root": expected_root, "unit": unit})
            publish_unit(
                root,
                unit,
                expected,
                lambda directory: write_rows(
                    directory / "rows.parquet", [], STAGE_SCHEMA
                ),
            )
            paths.append(root / "units" / unit / "rows.parquet")
            units.append({"id": unit, "signature": expected})
        atomic_json(
            completed_path,
            seal(
                {
                    "signature": expected_root,
                    "required_units": units,
                    "selected_rows": processed,
                }
            ),
        )
        return paths
    finally:
        if executor is not None:
            executor.shutdown(wait=True)


@dataclass(frozen=True)
class MoleculeIndex:
    """Immutable on-disk index handle. Opening checks completion and chemistry versions."""

    path: Path
    manifest: dict

    @classmethod
    def open(cls, path: str | Path) -> "MoleculeIndex":
        path = Path(path).resolve()
        manifest = read_json(path / "manifest.json")
        check_preparation(manifest)
        check_document(manifest)
        if signature(index_descriptor(manifest)) != manifest.get("signature"):
            raise ValueError("index manifest configuration was changed")
        # Metadata-only opening; scientific readers checksum the units they consume.
        for unit in manifest["required_units"]:
            receipt_path = path / "receipts" / f"{unit['id']}.json"
            if not receipt_path.is_file():
                raise ValueError(f"incomplete index unit: {unit['id']}")
            receipt = read_json(receipt_path)
            if (
                receipt.get("signature") != unit["signature"]
                or not receipt.get("files")
                or not all((path / item["path"]).is_file() for item in receipt["files"])
            ):
                raise ValueError(f"incomplete or incompatible index unit: {unit['id']}")
        if manifest["rdkit_version"] != RDKIT_VERSION:
            raise ValueError(
                "index RDKit version differs from runtime; rebuild the index"
            )
        return cls(path, manifest)

    @property
    def chemistry(self) -> ChemistryConfig:
        return ChemistryConfig.from_dict(self.manifest["chemistry"])

    @property
    def fingerprint(self) -> FingerprintConfig:
        return FingerprintConfig.from_dict(self.manifest["fingerprint"])

    @property
    def size(self) -> int:
        return self.manifest["molecule_count"]

    @property
    def molecules_path(self) -> Path:
        return self.path / "units" / "catalog" / "molecules.parquet"

    @property
    def aliases_path(self) -> Path:
        return self.path / "units" / "catalog" / "aliases.parquet"

    @property
    def invalid_path(self) -> Path:
        return self.path / "units" / "catalog" / "invalid.parquet"

    def molecules(self, batch_size: int = 8192) -> Iterator[pa.RecordBatch]:
        """Iterate Arrow batches without materializing the catalog."""
        return pq.ParquetFile(self.molecules_path).iter_batches(batch_size=batch_size)


def index_descriptor(manifest):
    return {
        key: value
        for key, value in manifest.items()
        if key
        not in {
            "signature",
            "integrity",
            "toolkit_version",
            "molecule_count",
            "alias_count",
            "invalid_count",
            "stage_units",
            "build_inputs",
            "shards",
            "required_units",
            "compact",
        }
    }


DEFAULT_CHEMISTRY = ChemistryConfig()
DEFAULT_FINGERPRINT = FingerprintConfig()


def _prepare_build(
    sources,
    out: str | Path,
    *,
    chemistry: ChemistryConfig = DEFAULT_CHEMISTRY,
    fingerprint: FingerprintConfig = DEFAULT_FINGERPRINT,
    shard_size: int = 1_000_000,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: str | Path | None = None,
    resume: bool = True,
    chunk_size: int = 8192,
) -> Path:
    """Build/resume an immutable index; changed source bytes or policies are rejected.

    ``sources`` is a Source or iterable of Sources. Stage fingerprints are reused
    on resume. Deduplication sorts on disk through DuckDB with a memory limit.
    Each shard has at most one million local integer IDs, irrespective of corpus
    size. Invalid source rows remain in ``invalid.parquet`` with provenance.
    """
    validate_configs(chemistry, fingerprint)
    if (
        not 1 <= shard_size <= 1_000_000
        or memory_gb <= 0
        or workers < 1
        or chunk_size < 1
    ):
        raise ValueError("invalid shard_size, workers, memory_gb or chunk_size")
    sources = [sources] if isinstance(sources, Source) else list(sources)
    root = Path(out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    source_stats = {
        str(source.path): Path(source.path).stat()
        for source in sources
        if source.path is not None
    }
    descriptor = {
        "format_version": FORMAT_VERSION,
        "preparation_version": PREPARATION_VERSION,
        "sources": source_descriptors(sources),
        "chemistry": chemistry.to_dict(),
        "fingerprint": fingerprint.to_dict(),
        "rdkit_version": RDKIT_VERSION,
        "dependency_versions": {
            name: version(name)
            for name in ("FPSim2", "rdkit", "numpy", "pyarrow", "duckdb")
        },
        "shard_size": shard_size,
        "chunk_size": chunk_size,
    }
    expected = signature(descriptor)
    if (root / "manifest.json").exists() and read_json(root / "manifest.json").get(
        "compact"
    ):
        raise ValueError(
            "compact indexes are immutable search artifacts; resume the original build"
        )
    if (root / "build.json").exists():
        if not resume or read_json(root / "build.json")["signature"] != expected:
            raise ValueError(
                "resume refused: inputs or configuration changed, or resume disabled"
            )
    else:
        atomic_json(
            root / "build.json",
            {**descriptor, "signature": expected, "toolkit_version": __version__},
        )
    paths = stage_sources(
        sources,
        root,
        chemistry,
        fingerprint,
        chunk_size=chunk_size,
        workers=workers,
        memory_gb=memory_gb,
    )

    for source in sources:
        if source.path is not None:
            before, after = source_stats[str(source.path)], Path(source.path).stat()
            if (before.st_size, before.st_mtime_ns, before.st_ino) != (
                after.st_size,
                after.st_mtime_ns,
                after.st_ino,
            ):
                raise ValueError(
                    "input changed while staging; use a stable input snapshot"
                )

    def catalog(directory):
        con = connect(memory_gb, scratch or root / "scratch")
        try:
            con.execute(
                "CREATE VIEW staged AS SELECT * FROM read_parquet(["
                + ",".join(sql_path(p) for p in paths)
                + "])"
            )
            con.execute(
                f"COPY (SELECT *,row_number() OVER (ORDER BY molecule_id) AS ordinal FROM (SELECT DISTINCT molecule_id,canonical_smiles,fingerprint,popcount FROM staged WHERE error IS NULL) ORDER BY molecule_id) TO {sql_path(directory / 'molecules.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
            )
            con.execute(
                f"COPY (SELECT row_id,source,source_row,original_id,original_smiles,split,metadata,molecule_id FROM staged WHERE error IS NULL ORDER BY source,source_row) TO {sql_path(directory / 'aliases.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
            )
            con.execute(
                f"COPY (SELECT row_id,source,source_row,original_id,original_smiles,split,metadata,error FROM staged WHERE error IS NOT NULL ORDER BY source,source_row) TO {sql_path(directory / 'invalid.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
            )
        finally:
            con.close()

    publish_unit(root, "catalog", expected, catalog)
    catalog_path = root / "units" / "catalog" / "molecules.parquet"
    count = pq.ParquetFile(catalog_path).metadata.num_rows
    shards = [
        {
            "id": f"shard-{number:08d}",
            "count": min(shard_size, count - number * shard_size),
            "signature": signature([expected, f"shard-{number:08d}"]),
        }
        for number in range(math.ceil(count / shard_size))
    ]
    build_inputs = []
    for number, shard in enumerate(shards):
        unit = f"build-input-{number:08d}"
        unit_signature = signature([expected, unit])

        def partition(directory, number=number):
            con = connect(memory_gb, scratch or root / "scratch")
            try:
                con.execute(
                    f"COPY (SELECT molecule_id,canonical_smiles,fingerprint,popcount,ordinal - {number * shard_size} AS local_id FROM read_parquet({sql_path(catalog_path)}) WHERE ordinal > {number * shard_size} AND ordinal <= {(number + 1) * shard_size} ORDER BY ordinal) TO {sql_path(directory / 'molecules.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
                )
            finally:
                con.close()

        publish_unit(root, unit, unit_signature, partition)
        build_inputs.append({"id": unit, "signature": unit_signature})
    manifest = {
        **descriptor,
        "toolkit_version": read_json(root / "build.json")["toolkit_version"],
        "compact": False,
        "build_inputs": build_inputs,
        "required_units": [{"id": "catalog", "signature": expected}]
        + [{"id": sh["id"], "signature": sh["signature"]} for sh in shards],
        "stage_units": [
            {
                "id": p.parent.name,
                "signature": read_json(root / "receipts" / f"{p.parent.name}.json")[
                    "signature"
                ],
            }
            for p in paths
        ],
        "signature": expected,
        "molecule_count": count,
        "alias_count": pq.ParquetFile(
            root / "units/catalog/aliases.parquet"
        ).metadata.num_rows,
        "invalid_count": pq.ParquetFile(
            root / "units/catalog/invalid.parquet"
        ).metadata.num_rows,
        "shards": shards,
    }
    manifest = seal(manifest)
    plan = {
        "kind": "build",
        "required_units": manifest["stage_units"]
        + [{"id": "catalog", "signature": expected}]
        + build_inputs
        + [{"id": sh["id"], "signature": sh["signature"]} for sh in shards],
        "signature": expected,
        "manifest": manifest,
        "memory_gb": memory_gb,
        "scratch": str(scratch or root / "scratch"),
        "tasks": [
            {"id": number, "shard": shard["id"]} for number, shard in enumerate(shards)
        ],
    }
    plan["tasks_signature"] = signature(plan["tasks"])
    atomic_json(root / "build-plan.json", seal(plan))
    return root / "build-plan.json"


def build_index(
    sources: Source | Iterable[Source],
    out: str | Path,
    *,
    chemistry: ChemistryConfig = DEFAULT_CHEMISTRY,
    fingerprint: FingerprintConfig = DEFAULT_FINGERPRINT,
    shard_size: int = 1_000_000,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: str | Path | None = None,
    resume: bool = True,
    chunk_size: int = 8192,
) -> MoleculeIndex:
    """Build/resume a structure-deduplicated index from streaming sources.

    Defaults keep charges, salts, isotope labels and tautomer forms while
    removing stereochemistry. Invalid input molecules are quarantined with
    provenance. Successful chunks reuse stored fingerprints on resume; changed
    input bytes, chemistry, fingerprints or shard layout are refused.

    The one staging coordinator limits DuckDB memory and spills to ``scratch``.
    Chemistry uses bounded query threads; shard writers execute sequentially.
    For distributed shard builds use :func:`plan_build`, ``run_task`` and
    ``merge``. Each HDF5 shard uses IDs in [1, 1,000,000], independently of the
    total corpus size. Source aliases retain duplicate IDs and memberships.
    """
    path = _prepare_build(
        sources,
        out,
        chemistry=chemistry,
        fingerprint=fingerprint,
        shard_size=shard_size,
        workers=workers,
        memory_gb=memory_gb,
        scratch=scratch,
        resume=resume,
        chunk_size=chunk_size,
    )
    context = _build_plan(path)
    for task in context[1]["tasks"]:
        run_build_task(path, task["id"], _context=context)
    return merge_build(path)


def plan_build(
    sources: Source | Iterable[Source],
    out: str | Path,
    *,
    chemistry: ChemistryConfig = DEFAULT_CHEMISTRY,
    fingerprint: FingerprintConfig = DEFAULT_FINGERPRINT,
    shard_size: int = 1_000_000,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: str | Path | None = None,
    resume: bool = True,
    chunk_size: int = 8192,
) -> Path:
    """Coordinate streaming staging and dedup, then plan independent HDF5 shards.

    Only this preparation step writes staging/catalog tables. Distributed tasks
    write separate immutable shard directories and never share HDF5 writers.
    """
    return _prepare_build(
        sources,
        out,
        chemistry=chemistry,
        fingerprint=fingerprint,
        shard_size=shard_size,
        workers=workers,
        memory_gb=memory_gb,
        scratch=scratch,
        resume=resume,
        chunk_size=chunk_size,
    )


def _build_plan(plan_path):
    path = Path(plan_path).resolve()
    plan = read_json(path)
    check_document(plan)
    check_document(plan["manifest"])
    check_preparation(plan["manifest"])
    if (
        signature(index_descriptor(plan["manifest"])) != plan["signature"]
        or signature(plan["tasks"]) != plan["tasks_signature"]
    ):
        raise ValueError("build plan configuration or tasks changed")
    if [task["id"] for task in plan["tasks"]] != list(range(len(plan["tasks"]))):
        raise ValueError("invalid or duplicate build task IDs")
    return path.parent, plan


def run_build_task(plan_path: str | Path, task_id: int, *, _context=None) -> Path:
    """Execute one shard build from cached fingerprints and publish its receipt."""
    root, plan = _context or _build_plan(plan_path)
    if task_id < 0 or task_id >= len(plan["tasks"]):
        raise ValueError("task_id outside build plan")
    manifest = plan["manifest"]
    shard = manifest["shards"][task_id]
    if plan["tasks"][task_id]["shard"] != shard["id"]:
        raise ValueError("build task does not match shard")

    dependency = manifest["build_inputs"][task_id]
    check_units(root, [dependency])

    def make(directory):
        shutil.copyfile(
            root / "units" / dependency["id"] / "molecules.parquet",
            directory / "molecules.parquet",
        )
        write_shard(
            directory / "fingerprints.h5",
            batches(directory / "molecules.parquet"),
            FingerprintConfig.from_dict(manifest["fingerprint"]),
        )

    publish_unit(root, shard["id"], shard["signature"], make)
    return root / "receipts" / f"{shard['id']}.json"


def merge_build(plan_path: str | Path) -> MoleculeIndex:
    """Publish an index only after all immutable shard receipts validate."""
    root, plan = _build_plan(plan_path)
    check_units(root, plan["required_units"], exact=True)
    expected = {shard["id"] for shard in plan["manifest"]["shards"]}
    existing = {path.stem for path in (root / "receipts").glob("shard-*.json")}
    if existing - expected:
        raise ValueError("unexpected build shard receipts")
    for shard in plan["manifest"]["shards"]:
        if not receipt_valid(root, shard["id"], shard["signature"]):
            raise ValueError(f"incomplete or corrupt build shard {shard['id']}")
    atomic_json(root / "manifest.json", plan["manifest"])
    return MoleculeIndex.open(root)


def export_index(index: MoleculeIndex | str | Path, out: str | Path) -> MoleculeIndex:
    """Atomically copy a verified index without staging/build inputs or raw sources.

    The compact index supports search/filter/verify, but not rebuilding. Original
    source descriptors and all aliases remain as provenance. Never modifies inputs.
    """
    import os
    import tempfile

    index = MoleculeIndex.open(
        index.path if isinstance(index, MoleculeIndex) else index
    )
    target = Path(out).resolve()
    if target.exists():
        raise ValueError("export destination must not exist")
    check_units(index.path, index.manifest["required_units"])
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".compact-", dir=target.parent))
    try:
        for unit in index.manifest["required_units"]:
            shutil.copytree(
                index.path / "units" / unit["id"], temporary / "units" / unit["id"]
            )
            (temporary / "receipts").mkdir(exist_ok=True)
            shutil.copyfile(
                index.path / "receipts" / f"{unit['id']}.json",
                temporary / "receipts" / f"{unit['id']}.json",
            )
        manifest = seal(
            {**index.manifest, "stage_units": [], "build_inputs": [], "compact": True}
        )
        atomic_json(temporary / "manifest.json", manifest)
        os.rename(temporary, target)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)
    return MoleculeIndex.open(target)
