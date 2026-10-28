"""Resumable streaming normalization, disk-backed deduplication and immutable shards."""

import json
import math
from collections.abc import Iterable, Iterator
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .adapter import write_shard
from .chemistry import RDKIT_VERSION, prepare
from .config import ChemistryConfig, FingerprintConfig, Source, validate_configs
from .readers import selected_rows
from .resources import effective_workers
from .storage import (
    atomic_json,
    batches,
    connect,
    digest,
    publish_unit,
    read_json,
    receipt_valid,
    signature,
    sql_path,
    write_rows,
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
    paths = []
    for source_index, source in enumerate(sources):
        name = source.name or (
            Path(source.path).name if source.path else f"source-{source_index}"
        )
        pending = []
        chunk = 0

        def flush(items, number, source_index=source_index, name=name, source=source):
            unit = f"stage-{source_index:04d}-{number:08d}"
            expected = signature(
                {"root": read_json(root / "build.json")["signature"], "unit": unit}
            )

            def make(directory):
                def normalize(item):
                    source_row, original = item
                    row = {
                        "row_id": signature([source_index, source_row]),
                        "source": name,
                        "source_row": source_row,
                        "original_id": str(original.get(source.id_column, ""))
                        if original.get(source.id_column) is not None
                        else "",
                        "original_smiles": str(original.get(source.smiles_column, ""))
                        if original.get(source.smiles_column) is not None
                        else "",
                        "split": None
                        if source.split_column is None
                        else str(original.get(source.split_column, "")),
                        "metadata": json.dumps(
                            {key: original.get(key) for key in source.metadata_columns},
                            sort_keys=True,
                            default=str,
                        ),
                    }
                    try:
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
                    return row

                with ThreadPoolExecutor(
                    max_workers=effective_workers(workers, memory_gb)
                ) as pool:
                    normalized = list(pool.map(normalize, items))
                write_rows(directory / "rows.parquet", normalized, STAGE_SCHEMA)

            publish_unit(root, unit, expected, make)
            paths.append(root / "units" / unit / "rows.parquet")

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
        expected = signature(
            {"root": read_json(root / "build.json")["signature"], "unit": unit}
        )
        publish_unit(
            root,
            unit,
            expected,
            lambda directory: write_rows(directory / "rows.parquet", [], STAGE_SCHEMA),
        )
        paths.append(root / "units" / unit / "rows.parquet")
    return paths


@dataclass(frozen=True)
class MoleculeIndex:
    """Immutable on-disk index handle. Opening checks completion and chemistry versions."""

    path: Path
    manifest: dict

    @classmethod
    def open(cls, path: str | Path) -> "MoleculeIndex":
        path = Path(path).resolve()
        manifest = read_json(path / "manifest.json")
        if manifest.get("format_version") != 1:
            raise ValueError("unsupported index format version")
        descriptor = {
            key: value
            for key, value in manifest.items()
            if key
            not in {
                "signature",
                "molecule_count",
                "alias_count",
                "invalid_count",
                "stage_units",
                "shards",
            }
        }
        if signature(descriptor) != manifest.get("signature"):
            raise ValueError("index manifest configuration was changed")
        plan = read_json(path / "build-plan.json")
        if (
            plan["manifest"] != manifest
            or signature(plan["tasks"]) != plan["tasks_signature"]
        ):
            raise ValueError(
                "index manifest is incomplete or incompatible with its build plan"
            )
        required = (
            [("catalog", manifest["signature"])]
            + [(unit["id"], unit["signature"]) for unit in manifest["stage_units"]]
            + [(shard["id"], shard["signature"]) for shard in manifest["shards"]]
        )
        for unit, expected in required:
            receipt_path = path / "receipts" / f"{unit}.json"
            if not receipt_path.is_file():
                raise ValueError(f"incomplete index unit: {unit}")
            receipt = read_json(receipt_path)
            if (
                receipt["signature"] != expected
                or not receipt["files"]
                or not all((path / item["path"]).is_file() for item in receipt["files"])
            ):
                raise ValueError(f"incomplete or incompatible index unit: {unit}")
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
        "format_version": 1,
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
    if (root / "build.json").exists():
        if not resume or read_json(root / "build.json")["signature"] != expected:
            raise ValueError(
                "resume refused: inputs or configuration changed, or resume disabled"
            )
    else:
        atomic_json(root / "build.json", {**descriptor, "signature": expected})
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
    manifest = {
        **descriptor,
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
    plan = {
        "kind": "build",
        "signature": expected,
        "manifest": manifest,
        "memory_gb": memory_gb,
        "scratch": str(scratch or root / "scratch"),
        "tasks": [
            {"id": number, "shard": shard["id"]} for number, shard in enumerate(shards)
        ],
    }
    plan["tasks_signature"] = signature(plan["tasks"])
    atomic_json(root / "build-plan.json", plan)
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
    descriptor = {
        key: value
        for key, value in plan["manifest"].items()
        if key
        not in {
            "signature",
            "molecule_count",
            "alias_count",
            "invalid_count",
            "stage_units",
            "shards",
        }
    }
    if (
        signature(descriptor) != plan["signature"]
        or signature(plan["tasks"]) != plan["tasks_signature"]
    ):
        raise ValueError("build plan configuration or tasks changed")
    if [task["id"] for task in plan["tasks"]] != list(range(len(plan["tasks"]))):
        raise ValueError("invalid or duplicate build task IDs")
    if not receipt_valid(path.parent, "catalog", plan["signature"]):
        raise ValueError("build catalog is incomplete or corrupt")
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

    def make(directory):
        con = connect(plan["memory_gb"], plan["scratch"])
        try:
            catalog_path = root / "units/catalog/molecules.parquet"
            shard_size = manifest["shard_size"]
            con.execute(
                f"COPY (SELECT molecule_id,canonical_smiles,fingerprint,popcount,ordinal - {task_id * shard_size} AS local_id FROM read_parquet({sql_path(catalog_path)}) WHERE ordinal > {task_id * shard_size} AND ordinal <= {(task_id + 1) * shard_size} ORDER BY ordinal) TO {sql_path(directory / 'molecules.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
            )
            write_shard(
                directory / "fingerprints.h5",
                batches(directory / "molecules.parquet"),
                FingerprintConfig.from_dict(manifest["fingerprint"]),
            )
        finally:
            con.close()

    publish_unit(root, shard["id"], shard["signature"], make)
    return root / "receipts" / f"{shard['id']}.json"


def merge_build(plan_path: str | Path) -> MoleculeIndex:
    """Publish an index only after all immutable shard receipts validate."""
    root, plan = _build_plan(plan_path)
    expected = {shard["id"] for shard in plan["manifest"]["shards"]}
    existing = {path.stem for path in (root / "receipts").glob("shard-*.json")}
    if existing - expected:
        raise ValueError("unexpected build shard receipts")
    for shard in plan["manifest"]["shards"]:
        if not receipt_valid(root, shard["id"], shard["signature"]):
            raise ValueError(f"incomplete or corrupt build shard {shard['id']}")
    atomic_json(root / "manifest.json", plan["manifest"])
    return MoleculeIndex.open(root)
