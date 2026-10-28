"""Exact all-shard search and leakage filtering with immutable scheduler tasks."""

import heapq
import itertools
import math
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .adapter import ShardEngine
from .chemistry import intersection_union
from .config import Source
from .index import MoleculeIndex, _build_plan, build_index, merge_build, run_build_task
from .resources import effective_workers
from .results import Result
from .storage import (
    atomic_json,
    batches,
    check_document,
    check_units,
    connect,
    failpoint,
    publish_unit,
    read_json,
    receipt_valid,
    seal,
    signature,
    sql_path,
    write_rows,
)
from .versions import OPERATION_VERSION, PREPARATION_VERSION

MATCH_SCHEMA = pa.schema(
    [
        ("query_id", pa.string()),
        ("query_original_id", pa.string()),
        ("query_molecule_id", pa.string()),
        ("molecule_id", pa.string()),
        ("rank", pa.int64()),
        ("similarity", pa.float64()),
        ("intersection", pa.int64()),
        ("union", pa.int64()),
    ]
)
HIT_SCHEMA = pa.schema(
    [
        ("molecule_id", pa.string()),
        ("protected_id", pa.string()),
        ("similarity", pa.float64()),
        ("intersection", pa.int64()),
        ("union", pa.int64()),
    ]
)
QUERY_SCHEMA = pa.schema(
    [
        ("query_id", pa.string()),
        ("query_original_id", pa.string()),
        ("molecule_id", pa.string()),
        ("fingerprint", pa.binary()),
    ]
)


def _index(value, root: Path, reference: MoleculeIndex, workers, memory_gb, scratch):
    if isinstance(value, MoleculeIndex):
        result = MoleculeIndex.open(value.path)
    elif isinstance(value, (str, Path)) and (Path(value) / "manifest.json").exists():
        result = MoleculeIndex.open(value)
    else:
        if isinstance(value, (str, Path)):
            value = Source(value)
        result = build_index(
            value,
            root,
            chemistry=reference.chemistry,
            fingerprint=reference.fingerprint,
            workers=workers,
            memory_gb=memory_gb,
            scratch=scratch,
        )
    if not receipt_valid(result.path, "catalog", result.manifest["signature"]):
        raise ValueError("query/protected catalog is corrupt; resume its build")
    for key in ("chemistry", "fingerprint", "rdkit_version", "preparation_version"):
        if result.manifest[key] != reference.manifest[key]:
            raise ValueError(f"incompatible indexes: {key}")
    return result


def _reference(value):
    return MoleculeIndex.open(value.path if isinstance(value, MoleculeIndex) else value)


def plan_search(
    index: MoleculeIndex | str | Path,
    queries: Source | list[Source] | MoleculeIndex | str | Path,
    k: int,
    out: str | Path,
    *,
    exclude_identity: bool = False,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: str | Path | None = None,
    query_chunk_size: int = 256,
) -> Path:
    """Stage query rows and create immutable query-chunk × shard task manifest.

    Query IDs are generated source-row IDs; ``query_original_id`` preserves the
    provided ID even when original IDs repeat. Identical query rows each get
    their own result rows. Index-to-index searches require identical policies.
    """
    if k < 1 or query_chunk_size < 1:
        raise ValueError("k and query_chunk_size must be positive")
    reference = _reference(index)
    if not receipt_valid(reference.path, "catalog", reference.manifest["signature"]):
        raise ValueError("pool catalog is corrupt; resume its build")
    root = Path(out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    query_index = _index(
        queries, root / "queries", reference, workers, memory_gb, scratch
    )
    descriptor = {
        "kind": "search",
        "operation_version": OPERATION_VERSION,
        "preparation_version": PREPARATION_VERSION,
        "index": str(reference.path),
        "index_signature": reference.manifest["signature"],
        "queries": str(query_index.path),
        "query_signature": query_index.manifest["signature"],
        "k": k,
        "exclude_identity": exclude_identity,
        "query_chunk_size": query_chunk_size,
    }
    return _plan(root, descriptor, reference, query_index, workers, memory_gb, scratch)


def plan_filter(
    index: MoleculeIndex | str | Path,
    protected: Source | list[Source] | MoleculeIndex | str | Path,
    out: str | Path,
    *,
    threshold: float | str = 0.8,
    comparison: str = ">",
    invalid_protected: str = "error",
    workers: int = 4,
    memory_gb: float = 4,
    scratch: str | Path | None = None,
) -> Path:
    """Create one pool-shard task per unit; protected rows must contain molecules.

    The decimal threshold is compared as an exact rational number. Rejection
    witnesses use the lowest protected molecule ID satisfying the comparison,
    irrespective of which witness has the highest similarity.
    """
    value = Fraction(str(threshold))
    if not 0 <= value <= 1 or comparison not in (">", ">="):
        raise ValueError("threshold must be in [0,1]; comparison must be > or >=")
    if invalid_protected not in ("error", "quarantine"):
        raise ValueError("invalid_protected must be error or quarantine")
    reference = _reference(index)
    if not receipt_valid(reference.path, "catalog", reference.manifest["signature"]):
        raise ValueError("pool catalog is corrupt; resume its build")
    root = Path(out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    protected_index = _index(
        protected, root / "protected", reference, workers, memory_gb, scratch
    )
    if not protected_index.size:
        raise ValueError("protected set is empty")
    if protected_index.manifest["invalid_count"] and invalid_protected == "error":
        raise ValueError(
            "protected set contains invalid molecules; inspect its invalid.parquet or explicitly quarantine"
        )
    descriptor = {
        "kind": "filter",
        "operation_version": OPERATION_VERSION,
        "preparation_version": PREPARATION_VERSION,
        "index": str(reference.path),
        "index_signature": reference.manifest["signature"],
        "queries": str(protected_index.path),
        "query_signature": protected_index.manifest["signature"],
        "threshold": str(threshold),
        "comparison": comparison,
        "invalid_protected": invalid_protected,
    }
    return _plan(
        root, descriptor, reference, protected_index, workers, memory_gb, scratch
    )


def _plan(root, descriptor, reference, query_index, workers, memory_gb, scratch):
    expected = signature(descriptor)
    path = root / "plan.json"
    if path.exists() and read_json(path)["signature"] != expected:
        raise ValueError("resume refused: operation inputs or configuration changed")

    def query_rows(directory):
        con = connect(memory_gb, scratch or root / "scratch")
        try:
            if descriptor["kind"] == "search":
                con.execute(
                    f"COPY (SELECT a.row_id AS query_id,a.original_id AS query_original_id,m.molecule_id,m.fingerprint FROM read_parquet({sql_path(query_index.aliases_path)}) a JOIN read_parquet({sql_path(query_index.molecules_path)}) m USING(molecule_id) ORDER BY a.row_id) TO {sql_path(directory / 'queries.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
                )
            else:
                con.execute(
                    f"COPY (SELECT molecule_id AS query_id,molecule_id AS query_original_id,molecule_id,fingerprint FROM read_parquet({sql_path(query_index.molecules_path)}) ORDER BY molecule_id) TO {sql_path(directory / 'queries.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
                )
        finally:
            con.close()

    publish_unit(root, "queries", expected, query_rows)
    query_path = root / "units/queries/queries.parquet"
    query_count = pq.ParquetFile(query_path).metadata.num_rows
    tasks = []
    if descriptor["kind"] == "search":
        chunks = math.ceil(query_count / descriptor["query_chunk_size"])
        query_batches = pq.ParquetFile(query_path).iter_batches(
            batch_size=descriptor["query_chunk_size"]
        )
        for chunk in range(chunks):
            batch = next(query_batches)
            chunk_unit = f"query-chunk-{chunk:08d}"
            chunk_signature = signature([expected, chunk_unit])

            def make_chunk(directory, batch=batch):
                pq.write_table(
                    pa.Table.from_batches([batch]),
                    directory / "queries.parquet",
                    compression="zstd",
                )

            publish_unit(root, chunk_unit, chunk_signature, make_chunk)
            for shard in reference.manifest["shards"]:
                tasks.append(
                    {
                        "id": len(tasks),
                        "shard": shard["id"],
                        "query_chunk": chunk_unit,
                        "query_chunk_signature": chunk_signature,
                        "query_start": 0,
                        "query_count": min(
                            descriptor["query_chunk_size"],
                            query_count - chunk * descriptor["query_chunk_size"],
                        ),
                    }
                )
    else:
        tasks = [
            {
                "id": i,
                "shard": shard["id"],
                "query_start": 0,
                "query_count": query_count,
            }
            for i, shard in enumerate(reference.manifest["shards"])
        ]
    plan = {
        **descriptor,
        "descriptor": descriptor,
        "signature": expected,
        "tasks": tasks,
        "tasks_signature": signature(tasks),
        "query_count": query_count,
        "workers": effective_workers(workers, memory_gb),
        "memory_gb": memory_gb,
        "scratch": str(scratch or root / "scratch"),
    }
    plan["required_units"] = _required_units(plan)
    for task in tasks:
        atomic_json(
            root / "task-specs" / f"{task['id']:08d}.json",
            seal({"plan_signature": expected, "task": task}),
        )
    atomic_json(path, seal(plan))
    atomic_json(root / "worker-plan.json", _worker_header(plan, path))
    return path


def _worker_header(plan, path):
    keys = [
        *plan["descriptor"],
        "descriptor",
        "signature",
        "query_count",
        "workers",
        "memory_gb",
        "scratch",
    ]
    stamp = Path(path).stat()
    return seal(
        {
            **{key: plan[key] for key in keys},
            "task_count": len(plan["tasks"]),
            "plan_stamp": [stamp.st_size, stamp.st_mtime_ns],
        }
    )


def _validate_operation(plan):
    check_document(plan)
    descriptor = plan["descriptor"]
    if (
        descriptor.get("operation_version") != OPERATION_VERSION
        or descriptor.get("preparation_version") != PREPARATION_VERSION
    ):
        raise ValueError("unsupported operation/preparation version; create a new run")
    if signature(descriptor) != plan["signature"] or any(
        plan.get(key) != value for key, value in descriptor.items()
    ):
        raise ValueError(
            "plan signature does not match immutable operation configuration"
        )
    reference = MoleculeIndex.open(plan["index"])
    queries = MoleculeIndex.open(plan["queries"])
    if (
        reference.manifest["signature"] != plan["index_signature"]
        or queries.manifest["signature"] != plan["query_signature"]
    ):
        raise ValueError("plan index signatures changed")
    for key in ("chemistry", "fingerprint", "rdkit_version", "preparation_version"):
        if reference.manifest[key] != queries.manifest[key]:
            raise ValueError(f"incompatible plan indexes: {key}")
    count = (
        queries.manifest["alias_count"] if plan["kind"] == "search" else queries.size
    )
    if count != plan["query_count"]:
        raise ValueError("prepared query count changed")
    return reference


def _expected_task(plan, task_id, reference):
    shards = reference.manifest["shards"]
    count = (
        math.ceil(plan["query_count"] / plan["query_chunk_size"]) * len(shards)
        if plan["kind"] == "search"
        else len(shards)
    )
    if task_id < 0 or task_id >= count:
        raise ValueError("task_id outside plan")
    if plan["kind"] == "filter":
        return {
            "id": task_id,
            "shard": shards[task_id]["id"],
            "query_start": 0,
            "query_count": plan["query_count"],
        }
    chunk, shard = divmod(task_id, len(shards))
    unit = f"query-chunk-{chunk:08d}"
    return {
        "id": task_id,
        "shard": shards[shard]["id"],
        "query_chunk": unit,
        "query_chunk_signature": signature([plan["signature"], unit]),
        "query_start": 0,
        "query_count": min(
            plan["query_chunk_size"],
            plan["query_count"] - chunk * plan["query_chunk_size"],
        ),
    }


def _required_units(plan):
    units = [{"id": "queries", "signature": plan["signature"]}]
    if plan["kind"] == "search":
        for chunk in range(math.ceil(plan["query_count"] / plan["query_chunk_size"])):
            unit = f"query-chunk-{chunk:08d}"
            units.append(
                {"id": unit, "signature": signature([plan["signature"], unit])}
            )
    units.extend(
        {"id": f"task-{task['id']:08d}", "signature": _task_signature(plan, task)}
        for task in plan["tasks"]
    )
    units.append({"id": "merge", "signature": plan["signature"]})
    return units


def _load_plan(plan_path):
    path = Path(plan_path).resolve()
    plan = read_json(path)
    reference = _validate_operation(plan)
    if signature(plan["tasks"]) != plan["tasks_signature"]:
        raise ValueError("task manifest changed")
    if plan["tasks"] != [
        _expected_task(plan, i, reference) for i in range(len(plan["tasks"]))
    ]:
        raise ValueError("duplicated or invalid task IDs")
    expected_count = (
        math.ceil(plan["query_count"] / plan["query_chunk_size"])
        * len(reference.manifest["shards"])
        if plan["kind"] == "search"
        else len(reference.manifest["shards"])
    )
    if len(plan["tasks"]) != expected_count or plan[
        "required_units"
    ] != _required_units(plan):
        raise ValueError("incomplete task/dependency manifest")
    if read_json(path.parent / "worker-plan.json") != _worker_header(plan, path):
        raise ValueError("worker manifest changed")
    for task in plan["tasks"]:
        expected = seal({"plan_signature": plan["signature"], "task": task})
        if read_json(path.parent / "task-specs" / f"{task['id']:08d}.json") != expected:
            raise ValueError("task specification changed or incomplete")
    return path.parent, plan, reference


def _load_worker(plan_path, task_id):
    root = Path(plan_path).resolve().parent
    plan = read_json(root / "worker-plan.json")
    reference = _validate_operation(plan)
    stamp = Path(plan_path).stat()
    if plan["plan_stamp"] != [stamp.st_size, stamp.st_mtime_ns]:
        raise ValueError("coordinator plan changed; replan before executing workers")
    task = _expected_task(plan, task_id, reference)
    if task_id >= plan["task_count"]:
        raise ValueError("task_id outside plan")
    spec = read_json(root / "task-specs" / f"{task_id:08d}.json")
    if spec != seal({"plan_signature": plan["signature"], "task": task}):
        raise ValueError("task specification changed")
    return root, plan, reference, task


def _task_signature(plan, task):
    return signature([plan["signature"], task])


def _queries(root, task):
    iterator = batches(
        root / "units" / task.get("query_chunk", "queries") / "queries.parquet"
    )
    yield from itertools.islice(
        iterator, task["query_start"], task["query_start"] + task["query_count"]
    )


def _best(rows, k):
    # Exact rational primary score and stable molecule ID tie break.
    return heapq.nsmallest(
        k,
        rows,
        key=lambda row: (
            -Fraction(row["intersection"], row["union"]),
            row["molecule_id"],
        ),
    )


def run_task(
    plan_path: str | Path, task_id: int, *, workers: int | None = None
) -> Path:
    """Run one planned task and return its checksummed completion-receipt path.

    Accepts build, search, or filter plans. A worker override is capped by the
    scheduler allocation and plan memory budget. Missing or damaged work is
    repaired; a valid completed unit is reused. Invalid IDs, incompatible plans,
    and corrupt inputs raise before publication. Run only one copy of each task
    at a time; sequential retries are safe.

    Example: ``mst.run_task('runs/candidates/plan.json', 0, workers=4)``.
    """
    return _run_task(plan_path, task_id, workers=workers)


def _run_task(
    plan_path: str | Path,
    task_id: int,
    *,
    workers: int | None = None,
    _resident=None,
    _context=None,
) -> Path:
    """Execute one idempotent scheduler-independent task and publish its receipt.

    No shared HDF5 or database is written by workers. All inner FPSim2 kernels
    are serial; a read-only shard engine is shared across outer query threads.
    """
    worker_file = Path(plan_path).resolve().parent / "worker-plan.json"
    if _context is None and not worker_file.exists():
        if read_json(Path(plan_path))["kind"] == "build":
            return run_build_task(plan_path, task_id)
        raise ValueError("incomplete operation worker manifest")
    if _context is None:
        root, plan, reference, task = _load_worker(plan_path, task_id)
    else:
        root, plan, reference = _context
        if task_id < 0 or task_id >= len(plan["tasks"]):
            raise ValueError("task_id outside plan")
        task = plan["tasks"][task_id]
    unit = f"task-{task_id:08d}"
    expected = _task_signature(plan, task)
    shard_info = next(
        shard for shard in reference.manifest["shards"] if shard["id"] == task["shard"]
    )
    if plan["kind"] == "filter" and not receipt_valid(
        root, "queries", plan["signature"]
    ):
        raise ValueError("staged queries are missing or corrupt")
    if task.get("query_chunk") and not receipt_valid(
        root, task["query_chunk"], task["query_chunk_signature"]
    ):
        raise ValueError("query chunk is missing or corrupt")
    cache = {} if _resident is None else _resident
    if task["shard"] not in cache and not receipt_valid(
        reference.path, shard_info["id"], shard_info["signature"]
    ):
        raise ValueError("pool shard is missing or corrupt; resume its build")

    def make(directory):
        shard_dir = reference.path / "units" / task["shard"]
        # Reserve resident buffers, witness/heap objects and worst-case candidate
        # arrays before allocating the shard. Memory is an execution budget, not
        # an OS-enforced hard limit; allocator/runtime overhead retains headroom.
        n = shard_info["count"]
        fp_bytes = reference.fingerprint.n_bits // 8
        resident_bytes = n * (2 * fp_bytes + 160) + 128_000_000
        state_bytes = n * 800 if plan["kind"] == "filter" else 0
        query_bytes = max(
            16_000_000,
            n * 40 + (min(n, plan["k"]) * 900 if plan["kind"] == "search" else 0),
        )
        available = plan["memory_gb"] * 2**30 * 0.8 - resident_bytes - state_bytes
        if available < query_bytes:
            raise MemoryError(
                "shard/query state exceeds memory budget; rebuild with smaller shard_size or increase memory_gb"
            )
        worker_count = min(
            effective_workers(workers or plan["workers"], plan["memory_gb"]),
            max(1, int(available // query_bytes)),
        )
        if task["shard"] not in cache:
            table = pq.read_table(
                shard_dir / "molecules.parquet",
                columns=["molecule_id", "fingerprint"],
                use_threads=False,
            )
            engine = ShardEngine(shard_dir / "fingerprints.h5")
            cache[task["shard"]] = table, engine
        table, engine = cache[task["shard"]]

        def molecule(local_id):
            offset = int(local_id) - 1
            if not 0 <= offset < len(table):
                raise ValueError("FPSim2 returned a nonlocal molecule ID")
            return table["molecule_id"][offset].as_py(), table["fingerprint"][
                offset
            ].as_py()

        if plan["kind"] == "search":

            def query_one(query):
                candidates = engine.top_k_candidates(
                    query["fingerprint"],
                    min(shard_info["count"], plan["k"] + int(plan["exclude_identity"])),
                )

                def scored():
                    for local_id in candidates:
                        mol_id, fp = molecule(local_id)
                        if plan["exclude_identity"] and mol_id == query["molecule_id"]:
                            continue
                        a, b = intersection_union(query["fingerprint"], fp)
                        yield {
                            "query_id": query["query_id"],
                            "query_original_id": query["query_original_id"],
                            "query_molecule_id": query["molecule_id"],
                            "molecule_id": mol_id,
                            "rank": 0,
                            "similarity": a / b,
                            "intersection": a,
                            "union": b,
                        }

                return _best(scored(), plan["k"])

            with pq.ParquetWriter(
                directory / "matches.parquet", MATCH_SCHEMA, compression="zstd"
            ) as writer:
                query_iterator = _queries(root, task)
                pending = []
                with ThreadPoolExecutor(max_workers=worker_count) as pool:
                    # Bounded submissions: Executor.map eagerly consumes inputs before Python 3.14.
                    while block := list(itertools.islice(query_iterator, worker_count)):
                        for rows in pool.map(query_one, block):
                            pending.extend(rows)
                            if len(pending) >= 8192:
                                writer.write_table(
                                    pa.Table.from_pylist(pending, schema=MATCH_SCHEMA)
                                )
                                pending = []
                if pending:
                    writer.write_table(
                        pa.Table.from_pylist(pending, schema=MATCH_SCHEMA)
                    )
        else:
            threshold = Fraction(plan["threshold"])
            witnesses = {}

            def threshold_query(query):
                return query, engine.threshold_candidates(
                    query["fingerprint"], float(threshold)
                )

            query_iterator = _queries(root, task)
            with ThreadPoolExecutor(max_workers=worker_count) as pool:
                while len(witnesses) < n and (
                    block := list(itertools.islice(query_iterator, worker_count))
                ):
                    # map preserves protected ID order, hence the first witness
                    # is the lowest ID even when backend queries finish out of order.
                    for query, candidates in pool.map(threshold_query, block):
                        for local_id in candidates:
                            local_id = int(local_id)
                            if local_id in witnesses:
                                continue
                            mol_id, fp = molecule(local_id)
                            a, b = intersection_union(query["fingerprint"], fp)
                            hit = (
                                a * threshold.denominator > threshold.numerator * b
                                if plan["comparison"] == ">"
                                else a * threshold.denominator
                                >= threshold.numerator * b
                            )
                            if hit:
                                witnesses[local_id] = {
                                    "molecule_id": mol_id,
                                    "protected_id": query["molecule_id"],
                                    "similarity": a / b,
                                    "intersection": a,
                                    "union": b,
                                }
                        if len(witnesses) == n:
                            break
            with pq.ParquetWriter(
                directory / "rejected.parquet", HIT_SCHEMA, compression="zstd"
            ) as writer:
                ordered = (witnesses[key] for key in sorted(witnesses))
                while block := list(itertools.islice(ordered, 8192)):
                    writer.write_table(pa.Table.from_pylist(block, schema=HIT_SCHEMA))

    publish_unit(root, unit, expected, make)
    return root / "receipts" / f"{unit}.json"


def merge(plan_path: str | Path) -> Result | MoleculeIndex:
    """Validate every planned receipt then atomically reduce complete task outputs."""
    if read_json(Path(plan_path))["kind"] == "build":
        return merge_build(plan_path)
    root, plan, reference = _load_plan(plan_path)
    from .inspection import run_dependencies

    run_dependencies(root, plan, reference, include_merge=False)
    expected_units = {f"task-{task['id']:08d}" for task in plan["tasks"]}
    existing_units = {path.stem for path in (root / "receipts").glob("task-*.json")}
    if existing_units - expected_units:
        raise ValueError("unexpected or duplicated task receipts")
    for task in plan["tasks"]:
        if not receipt_valid(
            root, f"task-{task['id']:08d}", _task_signature(plan, task)
        ):
            raise ValueError(f"incomplete or corrupt task {task['id']}")

    def reduce(directory):
        con = connect(plan["memory_gb"], plan["scratch"])
        try:
            task_files = [
                root
                / "units"
                / f"task-{task['id']:08d}"
                / (
                    "matches.parquet"
                    if plan["kind"] == "search"
                    else "rejected.parquet"
                )
                for task in plan["tasks"]
            ]
            output_count = 0
            if plan["kind"] == "search":
                if task_files:
                    relation = (
                        "read_parquet(["
                        + ",".join(sql_path(path) for path in task_files)
                        + "])"
                    )
                    cursor = con.execute(
                        f"SELECT * FROM {relation} ORDER BY query_id,molecule_id"
                    )
                    with pq.ParquetWriter(
                        directory / "matches.parquet", MATCH_SCHEMA, compression="zstd"
                    ) as writer:

                        def ordered_rows():
                            while block := cursor.fetchmany(8192):
                                for values in block:
                                    yield dict(zip(MATCH_SCHEMA.names, values))

                        pending = []
                        for query_id, rows in itertools.groupby(
                            ordered_rows(), key=lambda row: row["query_id"]
                        ):
                            selected = _best(rows, plan["k"])
                            for rank, row in enumerate(selected, 1):
                                row["rank"] = rank
                            output_count += len(selected)
                            pending.extend(selected)
                            if len(pending) >= 8192:
                                writer.write_table(
                                    pa.Table.from_pylist(pending, schema=MATCH_SCHEMA)
                                )
                                pending = []
                        if pending:
                            writer.write_table(
                                pa.Table.from_pylist(pending, schema=MATCH_SCHEMA)
                            )
                else:
                    write_rows(directory / "matches.parquet", [], MATCH_SCHEMA)
                # Invalid query rows are retained rather than silently dropped.
                query_index = MoleculeIndex.open(plan["queries"])
                con.execute(
                    f"COPY (SELECT * FROM read_parquet({sql_path(query_index.invalid_path)})) TO {sql_path(directory / 'invalid.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
                )
                summary = {
                    "kind": "search",
                    "query_count": plan["query_count"],
                    "pool_count": reference.size,
                    "match_count": output_count,
                    "invalid_query_count": query_index.manifest["invalid_count"],
                    "k": plan["k"],
                    "exclude_identity": plan["exclude_identity"],
                }
            else:
                if task_files:
                    relation = (
                        "read_parquet(["
                        + ",".join(sql_path(path) for path in task_files)
                        + "])"
                    )
                    con.execute(
                        f"COPY (SELECT * FROM {relation} ORDER BY molecule_id) TO {sql_path(directory / 'rejected.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
                    )
                else:
                    write_rows(directory / "rejected.parquet", [], HIT_SCHEMA)
                con.execute(
                    f"COPY (SELECT m.* EXCLUDE(fingerprint,ordinal) FROM read_parquet({sql_path(reference.molecules_path)}) m ANTI JOIN read_parquet({sql_path(directory / 'rejected.parquet')}) r USING(molecule_id) ORDER BY molecule_id) TO {sql_path(directory / 'eligible.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
                )
                con.execute(
                    f"COPY (SELECT * FROM read_parquet({sql_path(reference.invalid_path)})) TO {sql_path(directory / 'invalid.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
                )
                con.execute(
                    f"COPY (SELECT a.*, CASE WHEN r.molecule_id IS NULL THEN 'eligible' ELSE 'rejected' END AS disposition FROM read_parquet({sql_path(reference.aliases_path)}) a LEFT JOIN read_parquet({sql_path(directory / 'rejected.parquet')}) r USING(molecule_id) ORDER BY source,source_row) TO {sql_path(directory / 'aliases.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)"
                )
                rejected = pq.ParquetFile(
                    directory / "rejected.parquet"
                ).metadata.num_rows
                summary = {
                    "kind": "filter",
                    "pool_count": reference.size,
                    "eligible_count": reference.size - rejected,
                    "rejected_count": rejected,
                    "invalid_count": reference.manifest["invalid_count"],
                    "protected_count": plan["query_count"],
                    "protected_invalid_count": MoleculeIndex.open(
                        plan["queries"]
                    ).manifest["invalid_count"],
                    "threshold": plan["threshold"],
                    "comparison": plan["comparison"],
                }
            atomic_json(directory / "summary.json", summary)
            failpoint("merge:after_write")
        finally:
            con.close()

    publish_unit(root, "merge", plan["signature"], reduce)
    return Result(root, plan["kind"])


def execute(
    plan_path: str | Path, *, workers: int | None = None
) -> Result | MoleculeIndex:
    """Execute an immutable build/search/filter plan locally and publish its merge.

    Uses the same receipts and atomic task outputs as external schedulers. Search
    and filter execution retains one read-only resident shard at a time, shared
    by bounded query threads. Missing/corrupt task outputs are rerun; successful
    units are reused. ``memory_gb`` in the plan is measured in GiB.
    """
    if read_json(Path(plan_path))["kind"] == "build":
        context = _build_plan(plan_path)
        for task in context[1]["tasks"]:
            run_build_task(plan_path, task["id"], _context=context)
        return merge_build(plan_path)
    context = _load_plan(plan_path)
    _root, plan, _ = context
    check_units(_root, [u for u in plan["required_units"] if u["id"] == "queries"])
    resident = {}
    shard_id = None
    for task in sorted(plan["tasks"], key=lambda task: (task["shard"], task["id"])):
        if task["shard"] != shard_id:
            resident.clear()
            shard_id = task["shard"]
        _run_task(
            plan_path, task["id"], workers=workers, _resident=resident, _context=context
        )
    resident.clear()
    return merge(plan_path)


def search(
    index: MoleculeIndex | str | Path,
    queries: Source | list[Source] | MoleculeIndex | str | Path,
    k: int = 10,
    out: str | Path = "search-results",
    *,
    exclude_identity: bool = False,
    workers: int = 4,
    memory_gb: float = 4,
    scratch: str | Path | None = None,
    query_chunk_size: int = 256,
) -> Result:
    """Exact top-k binary Morgan Tanimoto search across every shard.

    Scores are ranked using integer intersection/union ratios, with molecule ID
    tie breaks; canonical identity exclusion happens before final truncation.
    Results are streamed Parquet with query IDs, original IDs, molecule IDs,
    rank, similarity, intersection and union. No query × corpus matrix is built.
    """
    return execute(
        plan_search(
            index,
            queries,
            k,
            out,
            exclude_identity=exclude_identity,
            workers=workers,
            memory_gb=memory_gb,
            scratch=scratch,
            query_chunk_size=query_chunk_size,
        ),
        workers=workers,
    )


def filter_pool(
    index: MoleculeIndex | str | Path,
    protected: Source | list[Source] | MoleculeIndex | str | Path,
    threshold: float | str = 0.8,
    comparison: str = ">",
    out: str | Path = "filter-results",
    *,
    invalid_protected: str = "error",
    workers: int = 4,
    memory_gb: float = 4,
    scratch: str | Path | None = None,
) -> Result:
    """Reject each pool molecule matching any protected structure at the boundary.

    Default rejection means strictly greater than 0.8, evaluated exactly with
    integer counts. Protected invalids fail by default; empty protected sets
    always fail. Output tables are eligible, rejected, invalid and aliases, plus
    summary JSON. The rejected table records a deterministic protected witness.
    """
    return execute(
        plan_filter(
            index,
            protected,
            out,
            threshold=threshold,
            comparison=comparison,
            invalid_protected=invalid_protected,
            workers=workers,
            memory_gb=memory_gb,
            scratch=scratch,
        ),
        workers=workers,
    )
