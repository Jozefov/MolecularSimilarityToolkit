"""Atomic publications and checksummed receipts shared by local and batch work."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import pyarrow as pa
import pyarrow.parquet as pq


def digest(path: str | Path) -> str:
    checksum = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def signature(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def atomic_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")
    os.replace(temporary, path)


def read_json(path: Path):
    return json.loads(path.read_text())


def failpoint(name: str):
    """Test hook: MST_FAILPOINT=name injects an interruption at an atomic boundary."""
    if os.environ.get("MST_FAILPOINT") == name:
        raise RuntimeError(f"injected interruption: {name}")


def receipt_valid(root: Path, unit: str, expected: str) -> bool:
    receipt = root / "receipts" / f"{unit}.json"
    if not receipt.exists():
        return False
    try:
        value = read_json(receipt)
        return value["signature"] == expected and all((root / item["path"]).is_file() and digest(root / item["path"]) == item["sha256"] for item in value["files"])
    except (ValueError, KeyError, OSError):
        return False


def publish_unit(root: Path, unit: str, expected: str, make):
    """Publish a unit directory then a receipt. Missing receipts safely rerun a unit."""
    if receipt_valid(root, unit, expected):
        return
    target = root / "units" / unit
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{unit}-", dir=target.parent))
    try:
        make(temporary)
        failpoint(f"{unit}:before_publish")
        if target.exists():
            shutil.rmtree(target)
        os.replace(temporary, target)
        failpoint(f"{unit}:after_publish")
        files = [{"path": str(p.relative_to(root)), "sha256": digest(p)} for p in sorted(target.rglob("*")) if p.is_file()]
        atomic_json(root / "receipts" / f"{unit}.json", {"signature": expected, "files": files})
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def write_rows(path: Path, rows: list[dict], schema: pa.Schema):
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path, compression="zstd")


def batches(path: Path, batch_size=8192):
    for batch in pq.ParquetFile(path).iter_batches(batch_size=batch_size):
        yield from batch.to_pylist()


def connect(memory_gb=4, scratch=None):
    import duckdb
    connection = duckdb.connect()
    connection.execute(f"SET memory_limit='{float(memory_gb)}GB'")
    connection.execute("SET threads=1")
    connection.execute("SET preserve_insertion_order=false")
    if scratch:
        Path(scratch).mkdir(parents=True, exist_ok=True)
        connection.execute("SET temp_directory=?", [str(Path(scratch).resolve())])
    return connection


def sql_path(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"
