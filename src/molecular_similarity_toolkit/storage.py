"""Atomic publications and checksummed receipts shared by local and batch work."""

import fcntl
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq


def digest(path: str | Path) -> str:
    checksum = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            checksum.update(block)
    return checksum.hexdigest()


def signature(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()
    ).hexdigest()


def atomic_json(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True, default=str) + "\n"
    )
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
        files = value["files"]
        if value["signature"] != expected or not isinstance(files, list) or not files:
            return False
        target = root / "units" / unit
        actual = {
            p.relative_to(root).as_posix() for p in target.rglob("*") if p.is_file()
        }
        recorded = [item["path"] for item in files]
        if len(set(recorded)) != len(recorded) or set(recorded) != actual:
            return False
        for item in files:
            path = root / item["path"]
            checksum = item["sha256"]
            if not path.resolve().is_relative_to(target.resolve()):
                return False
            if (
                not isinstance(checksum, str)
                or len(checksum) != 64
                or any(c not in "0123456789abcdef" for c in checksum)
            ):
                return False
            if digest(path) != checksum:
                return False
        return True
    except (ValueError, KeyError, TypeError, OSError):
        return False


def _publish_locked(root: Path, unit: str, expected: str, make):
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
        files = [
            {"path": str(p.relative_to(root)), "sha256": digest(p)}
            for p in sorted(target.rglob("*"))
            if p.is_file()
        ]
        atomic_json(
            root / "receipts" / f"{unit}.json", {"signature": expected, "files": files}
        )
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def publish_unit(root: Path, unit: str, expected: str, make):
    """Publish one immutable unit under a crash-released POSIX advisory lock."""
    locks = root / ".locks"
    locks.mkdir(parents=True, exist_ok=True)
    with (locks / f"{unit}.lock").open("a") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(
                f"unit {unit} is already running; retry after it completes"
            ) from error
        try:
            return _publish_locked(root, unit, expected, make)
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def write_rows(path: Path, rows: list[dict], schema: pa.Schema):
    pq.write_table(pa.Table.from_pylist(rows, schema=schema), path, compression="zstd")


def batches(path: Path, batch_size=8192):
    for batch in pq.ParquetFile(path).iter_batches(
        batch_size=batch_size, use_threads=False
    ):
        yield from batch.to_pylist()


def connect(memory_gb=4, scratch=None):
    import duckdb

    connection = duckdb.connect()
    connection.execute(f"SET memory_limit='{float(memory_gb)}GiB'")
    connection.execute("SET threads=1")
    connection.execute("SET preserve_insertion_order=false")
    if scratch:
        Path(scratch).mkdir(parents=True, exist_ok=True)
        connection.execute("SET temp_directory=?", [str(Path(scratch).resolve())])
    return connection


def sql_path(path: Path) -> str:
    return "'" + str(path).replace("'", "''") + "'"


def seal(value: dict) -> dict:
    """Checksum an immutable JSON document (integrity, not authentication)."""
    value = {key: item for key, item in value.items() if key != "integrity"}
    return {**value, "integrity": signature(value)}


def check_document(value: dict):
    if value.get("integrity") != seal(value)["integrity"]:
        raise ValueError("manifest changed or incomplete document; replan/rebuild")


def check_units(root: Path, required: list[dict], *, skip=(), exact=False):
    """One dependency contract for execution, reduction and full verification."""
    identifiers = [unit["id"] for unit in required]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("duplicated required units")
    if exact and {p.stem for p in (root / "receipts").glob("*.json")} - set(
        identifiers
    ):
        raise ValueError("unexpected or duplicated unit receipts")
    for unit in required:
        if unit["id"] not in skip and not receipt_valid(
            root, unit["id"], unit["signature"]
        ):
            raise ValueError(f"required unit {unit['id']} is incomplete or corrupt")
