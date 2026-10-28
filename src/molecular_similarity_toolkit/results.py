"""Lazy Parquet result handles and bounded CSV/TSV export."""

import csv
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .storage import check_document, read_json, receipt_valid


@dataclass(frozen=True)
class Result:
    """Disk-backed result; ``batches`` returns an Arrow RecordBatch iterator."""

    path: Path
    kind: str

    @classmethod
    def open(cls, path: str | Path) -> "Result":
        """Open a completed result, validating its published table receipt.

        Use ``verify`` to include every task and external input dependency.
        """
        path = Path(path).resolve()
        plan = read_json(path / "plan.json")
        check_document(plan)
        if plan["kind"] not in {"search", "filter"} or not receipt_valid(
            path, "merge", plan["signature"]
        ):
            raise ValueError("result is incomplete or corrupt")
        return cls(path, plan["kind"])

    @property
    def summary(self) -> dict:
        return read_json(self.path / "units/merge/summary.json")

    def table_path(self, table: str | None = None) -> Path:
        table = table or ("matches" if self.kind == "search" else "eligible")
        allowed = (
            {"matches", "invalid"}
            if self.kind == "search"
            else {"eligible", "rejected", "invalid", "aliases"}
        )
        if table not in allowed:
            raise ValueError(
                f"unknown {self.kind} table {table!r}; choose from {sorted(allowed)}"
            )
        return self.path / "units/merge" / f"{table}.parquet"

    def batches(
        self, table: str | None = None, *, batch_size: int = 8192
    ) -> Iterator[pa.RecordBatch]:
        """Stream search matches or filter eligible/rejected/invalid rows."""
        return pq.ParquetFile(self.table_path(table)).iter_batches(
            batch_size=batch_size
        )

    def export(
        self,
        path: str | Path,
        *,
        table: str | None = None,
        delimiter: str | None = None,
    ) -> Path:
        """Export one result table to CSV or TSV with bounded memory."""
        path = Path(path)
        if path.suffix.lower() not in {".csv", ".tsv"}:
            raise ValueError("export path must end in .csv or .tsv")
        if delimiter is not None and (len(delimiter) != 1 or delimiter in "\r\n"):
            raise ValueError("delimiter must be one non-newline character")
        self.table_path(table)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = None
            for batch in self.batches(table):
                if writer is None:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=batch.schema.names,
                        delimiter=delimiter
                        or ("\t" if path.suffix.lower() == ".tsv" else ","),
                    )
                    writer.writeheader()
                writer.writerows(batch.to_pylist())
            if writer is None:
                schema = pq.ParquetFile(self.table_path(table)).schema_arrow
                csv.writer(
                    handle,
                    delimiter=delimiter
                    or ("\t" if path.suffix.lower() == ".tsv" else ","),
                ).writerow(schema.names)
        return path
