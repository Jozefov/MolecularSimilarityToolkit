"""Lazy Parquet result handles and bounded CSV/TSV export."""

import csv
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .storage import read_json


@dataclass(frozen=True)
class Result:
    """Disk-backed result; ``batches`` returns an Arrow RecordBatch iterator."""

    path: Path
    kind: str

    @property
    def summary(self) -> dict:
        return read_json(self.path / "units/merge/summary.json")

    def table_path(self, table: str | None = None) -> Path:
        return (
            self.path
            / "units/merge"
            / f"{table or ('matches' if self.kind == 'search' else 'eligible')}.parquet"
        )

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
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = None
            for batch in self.batches(table):
                if writer is None:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=batch.schema.names,
                        delimiter=delimiter or ("\t" if path.suffix == ".tsv" else ","),
                    )
                    writer.writeheader()
                writer.writerows(batch.to_pylist())
            if writer is None:
                schema = pq.ParquetFile(self.table_path(table)).schema_arrow
                csv.writer(
                    handle,
                    delimiter=delimiter or ("\t" if path.suffix == ".tsv" else ","),
                ).writerow(schema.names)
        return path
