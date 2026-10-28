"""Streaming source readers; no corpus-sized dataframe is constructed."""

import csv
import gzip
from pathlib import Path

import pyarrow.parquet as pq
from rdkit import Chem

from .config import Source


def records(source: Source):
    """Yield original rows, retaining malformed SDF records for quarantine."""
    if source._records is not None:
        yield from source._records
        return
    path = Path(source.path)
    suffix = path.name.lower()
    if suffix.endswith(".parquet"):
        file = pq.ParquetFile(path)
        _columns(source, file.schema_arrow.names)
        for batch in file.iter_batches(batch_size=8192, use_threads=False):
            yield from batch.to_pylist()
        return
    opener = gzip.open if suffix.endswith(".gz") else open
    plain = suffix.removesuffix(".gz")
    if plain.endswith((".sdf", ".sd")):
        with opener(path, "rb") as handle:
            supplier = Chem.ForwardSDMolSupplier(handle, removeHs=False)
            for mol in supplier:
                if mol is None:
                    yield {
                        source.smiles_column: "",
                        "_reader_error": "invalid SDF record",
                    }
                else:
                    row = mol.GetPropsAsDict(autoConvertStrings=False)
                    row[source.smiles_column] = Chem.MolToSmiles(mol)
                    if source.id_column is None:
                        row.setdefault(
                            "id", mol.GetProp("_Name") if mol.HasProp("_Name") else ""
                        )
                    yield row
        return
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        if plain.endswith((".smi", ".smiles")):
            for line in handle:
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                fields = line.rstrip("\n").split(maxsplit=1)
                yield {
                    source.smiles_column: fields[0],
                    (source.id_column or "id"): fields[1] if len(fields) > 1 else "",
                }
        elif plain.endswith((".csv", ".tsv", ".txt")):
            reader = csv.DictReader(
                handle, delimiter="\t" if plain.endswith((".tsv", ".txt")) else ","
            )
            _columns(source, reader.fieldnames or [])
            yield from reader
        else:
            raise ValueError(f"unsupported input format: {path}")


def selected_rows(source: Source):
    validated = False
    for number, row in enumerate(records(source), 1):
        # Unknown fold membership is never a reason to discard a reader failure.
        if "_reader_error" in row:
            yield number, row
            continue
        required = [source.smiles_column, *source.metadata_columns]
        if source.split_column:
            required.append(source.split_column)
        if source.id_column is not None:
            required.append(source.id_column)
        missing = [name for name in required if name not in row]
        if (not validated and missing) or (
            source.id_column is not None and source.id_column not in row
        ):
            raise ValueError(
                f"source {source.name or source.path}: missing mapped columns {missing}"
            )
        validated = True
        if source.splits is not None:
            if (
                source.split_column not in row
                or row[source.split_column] is None
                or not str(row[source.split_column]).strip()
            ):
                yield number, {**row, "_reader_error": "missing split membership"}
                continue
            if str(row[source.split_column]) not in source.splits:
                continue
        yield number, row


def _columns(source, available):
    required = [source.smiles_column, *source.metadata_columns]
    if source.split_column:
        required.append(source.split_column)
    if source.id_column:
        required.append(source.id_column)
    missing = [column for column in required if column not in available]
    if missing:
        raise ValueError(
            f"source {source.name or source.path}: missing mapped columns {missing}"
        )
