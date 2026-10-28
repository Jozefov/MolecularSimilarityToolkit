"""Streaming source readers; no corpus-sized dataframe is constructed."""
import csv
import gzip
import json
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
        for batch in pq.ParquetFile(path).iter_batches(batch_size=8192):
            yield from batch.to_pylist()
        return
    opener = gzip.open if suffix.endswith(".gz") else open
    plain = suffix[:-3] if suffix.endswith(".gz") else suffix
    if plain.endswith((".sdf", ".sd")):
        with opener(path, "rb") as handle:
            supplier = Chem.ForwardSDMolSupplier(handle, removeHs=False)
            for mol in supplier:
                if mol is None:
                    yield {source.smiles_column: "", "_reader_error": "invalid SDF record"}
                else:
                    row = mol.GetPropsAsDict()
                    row[source.smiles_column] = Chem.MolToSmiles(mol)
                    row.setdefault(source.id_column, mol.GetProp("_Name") if mol.HasProp("_Name") else "")
                    yield row
        return
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        if plain.endswith((".smi", ".smiles")):
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split(maxsplit=1)
                yield {source.smiles_column: fields[0], source.id_column: fields[1] if len(fields) > 1 else ""}
        elif plain.endswith((".csv", ".tsv", ".txt")):
            yield from csv.DictReader(handle, delimiter="\t" if plain.endswith((".tsv", ".txt")) else ",")
        else:
            raise ValueError(f"unsupported input format: {path}")


def selected_rows(source: Source):
    for number, row in enumerate(records(source), 1):
        if source.splits is not None and str(row.get(source.split_column)) not in source.splits:
            continue
        yield number, row
