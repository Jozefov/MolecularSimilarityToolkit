"""Resumable streaming normalization, disk-backed deduplication and immutable shards."""
from dataclasses import dataclass
import json
import math
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from .adapter import write_shard
from .chemistry import prepare, RDKIT_VERSION
from .config import ChemistryConfig, FingerprintConfig, Source, validate_configs
from .readers import selected_rows
from .storage import atomic_json, batches, connect, digest, publish_unit, read_json, signature, sql_path, write_rows

STAGE_SCHEMA = pa.schema([("row_id", pa.string()), ("source", pa.string()), ("source_row", pa.int64()), ("original_id", pa.string()), ("original_smiles", pa.string()), ("split", pa.string()), ("metadata", pa.string()), ("canonical_smiles", pa.string()), ("molecule_id", pa.string()), ("fingerprint", pa.binary()), ("popcount", pa.int64()), ("error", pa.string())])
MOLECULE_SCHEMA = pa.schema([("molecule_id", pa.string()), ("canonical_smiles", pa.string()), ("fingerprint", pa.binary()), ("popcount", pa.int64()), ("local_id", pa.int64())])


def source_descriptors(sources):
    values = []
    for index, source in enumerate(sources):
        value = source.to_dict()
        value["name"] = source.name or (Path(source.path).name if source.path else f"source-{index}")
        value["checksum"] = digest(source.path) if source.path else signature(value["_records"])
        values.append(value)
    names = [v["name"] for v in values]
    if len(set(names)) != len(names):
        raise ValueError("source names must be unique")
    return values


def stage_sources(sources, root: Path, chemistry, fingerprint, *, chunk_size=8192):
    paths = []
    for source_index, source in enumerate(sources):
        name = source.name or (Path(source.path).name if source.path else f"source-{source_index}")
        pending = []
        chunk = 0
        def flush(items, number):
            unit = f"stage-{source_index:04d}-{number:08d}"
            expected = signature({"root": read_json(root / "build.json")["signature"], "unit": unit})
            def make(directory):
                normalized = []
                for source_row, original in items:
                    row = {"row_id": signature([source_index, source_row]), "source": name, "source_row": source_row,
                           "original_id": str(original.get(source.id_column, "") or ""), "original_smiles": str(original.get(source.smiles_column, "") or ""),
                           "split": None if source.split_column is None else str(original.get(source.split_column, "")),
                           "metadata": json.dumps({key: original.get(key) for key in source.metadata_columns}, sort_keys=True, default=str)}
                    try:
                        canonical, mol_id, fp, count = prepare(row["original_smiles"], chemistry, fingerprint)
                        row.update(canonical_smiles=canonical, molecule_id=mol_id, fingerprint=fp, popcount=count, error=None)
                    except Exception as error:
                        row.update(canonical_smiles=None, molecule_id=None, fingerprint=None, popcount=None, error=str(error))
                    normalized.append(row)
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
        empty = root / "empty.parquet"
        write_rows(empty, [], STAGE_SCHEMA)
        paths.append(empty)
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
        if manifest["rdkit_version"] != RDKIT_VERSION:
            raise ValueError("index RDKit version differs from runtime; rebuild the index")
        return cls(path, manifest)

    @property
    def chemistry(self):
        return ChemistryConfig.from_dict(self.manifest["chemistry"])

    @property
    def fingerprint(self):
        return FingerprintConfig.from_dict(self.manifest["fingerprint"])

    @property
    def size(self):
        return self.manifest["molecule_count"]

    @property
    def molecules_path(self):
        return self.path / "units" / "catalog" / "molecules.parquet"

    @property
    def aliases_path(self):
        return self.path / "units" / "catalog" / "aliases.parquet"

    @property
    def invalid_path(self):
        return self.path / "units" / "catalog" / "invalid.parquet"

    def molecules(self, batch_size=8192):
        """Iterate Arrow batches without materializing the catalog."""
        return pq.ParquetFile(self.molecules_path).iter_batches(batch_size=batch_size)


def build_index(sources, out: str | Path, *, chemistry: ChemistryConfig = ChemistryConfig(), fingerprint: FingerprintConfig = FingerprintConfig(), shard_size: int = 1_000_000, workers: int = 4, memory_gb: float = 4, scratch: str | Path | None = None, resume: bool = True, chunk_size: int = 8192) -> MoleculeIndex:
    """Build/resume an immutable index; changed source bytes or policies are rejected.

    ``sources`` is a Source or iterable of Sources. Stage fingerprints are reused
    on resume. Deduplication sorts on disk through DuckDB with a memory limit.
    Each shard has at most one million local integer IDs, irrespective of corpus
    size. Invalid source rows remain in ``invalid.parquet`` with provenance.
    """
    validate_configs(chemistry, fingerprint)
    if not 1 <= shard_size <= 1_000_000 or memory_gb <= 0 or workers < 1 or chunk_size < 1:
        raise ValueError("invalid shard_size, workers, memory_gb or chunk_size")
    sources = [sources] if isinstance(sources, Source) else list(sources)
    root = Path(out).resolve()
    root.mkdir(parents=True, exist_ok=True)
    descriptor = {"format_version": 1, "sources": source_descriptors(sources), "chemistry": chemistry.to_dict(), "fingerprint": fingerprint.to_dict(), "rdkit_version": RDKIT_VERSION, "shard_size": shard_size, "chunk_size": chunk_size}
    expected = signature(descriptor)
    if (root / "build.json").exists():
        if not resume or read_json(root / "build.json")["signature"] != expected:
            raise ValueError("resume refused: inputs or configuration changed, or resume disabled")
    else:
        atomic_json(root / "build.json", {**descriptor, "signature": expected})
    paths = stage_sources(sources, root, chemistry, fingerprint, chunk_size=chunk_size)
    def catalog(directory):
        con = connect(memory_gb, scratch or root / "scratch")
        try:
            con.execute("CREATE VIEW staged AS SELECT * FROM read_parquet([" + ",".join(sql_path(p) for p in paths) + "])")
            con.execute(f"COPY (SELECT DISTINCT molecule_id,canonical_smiles,fingerprint,popcount FROM staged WHERE error IS NULL ORDER BY molecule_id) TO {sql_path(directory / 'molecules.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)")
            con.execute(f"COPY (SELECT row_id,source,source_row,original_id,original_smiles,split,metadata,molecule_id FROM staged WHERE error IS NULL ORDER BY source,source_row) TO {sql_path(directory / 'aliases.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)")
            con.execute(f"COPY (SELECT row_id,source,source_row,original_id,original_smiles,split,metadata,error FROM staged WHERE error IS NOT NULL ORDER BY source,source_row) TO {sql_path(directory / 'invalid.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)")
        finally:
            con.close()
    publish_unit(root, "catalog", expected, catalog)
    catalog_path = root / "units" / "catalog" / "molecules.parquet"
    count = pq.ParquetFile(catalog_path).metadata.num_rows
    shards = []
    con = connect(memory_gb, scratch or root / "scratch")
    try:
        for number in range(math.ceil(count / shard_size)):
            unit = f"shard-{number:08d}"
            shard_signature = signature([expected, unit])
            def make(directory, number=number):
                con.execute(f"COPY (SELECT *,row_number() OVER (ORDER BY molecule_id) AS local_id FROM (SELECT * FROM read_parquet({sql_path(catalog_path)}) ORDER BY molecule_id LIMIT {shard_size} OFFSET {number * shard_size})) TO {sql_path(directory / 'molecules.parquet')} (FORMAT PARQUET, COMPRESSION ZSTD)")
                write_shard(directory / "fingerprints.h5", batches(directory / "molecules.parquet"), fingerprint)
            publish_unit(root, unit, shard_signature, make)
            shards.append({"id": unit, "count": min(shard_size, count - number * shard_size), "signature": shard_signature})
    finally:
        con.close()
    manifest = {**descriptor, "signature": expected, "molecule_count": count, "alias_count": pq.ParquetFile(root / 'units/catalog/aliases.parquet').metadata.num_rows, "invalid_count": pq.ParquetFile(root / 'units/catalog/invalid.parquet').metadata.num_rows, "shards": shards}
    atomic_json(root / "manifest.json", manifest)
    return MoleculeIndex.open(root)
