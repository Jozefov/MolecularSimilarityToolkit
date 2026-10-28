# Input and output contract

## Source files

`Source` reads CSV, TSV, Parquet, SMILES (`.smi` / `.smiles`), or SDF (`.sdf` /
`.sd`). Text and SDF can be gzip-compressed. CSV/TSV need a header; `.smi` uses
`SMILES [identifier]` per data line. SDF uses its structure and properties; the
record name is the fallback ID. Comments and blank lines in `.smi` are skipped.

```python
source = mst.Source(
    'benchmark.tsv', name='benchmark',
    id_column='compound_id', smiles_column='canonical_smiles',
    split_column='fold', splits=['val', 'test'],
    metadata_columns=['instrument', 'license'],
)
```

Column mappings are explicit and case-sensitive. Source names must be unique
within a build. Source-row numbers are one-based **data-record** positions,
not physical line numbers in multi-line CSV or SDF files. Missing IDs get a
source-row identity; duplicate original IDs remain distinguishable.

`Source.from_records(...)` snapshots an iterable into memory. It is intended
for tiny examples and interactive use; pass file sources for large collections.

## Choosing protected folds

The quickstarts select `['val', 'test']`. A source without a fold selector uses
all of its rows. No dataset-specific fold naming is guessed. Use `splits=['test']`
to protect only test, or build a protected index from multiple named sources.

## Index layout

An index has JSON configuration/manifests, staged Parquet chunks, a deduplicated
molecule catalog, an alias table, invalid-input records, immutable HDF5 fingerprint
shards, and completion receipts. Keep the entire directory together.

| Table | Meaning |
|---|---|
| Molecules | One row per canonical structure; stable ID, canonical SMILES, packed fingerprint and bit count |
| Aliases | Every valid original source row; source, source row, original ID/SMILES, selected split/metadata, molecule ID |
| Invalid | Quarantined original rows and reasons |

Selected metadata is stored as JSON text, preserving values without expanding
arbitrary source columns into a shared wide schema.

```python
pool = mst.MoleculeIndex.open('indexes/pool')
for batch in pool.molecules(batch_size=8192):
    print(batch.num_rows)
```

`pool.aliases_path` and `pool.invalid_path` are Parquet paths for Arrow or DuckDB.

## Search results

The default table is `matches`:

| Column | Meaning |
|---|---|
| `query_id` | Unique query source-row ID, including repeated original IDs |
| `query_original_id` | Original user-provided query ID |
| `query_molecule_id` | Canonical query structure ID |
| `molecule_id` | Matched indexed structure ID |
| `rank` | One-based rank within the query |
| `similarity` | Tanimoto as a float64 display value |
| `intersection`, `union` | Exact integer counts behind the score |

`invalid` contains quarantined query rows. A query with no matches has no rows
in `matches`. `summary` records counts. Pool source memberships are available
from the indexed alias table, without multiplying ranked hits by duplicate
source memberships.

## Filtering results

- `eligible`: unique kept molecule IDs, canonical SMILES and bit counts.
- `rejected`: molecule ID, `protected_id`, score, intersection and union.
- `aliases`: original pool rows plus `disposition` (`eligible` or `rejected`).
- `invalid`: quarantined pool input rows.
- `summary.json`: input, kept, rejected and invalid counts, cutoff and rule.

`protected_id` is a deterministic witness: the smallest stable protected ID
that causes rejection. It is not necessarily the nearest protected molecule.
The protected index retains that witness's aliases and invalid-input report.

## Stream or export

```python
for batch in hits.batches(batch_size=8192):
    print(batch.num_rows)
clean.export('rejected.tsv', table='rejected')
hits.export('matches.csv')
```

Parquet is the primary storage format. CSV/TSV export is explicit and streams
batches. Avoid converting an entire PubChem-sized result to pandas at once.
