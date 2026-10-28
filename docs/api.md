# Python API

```python
import molecular_similarity_toolkit as mst
```

Use file paths or index handles. Large results remain on disk. Functions raise
ordinary Python exceptions; the CLI exits nonzero on failure.

## Describe inputs

```python
source = mst.Source('pool.tsv', name='pool', id_column='cid', smiles_column='smiles')
tiny = mst.Source.from_records([{'id': 'ethanol', 'smiles': 'CCO'}], name='demo')
protected_source = mst.Source('benchmark.tsv', name='benchmark',
                              split_column='fold', splits=['val', 'test'])
```

`Source` is a frozen dataclass. Default `id_column=None` detects an `id` field
or generates row IDs; explicit column names must exist. File readers stream data; `from_records` snapshots
small records in memory. `metadata_columns` selects extra fields to retain.
Use `to_dict()` / `from_dict()` for JSON configuration. Invalid source mappings,
unsupported formats, or non-unique source names raise `ValueError`.

::: molecular_similarity_toolkit.Source

## Configure chemistry

Defaults are stereo-free Morgan radius 2 / 2,048 bits. The configuration objects
are frozen dataclasses and serialize with `to_dict()` / `from_dict()`.
`FingerprintConfig` requires a nonnegative radius and bit length divisible by
64. Chiral fingerprints cannot be combined with removed stereochemistry.

::: molecular_similarity_toolkit.ChemistryConfig

::: molecular_similarity_toolkit.FingerprintConfig

## Build an index

```python
pool = mst.build_index(source, out='indexes/pool', workers=4, memory_gb=4,
                       shard_size=1_000_000, scratch='work/spill')
```

Accepts a `Source` or iterable of sources, output directory, chemistry and
fingerprint configs. Defaults: up to 1M unique structures per shard, 8,192-row
staging chunks, four requested workers, 4 GiB memory budget, resume enabled.
The return value is a `MoleculeIndex`.

Inputs are checksummed; a rerun with the same settings validates and reuses
completed units. Different input bytes, chemistry, fingerprint, shard size, or
staging settings require a new output directory. `resume=False` refuses reuse.
Invalid source molecules are quarantined. An empty source makes a valid empty
index. Memory and CPU controls cap working units; they are not OS-enforced
process limits. Reserve native-library overhead in scheduler requests.

::: molecular_similarity_toolkit.build_index

```python
pool = mst.MoleculeIndex.open('indexes/pool')
print(pool.size, pool.aliases_path, pool.invalid_path)
for batch in pool.molecules():
    print(batch.num_rows)
```

Opening requires a completed compatible index. Missing units, an unsupported
format, or incompatible RDKit version raises an error. Full byte-by-byte
integrity validation is available through `verify`.

::: molecular_similarity_toolkit.MoleculeIndex

## Search the whole index

```python
hits = mst.search(pool, mst.Source('queries.tsv', name='queries'), k=100,
                  out='runs/candidates', exclude_identity=True, workers=4)
```

`queries` accepts sources or a compatible index. `k` is positive (default 10).
`exclude_identity` defaults to false; when true it excludes canonical identity
before truncation, not all score-1 fingerprint collisions. Fewer than k eligible
molecules produces fewer rows. Query row identities remain separate even when
structures or original IDs repeat. Invalid queries are quarantined.

`workers`, `memory_gb`, and `scratch` control resources as in index construction.
`query_chunk_size=256` bounds query task sizes. A shard too large for the requested
budget raises `MemoryError`; rebuild with a smaller shard size or raise the
budget. Incompatible indexes or changed job settings raise `ValueError`.

::: molecular_similarity_toolkit.search

## Filter a pool

```python
clean = mst.filter_pool(pool, protected=protected_source, threshold='0.8',
                        comparison='>', out='runs/clean', workers=4)
```

`protected` accepts sources or a compatible index. Every valid protected
molecule is used. Defaults reject strictly above 0.8. The `>=` rule also rejects
equality. Thresholds must lie in [0, 1]; decimal strings preserve the intended
boundary. Invalid protected molecules raise by default, and an empty protected
set always raises. An explicit `invalid_protected='quarantine'` opts into weaker
protection covering only valid reference molecules.

The return value is a `Result`; its default table contains every eligible
canonical structure. The rejected table records deterministic witnesses.
Resource arguments have the same meaning as search. No sampling occurs.

::: molecular_similarity_toolkit.filter_pool

## Read and export results

```python
for batch in hits.batches(batch_size=8192):
    print(batch.num_rows)
clean.export('rejections.tsv', table='rejected')
print(clean.summary)
```

`batches` yields Arrow record batches. `table_path` exposes a Parquet path.
`export` writes CSV/TSV based on the output suffix or explicit `delimiter`.
These methods stream data and do not create a corpus-sized pandas frame.
Unknown table names, incomplete output, or unsupported export choices fail
clearly. See [table schemas](data.md).

::: molecular_similarity_toolkit.results.Result

## Plan, run and merge cluster work

```python
plan = mst.plan_search(pool, mst.Source('queries.tsv'), 100, 'runs/cluster')
# In separate workers, use each task ID from plan.json exactly once at a time:
mst.run_task(plan, task_id=0)
# After every planned task succeeds:
# result = mst.merge(plan)
```

`plan_build`, `plan_search`, and `plan_filter` return a plan-file path. Planning
performs input preparation on the coordinator. `run_task` executes one immutable
unit; it is safe to retry after a failure. `merge` requires every task to pass
integrity checks. Build plans return an index; search/filter plans return a
result. See [scheduler examples and shared-filesystem requirements](hpc.md).

::: molecular_similarity_toolkit.plan_build

::: molecular_similarity_toolkit.plan_search

::: molecular_similarity_toolkit.plan_filter

::: molecular_similarity_toolkit.run_task

::: molecular_similarity_toolkit.merge

## Inspect progress and verify integrity

```python
print(mst.status('runs/candidates'))
report = mst.verify('runs/candidates')
assert report['ok'], report['failures']
```

`status` is a quick progress view, not a checksum guarantee. `verify` reads all
output bytes covered by receipts and checks plans; it returns a dictionary with
`ok`, failures and progress. It is read-only and can take substantial I/O on a
large corpus. Use the CLI `mol-sim verify PATH` when a nonzero exit code is needed.

::: molecular_similarity_toolkit.status

::: molecular_similarity_toolkit.verify

## Execute a whole plan locally

```python
plan = mst.plan_build(tiny, 'indexes/planned-demo')
index = mst.execute(plan)
```

`execute` runs the planned units locally and reduces them. It returns the same
index/result as a successful distributed run. Its optional worker override
controls query concurrency; memory and scratch settings come from the plan.
It propagates task and integrity failures without marking the result complete.

::: molecular_similarity_toolkit.execute

## Export a compact index

```python
compact = mst.export_index(pool, 'indexes/portable-pool')
assert mst.verify(compact.path)['ok']
```

The destination must not exist. The function verifies runtime units, copies them
in bounded streams and publishes the destination atomically. It preserves the
original build and provenance. The returned `MoleculeIndex` can be moved without
its staging files or original input files. It supports search/filter/verify;
rebuilding requires the original build. Corrupt input units raise before export.

::: molecular_similarity_toolkit.export_index

## Reopen and export a completed run

```python
saved = mst.Result.open('runs/candidates')
saved.export('candidates.tsv')
```

`Result.open()` verifies its published result-table receipt. Use `verify()` to
include all tasks and source dependencies. From the command line:

```bash
mol-sim export runs/candidates candidates.tsv
mol-sim export runs/clean rejected.tsv --table rejected
```

To read the full generated API locally, run `uv sync --locked --group dev --extra
docs`, then `uv run mkdocs serve`. CI also builds the HTML reference into the
private `private-release-candidate` artifact. Source Markdown on GitHub contains
API-generation directives; the built HTML contains the expanded signatures.
