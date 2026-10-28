# Restart and integrity

## Retry the same operation

Run the same Python call or CLI command with the same output directory. Completed
units are checked and reused. Missing or damaged units are rebuilt. Changed
inputs or scientific settings require a different directory.

```bash
mol-sim status runs/candidates
mol-sim verify runs/candidates
mol-sim search indexes/pool queries.tsv runs/candidates --k 100
```

## What gets saved

Input bytes and resolved settings identify a build. Staging saves normalization
and fingerprints in Parquet chunks. A coordinator performs disk-backed
canonical deduplication and writes the catalog and source memberships. Independent
shard tasks write their own HDF5 files. Search and filtering use separate task
outputs.

Each unit is written in a temporary directory. The directory is published with
an atomic rename, then a checksum receipt is published. A crash between these
steps leaves an output without a valid receipt; retry rebuilds that unit. It
does not append another copy of its rows. A reduction is published the same way.

Gzip input may need to be replayed to find a later chunk. Completed chemistry
work is still reused. Indexes are immutable; do not replace source files in
place and expect the same index identity.

## Integrity checks

`status` counts receipts and reports progress. It is intentionally cheaper than
reading every result byte. `verify` validates output checksums and required
receipts against the plan. Scientific jobs validate the data they consume.

A reducer refuses missing, damaged, unexpected, or incompatible task receipts.
An interrupted or damaged result is not a valid completed scientific output.
Keep the entire index/run directory, including receipts, rather than copying
only an HDF5 file.

## Filesystem and concurrency

All workers must see the same index and run paths. Atomic rename must work on
the destination filesystem. Use a normal local or shared POSIX filesystem;
object-store URLs are not supported as live run directories. Workers write
separate task units and never share an HDF5 writer. Preparation and reduction
have a single coordinator. Do not launch duplicate copies of the same task
concurrently; retries should happen after the earlier attempt has stopped.

The scratch directory must have enough space for DuckDB spills. On a cluster,
select a shared scratch path when workers need it, or pass a worker-local scratch
configuration appropriate to the operation. The run's durable outputs must
survive the compute job.

## Resources

The default memory budget is 4 GiB, and workers are capped by the requested CPU
count, scheduler allocation and estimated memory. No query-by-corpus matrix is
built. Fingerprints stay packed. A smaller shard size lowers peak resident
memory at the cost of more files and scheduling overhead.

The configured budget is a planning bound, not a hard OS memory limit; Python,
RDKit, Arrow and native libraries have overhead. Set scheduler memory above the
measured peak for your workload. See [performance evidence](performance.md).
