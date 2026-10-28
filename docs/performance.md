# Performance and reproducibility

Packed fingerprints save space: 2,048 bits occupy 256 raw bytes per molecule.
That is 25.6 GB for 100 million raw fingerprints, before IDs, catalogs, staging,
indexes, compression, temporary files, and runtime overhead. It is not a total
storage estimate.

FPSim2 uses packed machine words, CPU bit counting and mathematically safe
similarity bounds. High thresholds can prune work. Exact top-k needs boundary
tie handling, and large outputs can cost more than the search itself.

## Run the benchmark

Use a fresh output directory so restart validation is not mistaken for build
speed. Run only one benchmark at a time. The harness makes a fixed first-N-row
subset and records its checksum, dataset label, hardware, worker count, settings,
phase times and peak RSS in separate subprocesses.

```bash
uv run python benchmarks/benchmark.py real --source molecules.tsv \
  --smiles-column smiles --id-column id --dataset 'my dated dataset snapshot' \
  --limit 100000 --queries 10 --workers 4 --out benchmarks/runs/real-100k

uv run python benchmarks/benchmark.py real --source molecules.tsv \
  --dataset 'my dated dataset snapshot' --limit 1000000 --queries 10 \
  --workers 4 --out benchmarks/runs/real-1m

uv run python benchmarks/benchmark.py synthetic --rows 10000000 \
  --out benchmarks/runs/synthetic-10m
```

The real-data harness compares historical NumPy byte-popcount scanning, direct
FPSim2, and the full package workflow. The package includes query preparation,
checksum validation, deterministic ranking/witnesses, task publication, and
Parquet output. The raw backend does not provide all those guarantees by itself.

NumPy's mask-update timing reuses similarities computed during search; the
report labels this and it must not be compared with a complete independent
filter run. Direct FPSim2 scores are float32; its timing is a baseline, not the
independent correctness oracle. Tests use exhaustive RDKit comparisons, and the
benchmark compares all query IDs, ordered match IDs, integer scores and rejected
molecule-ID sets against packed NumPy and an independent direct-backend path.

The synthetic test uses uniform random fingerprint bytes, not chemical
structures. It exercises 10M rows through the actual FPSim2 adapter in bounded
shards, verifies exact scores against a chunked NumPy reference, and records
peak RSS. It does not measure molecular normalization, provenance joins, or a
real-world distribution of chemical fingerprints.

## Cache and memory interpretation

The harness never flushes a machine-wide filesystem cache. Default reports say
`uncontrolled`; freshly built files are likely cached. Record `warm` or `cold`
only when you have deliberately arranged and documented that condition. For a
controlled cold-cache measurement use an isolated host or fresh volume and
record how you did it. Do not treat an uncontrolled run as a cold-cache result.

The parent OS high-water mark and sampled process-tree RSS are both recorded.
Process-tree sampling includes preparation workers. Total disk space and memory
depend on provenance size, duplicates, query
counts, cutoff ties and native-library overhead. Shard sizes and worker limits
bound resident working sets; they are not a hard process memory sandbox.

## Measured evidence

Results and interpretation from the private release review are recorded in
[release review](release.md). These are local measurements, not a whole-PubChem
speed claim. Full-corpus claims require a separately measured run on the stated
snapshot and hardware.

Upstream: [FPSim2 search documentation](https://chembl.github.io/FPSim2/user_guide/similarity/)
and [FPSim2 implementation](https://github.com/chembl/FPSim2).

## Version 0.2 validation and larger query sets

Each benchmark now compares the **complete query-ID set**, ordered matched
molecule IDs and integer intersection/union counts. Filtering compares rejected
molecule-ID sets, not only counts. Missing queries or swapped identities cannot
pass the validation flags.

Process-tree RSS is sampled every 50 ms and includes preparation subprocesses;
the report also retains the parent process's OS high-water mark. Sampling may
miss very short peaks, while summed process RSS may count shared pages in more
than one process. Preparation CPU and wall times and worker PIDs are recorded.

Use `--queries 1000 --shard-size 250000 --skip-numpy` for a larger, four-shard
1M-molecule run. This validates every result against direct FPSim2 with independent
RDKit integer rescoring, while omitting the much slower NumPy baseline. The 100K
and 1M ten-query runs still exercise all three implementations. Reports state
which baselines ran. Raw FPSim2 kernel timings exclude exact boundary expansion
and validation I/O; full package timings include them.

The direct baseline uses **one native kernel worker**, matching the package's
serial native kernels. The package can still process several queries at once.
A preliminary four-shard direct benchmark with four internal FPSim2 workers
segfaulted once; an identical diagnostic rerun passed. The underlying cause was
not established. The supported package does not use that internal parallel path,
and the benchmark now uses the same serial-kernel mode explicitly. This is
recorded as a benchmark limitation, not silently counted as a passing run.
