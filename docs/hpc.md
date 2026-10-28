# Cluster jobs

Local runs and cluster jobs use the same task files. The package prepares plans;
it never submits jobs to a scheduler. All nodes need the same package environment
and access to the recorded absolute paths on a shared filesystem.

## Prepare inputs once

Describe datasets in JSON. Relative input paths are resolved against the JSON
file's directory, not the shell's working directory.

```json
{"sources": [{"path": "benchmark.tsv", "name": "benchmark",
              "id_column": "id", "smiles_column": "smiles",
              "split_column": "fold", "splits": ["val", "test"]}]}
```

```bash
mol-sim build examples/sources.json indexes/pool --workers 4 --memory-gb 4
mol-sim build examples/protected.json indexes/protected
mol-sim plan filter indexes/pool indexes/protected runs/cluster --threshold 0.8 --comparison '>'
```

For candidate search:

```bash
mol-sim plan search indexes/pool examples/queries.tsv runs/cluster-search --k 100
```

The resulting `plan.json` lists numeric task IDs. Inspect it before submitting.
Filtering has one task per pool shard; search tasks cover query chunks and pool
shards. The planner stages queries and validates inputs on the coordinator.

Index construction can also distribute fingerprint-file creation:

```bash
mol-sim plan build examples/sources.json indexes/cluster-pool
```

This performs bounded preparation and global deduplication first, then writes
`build-plan.json` for independent immutable HDF5 shard builds. The coordinator
handles input staging; the cluster does not write a shared DuckDB database.

## Slurm or PBS

The repository includes `examples/slurm-array.sh` and `examples/pbs-array.sh`.
Set an array range from zero to `len(plan['tasks']) - 1`. For an empty task list,
skip array submission and run the merge directly. Resource and time requests
are examples; size them from measured data.

```bash
export MST_PLAN=/shared/project/runs/cluster/plan.json
# Example only, if the plan contains exactly eight tasks:
# sbatch --array=0-7 --export=ALL examples/slurm-array.sh
# qsub -J 0-7 -v MST_PLAN="$MST_PLAN" examples/pbs-array.sh
```

The scripts execute:

```bash
mol-sim task /shared/project/runs/cluster/plan.json 0 --workers 4
```

One task writes only its own output and receipt. Do not run duplicate copies
of the same task simultaneously. After a task has stopped, a retry validates
and reuses completed data or repairs the incomplete unit.

## Reduce only after all tasks finish

```bash
mol-sim status /shared/project/runs/cluster
mol-sim merge /shared/project/runs/cluster/plan.json
mol-sim verify /shared/project/runs/cluster
```

A merge refuses missing or damaged task outputs. Array completion alone is not
a scientific correctness check. Keep the output and receipt directories.

CPU limits honor scheduler variables where available. The package avoids nested
worker pools and shares a resident shard among query threads. More CPUs do not
guarantee proportional speedup; memory bandwidth and output I/O can dominate.
