"""Run from a checkout: python examples/quickstart.py [OUTPUT_DIRECTORY]."""

import sys
from pathlib import Path

import molecular_similarity_toolkit as mst

HERE = Path(__file__).resolve().parent
out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("runs/quickstart")
pool = mst.build_index(
    mst.Source(HERE / "pool.tsv", name="demo", metadata_columns=["origin"]),
    out / "pool-index",
    shard_size=3,
)
protected = mst.build_index(
    mst.Source(
        HERE / "benchmark.tsv",
        name="benchmark",
        split_column="fold",
        splits=["val", "test"],
    ),
    out / "protected-index",
)
clean = mst.filter_pool(
    pool, protected=protected, threshold=0.8, comparison=">", out=out / "clean"
)
hits = mst.search(
    pool,
    queries=mst.Source(HERE / "queries.tsv", name="queries"),
    k=3,
    out=out / "candidates",
)
clean.export(out / "eligible.tsv")
hits.export(out / "candidates.tsv")
print("Pool:", pool.size, "unique molecules")
print("Filter:", clean.summary)
print("Search:", hits.summary)
print("Results:", out.resolve())
