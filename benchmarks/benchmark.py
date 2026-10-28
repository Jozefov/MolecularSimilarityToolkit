"""Reproducible, subprocess-isolated benchmarks; no external datasets are bundled."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import resource
import subprocess
import sys
import time
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
POP = np.array([i.bit_count() for i in range(256)], dtype=np.uint8)


def checksum(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def timed(call):
    before = time.perf_counter()
    value = call()
    return value, time.perf_counter() - before


def rss():
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value / (1024**2 if sys.platform == "darwin" else 1024)


def stage(args):
    from FPSim2 import FPSim2Engine
    from rdkit import DataStructs

    import molecular_similarity_toolkit as mst
    from molecular_similarity_toolkit.storage import batches

    out = Path(args.out)
    timings = {}
    if args.phase == "build":
        index, timings["index_build_seconds"] = timed(
            lambda: mst.build_index(
                mst.Source(out / "input.tsv", name=args.dataset),
                out / "index",
                shard_size=args.shard_size,
                workers=args.workers,
                memory_gb=args.memory_gb,
            )
        )
        selected = []
        for row in batches(index.molecules_path):
            selected.append(
                {"id": str(len(selected)), "smiles": row["canonical_smiles"]}
            )
            if len(selected) >= args.queries:
                break
        with open(out / "queries.tsv", "w") as f:
            writer = csv.DictWriter(f, fieldnames=["id", "smiles"], delimiter="\t")
            writer.writeheader()
            writer.writerows(selected)
        timings.update(
            unique_molecules=index.size,
            invalid_molecules=index.manifest["invalid_count"],
        )
    elif args.phase == "package":
        index = mst.MoleculeIndex.open(out / "index")
        queries = mst.Source(out / "queries.tsv", name="queries")
        hits, timings["search_seconds"] = timed(
            lambda: mst.search(
                index,
                queries,
                k=args.k,
                out=out / "search",
                workers=args.workers,
                memory_gb=args.memory_gb,
            )
        )
        clean, timings["filter_seconds"] = timed(
            lambda: mst.filter_pool(
                index,
                protected=queries,
                out=out / "filter",
                workers=args.workers,
                memory_gb=args.memory_gb,
            )
        )
        timings.update(search=hits.summary, filter=clean.summary)
        scores = {}
        for batch in hits.batches():
            for row in batch.to_pylist():
                scores.setdefault(row["query_original_id"], []).append(
                    row["similarity"]
                )
        timings["scores"] = scores
        timings["rejected_count"] = clean.summary["rejected_count"]
    else:
        from molecular_similarity_toolkit.chemistry import prepare

        index = mst.MoleculeIndex.open(out / "index")
        with open(out / "queries.tsv") as f:
            queries = [
                (r["id"], prepare(r["smiles"], index.chemistry, index.fingerprint)[2])
                for r in csv.DictReader(f, delimiter="\t")
            ]
        scores = {qid: [] for qid, _ in queries}
        load_seconds = search_seconds = filter_seconds = 0.0
        rejected_count = 0
        for shard in index.manifest["shards"]:
            path = index.path / "units" / shard["id"]
            if args.phase == "direct":
                engine, elapsed = timed(
                    lambda path=path: FPSim2Engine(str(path / "fingerprints.h5"))
                )
                load_seconds += elapsed
                rejected = set()
                for qid, packed in queries:
                    fp = DataStructs.CreateFromBinaryText(packed)
                    hits, elapsed = timed(
                        lambda engine=engine, fp=fp: engine.top_k(
                            fp, k=args.k, threshold=0.0, n_workers=args.workers
                        )
                    )
                    search_seconds += elapsed
                    scores[qid].extend(float(v) for v in hits["coeff"])
                    hits, elapsed = timed(
                        lambda engine=engine, fp=fp: engine.similarity(
                            fp, threshold=0.8, n_workers=args.workers
                        )
                    )
                    filter_seconds += elapsed
                    # Direct backend float32 threshold timings are reported separately.
                    # Compare counts using > 0.8 after rounding; exact-boundary oracle tests
                    # live in the test suite rather than using this baseline as an oracle.
                    rejected.update(
                        int(r["mol_id"]) for r in hits if r["coeff"] > np.float32(0.8)
                    )
                rejected_count += len(rejected)
                del engine
            else:
                start = time.perf_counter()
                packed_rows = pq.read_table(
                    path / "molecules.parquet", columns=["fingerprint"]
                )["fingerprint"].to_pylist()
                matrix = np.frombuffer(b"".join(packed_rows), dtype=np.uint8).reshape(
                    -1, index.fingerprint.n_bits // 8
                )
                del packed_rows
                load_seconds += time.perf_counter() - start
                rejected = np.zeros(len(matrix), dtype=bool)
                for qid, packed in queries:
                    query = np.frombuffer(packed, dtype=np.uint8)
                    start = time.perf_counter()
                    # Historical NumPy baseline: byte lookup + per-query reference popcount.
                    inter = POP[query & matrix].sum(axis=1)
                    similarity = inter / np.maximum(
                        int(POP[query].sum()) + POP[matrix].sum(axis=1) - inter, 1
                    )
                    top = np.partition(similarity, max(0, len(similarity) - args.k))[
                        -args.k :
                    ]
                    scores[qid].extend(float(v) for v in top)
                    search_seconds += time.perf_counter() - start
                    start = time.perf_counter()
                    rejected |= similarity > 0.8
                    filter_seconds += time.perf_counter() - start
                rejected_count += int(rejected.sum())
                del matrix
        timings.update(
            load_seconds=load_seconds,
            search_seconds=search_seconds,
            filter_seconds=filter_seconds,
            rejected_count=rejected_count,
            scores={
                qid: sorted(v, reverse=True)[: args.k] for qid, v in scores.items()
            },
        )
        if args.phase == "numpy":
            timings["filter_timing_note"] = (
                "Mask update only; similarity computation is included in search_seconds and reused."
            )
    timings["peak_rss_mib"] = rss()
    timings["phase"] = args.phase
    (out / f"{args.phase}.json").write_text(json.dumps(timings, indent=2) + "\n")


def real(args):
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(
        (out / f"{phase}.json").exists()
        for phase in ("build", "package", "direct", "numpy")
    ):
        raise SystemExit(
            "Use a fresh output directory for measured runs; resume timings would be misleading."
        )
    # The selected input is immutable and small relative to the upstream corpus.
    source = Path(args.source).resolve()
    read_rows = 0
    with (
        open(source, newline="") as src,
        open(out / "input.tsv", "w", newline="") as dst,
    ):
        reader = csv.DictReader(src, delimiter="\t" if source.suffix == ".tsv" else ",")
        writer = csv.DictWriter(dst, fieldnames=["id", "smiles"], delimiter="\t")
        writer.writeheader()
        for i, row in enumerate(reader):
            writer.writerow(
                {
                    "id": row.get(args.id_column, str(i + 1)),
                    "smiles": row[args.smiles_column],
                }
            )
            read_rows += 1
            if read_rows >= args.limit:
                break
    metadata = {
        "dataset_label": args.dataset,
        "source_basename": source.name,
        "selection": f"first {read_rows} data rows; no random sampling",
        "input_sha256": checksum(out / "input.tsv"),
        "input_rows": read_rows,
        "hardware": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu_count": os.cpu_count(),
        },
        "python": sys.version,
        "dependency_versions": {
            name: version(name)
            for name in (
                "FPSim2",
                "rdkit",
                "numpy",
                "pyarrow",
                "duckdb",
                "molecular-similarity-toolkit",
            )
        },
        "workers": args.workers,
        "memory_gb": args.memory_gb,
        "cache_condition": args.cache_condition,
        "cache_note": "OS page cache was not flushed. Build precedes search; comparisons run in the recorded order.",
        "fingerprint": {
            "type": "Morgan",
            "radius": 2,
            "n_bits": 2048,
            "use_chirality": False,
        },
        "chemistry": "stereo_free; charges, salts, tautomers and isotopes retained",
        "queries": args.queries,
        "k": args.k,
        "shard_size": args.shard_size,
        "threshold": 0.8,
        "comparison": ">",
        "phase_order": ["build", "package", "direct", "numpy"],
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    for phase in metadata["phase_order"]:
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "_phase",
            "--phase",
            phase,
            "--out",
            str(out),
            "--dataset",
            args.dataset,
            "--queries",
            str(args.queries),
            "--k",
            str(args.k),
            "--workers",
            str(args.workers),
            "--memory-gb",
            str(args.memory_gb),
            "--shard-size",
            str(args.shard_size),
        ]
        subprocess.run(command, check=True)
        print(f"Completed {phase}: {out}", flush=True)
    report = {
        **metadata,
        "measurements": [
            json.loads((out / f"{p}.json").read_text()) for p in metadata["phase_order"]
        ],
    }
    package, direct, numpy = report["measurements"][1:]
    report["checks"] = {
        "topk_scores_vs_numpy": all(
            np.allclose(v, numpy["scores"][q], rtol=0, atol=1e-14)
            for q, v in package["scores"].items()
        ),
        "topk_scores_vs_direct_float32": all(
            np.allclose(v, direct["scores"][q], rtol=0, atol=1e-7)
            for q, v in package["scores"].items()
        ),
        "rejected_count_vs_numpy": package["rejected_count"] == numpy["rejected_count"],
    }
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    if not all(report["checks"].values()):
        raise SystemExit("Benchmark correctness comparison failed")


def synthetic(args):
    """Stress the actual FPSim2 adapter with bounded, synthetic fingerprint shards."""
    from molecular_similarity_toolkit.adapter import ShardEngine, write_shard
    from molecular_similarity_toolkit.chemistry import intersection_union
    from molecular_similarity_toolkit.config import FingerprintConfig

    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    config = FingerprintConfig()
    rng = np.random.default_rng(args.seed)
    query = rng.integers(0, 256, size=256, dtype=np.uint8).tobytes()
    build_seconds = search_seconds = 0.0
    best = []
    shard_count = (args.rows + args.shard_size - 1) // args.shard_size
    for shard in range(shard_count):
        count = min(args.shard_size, args.rows - shard * args.shard_size)
        path = out / f"shard-{shard:04d}.h5"
        packed_path = out / f"shard-{shard:04d}.npy"
        if path.exists() or packed_path.exists():
            raise SystemExit("Use a fresh synthetic output directory.")
        matrix = np.lib.format.open_memmap(
            packed_path, mode="w+", dtype=np.uint8, shape=(count, 256)
        )
        for begin in range(0, count, 8192):
            matrix[begin : begin + 8192] = rng.integers(
                0, 256, size=matrix[begin : begin + 8192].shape, dtype=np.uint8
            )
        matrix.flush()
        rows = (
            {"local_id": i + 1, "fingerprint": matrix[i].tobytes()}
            for i in range(count)
        )
        _, elapsed = timed(lambda path=path, rows=rows: write_shard(path, rows, config))
        build_seconds += elapsed
        start = time.perf_counter()
        engine = ShardEngine(path)
        ids = engine.top_k_candidates(query, args.k)
        values = [
            (intersection_union(query, matrix[int(i) - 1].tobytes()), int(i))
            for i in ids
        ]
        chosen = sorted((a / b for (a, b), _i in values), reverse=True)[: args.k]
        best.extend(chosen)
        search_seconds += time.perf_counter() - start
        # Independent chunked packed-bit oracle; memory does not scale with row count.
        oracle = []
        q = np.frombuffer(query, dtype=np.uint8)
        for begin in range(0, count, 8192):
            chunk = matrix[begin : begin + 8192]
            inter = POP[q & chunk].sum(1)
            sim = inter / (POP[q].sum() + POP[chunk].sum(1) - inter)
            oracle.extend(
                np.partition(sim, max(0, len(sim) - args.k))[-args.k :].tolist()
            )
        if not np.allclose(
            chosen, sorted(oracle, reverse=True)[: args.k], rtol=0, atol=1e-14
        ):
            raise AssertionError("Synthetic adapter/NumPy mismatch")
        del engine, matrix
        print(f"Synthetic shard {shard + 1}/{shard_count} verified", flush=True)
    report = {
        "kind": "synthetic fingerprints, NOT real molecules",
        "rows": args.rows,
        "seed": args.seed,
        "distribution": "uniform random packed bytes (approximately 50% set bits)",
        "shard_size": args.shard_size,
        "shards": shard_count,
        "n_bits": 2048,
        "build_seconds": build_seconds,
        "search_seconds": search_seconds,
        "peak_rss_mib": rss(),
        "oracle_agreement": True,
        "global_topk_scores": sorted(best, reverse=True)[: args.k],
        "hardware": platform.platform(),
        "cache_condition": "uncontrolled; freshly written shards",
        "scope": "FPSim2 adapter, bounded shard storage and exact rescoring; chemistry, aliases and scheduling excluded",
    }
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")


def positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def positive_float(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be a positive finite number")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for command in ("real", "_phase"):
        p = sub.add_parser(command)
        p.add_argument("--out", required=True)
        p.add_argument("--dataset", default="user supplied dataset")
        p.add_argument("--queries", type=positive_int, default=10)
        p.add_argument("--k", type=positive_int, default=20)
        p.add_argument("--workers", type=positive_int, default=4)
        p.add_argument("--memory-gb", type=positive_float, default=4)
        p.add_argument("--shard-size", type=positive_int, default=1_000_000)
        if command == "real":
            p.add_argument("--source", required=True)
            p.add_argument("--smiles-column", default="smiles")
            p.add_argument("--id-column", default="id")
            p.add_argument("--limit", type=positive_int, default=100_000)
            p.add_argument(
                "--cache-condition",
                choices=["warm", "cold", "uncontrolled"],
                default="uncontrolled",
            )
        else:
            p.add_argument(
                "--phase",
                choices=["build", "package", "direct", "numpy"],
                required=True,
            )
    p = sub.add_parser("synthetic")
    p.add_argument("--out", required=True)
    p.add_argument("--rows", type=positive_int, default=10_000_000)
    p.add_argument("--shard-size", type=positive_int, default=1_000_000)
    p.add_argument("--seed", type=int, default=20261028)
    p.add_argument("--k", type=positive_int, default=20)
    args = parser.parse_args()
    if args.command == "real":
        real(args)
    elif args.command == "synthetic":
        synthetic(args)
    else:
        stage(args)


if __name__ == "__main__":
    main()
