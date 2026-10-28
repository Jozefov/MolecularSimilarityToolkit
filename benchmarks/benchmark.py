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
import psutil
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


def direct_top_k(engine, query, k):
    """Match the package's serial native-kernel mode, outside its adapter."""
    return engine.top_k(query, k=k, threshold=0.0, n_workers=1)


def direct_threshold(engine, query, threshold):
    return engine.similarity(query, threshold=threshold, n_workers=1)


def _write_validation(out, phase, matches, rejected):
    import pyarrow as pa
    import pyarrow.parquet as pq

    schema = pa.schema(
        [
            ("query_id", pa.string()),
            ("molecule_id", pa.string()),
            ("rank", pa.int64()),
            ("intersection", pa.int64()),
            ("union", pa.int64()),
        ]
    )
    rows = []
    for qid, candidates in sorted(matches.items()):
        for rank, (mol_id, a, b) in enumerate(candidates, 1):
            rows.append(
                {
                    "query_id": qid,
                    "molecule_id": mol_id,
                    "rank": rank,
                    "intersection": a,
                    "union": b,
                }
            )
    pq.write_table(
        pa.Table.from_pylist(rows, schema=schema), out / f"{phase}-matches.parquet"
    )
    pq.write_table(
        pa.table(
            {"molecule_id": sorted(rejected)},
            schema=pa.schema([("molecule_id", pa.string())]),
        ),
        out / f"{phase}-rejected.parquet",
    )


def compare_results(out, phases, k):
    """Compare coverage, exact ordered identities/counts, and rejection identity sets."""
    import duckdb

    out = Path(out)
    with open(out / "queries.tsv") as f:
        expected = {r["id"] for r in csv.DictReader(f, delimiter="\t")}
    if not expected:
        raise ValueError("benchmark requires valid queries")
    con = duckdb.connect()

    def relation(phase, table):
        return (
            "read_parquet('"
            + str(out / f"{phase}-{table}.parquet").replace("'", "''")
            + "')"
        )

    checks = {}
    for phase in phases:
        table = relation(phase, "matches")
        coverage = dict(
            con.execute(
                f"SELECT query_id,count(*) FROM {table} GROUP BY query_id"
            ).fetchall()
        )
        checks[f"{phase}_query_coverage"] = set(coverage) == expected and all(
            n == k for n in coverage.values()
        )
        if phase == "package":
            continue
        for table in ["matches", "rejected"]:
            a, b = relation("package", table), relation(phase, table)
            difference = con.execute(
                f"SELECT count(*) FROM ((SELECT * FROM {a} EXCEPT ALL SELECT * FROM {b}) UNION ALL (SELECT * FROM {b} EXCEPT ALL SELECT * FROM {a}))"
            ).fetchone()[0]
            checks[f"{phase}_{table}_exact_identities"] = difference == 0
    con.close()
    return checks


def stage(args):
    from fractions import Fraction

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
        if not selected:
            raise ValueError("benchmark requires valid molecules")
        with open(out / "queries.tsv", "w") as f:
            writer = csv.DictWriter(f, fieldnames=["id", "smiles"], delimiter="\t")
            writer.writeheader()
            writer.writerows(selected)
        timings.update(
            unique_molecules=index.size,
            invalid_molecules=index.manifest["invalid_count"],
        )
        progress = json.loads((index.path / "progress.json").read_text())
        timings.update(
            preparation_cpu_seconds=progress["preparation_cpu_seconds"],
            worker_pids=progress["worker_pids"],
            process_workers=progress["process_workers"],
            preparation_wall_seconds=progress["preparation_wall_seconds"],
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
        con = __import__("duckdb").connect()
        con.execute(
            f"COPY (SELECT query_original_id AS query_id,molecule_id,rank,intersection,\"union\" FROM read_parquet('{hits.table_path()}')) TO '{out / 'package-matches.parquet'}' (FORMAT PARQUET)"
        )
        con.execute(
            f"COPY (SELECT molecule_id FROM read_parquet('{clean.table_path('rejected')}')) TO '{out / 'package-rejected.parquet'}' (FORMAT PARQUET)"
        )
        con.close()
        timings.update(
            search=hits.summary,
            filter=clean.summary,
            rejected_count=clean.summary["rejected_count"],
        )
    else:
        from molecular_similarity_toolkit.chemistry import prepare

        index = mst.MoleculeIndex.open(out / "index")
        with open(out / "queries.tsv") as f:
            queries = [
                (r["id"], prepare(r["smiles"], index.chemistry, index.fingerprint)[2])
                for r in csv.DictReader(f, delimiter="\t")
            ]
        matches = {qid: [] for qid, _ in queries}
        rejected = set()
        load_seconds = search_seconds = filter_seconds = 0.0
        for shard in index.manifest["shards"]:
            path = index.path / "units" / shard["id"]
            start = time.perf_counter()
            table = pq.read_table(
                path / "molecules.parquet", columns=["molecule_id", "fingerprint"]
            )
            ids = table["molecule_id"].to_pylist()
            packed = table["fingerprint"].to_pylist()
            matrix = np.frombuffer(b"".join(packed), dtype=np.uint8).reshape(
                -1, index.fingerprint.n_bits // 8
            )
            load_seconds += time.perf_counter() - start
            if args.phase == "direct":
                engine, elapsed = timed(
                    lambda path=path: FPSim2Engine(str(path / "fingerprints.h5"))
                )
                load_seconds += elapsed
                for qid, qbytes in queries:
                    q = DataStructs.CreateFromBinaryText(qbytes)
                    initial, elapsed = timed(
                        lambda q=q, engine=engine: direct_top_k(engine, q, args.k)
                    )
                    search_seconds += elapsed
                    boundary = max(
                        0.0,
                        float(
                            np.nextafter(
                                np.nextafter(
                                    np.float32(initial["coeff"].min()),
                                    np.float32(-np.inf),
                                )
                                if len(initial)
                                else np.float32(0),
                                np.float32(-np.inf),
                            )
                        ),
                    )
                    candidates = direct_threshold(engine, q, boundary)
                    local = []
                    for row in candidates:
                        i = int(row["mol_id"]) - 1
                        fp = DataStructs.CreateFromBinaryText(packed[i])
                        a = (q & fp).GetNumOnBits()
                        b = (q | fp).GetNumOnBits()
                        local.append((ids[i], a, b))
                    local.sort(key=lambda r: (-Fraction(r[1], r[2]), r[0]))
                    matches[qid].extend(local[: args.k])
                    hits, elapsed = timed(
                        lambda q=q, engine=engine: direct_threshold(
                            engine,
                            q,
                            float(np.nextafter(np.float32(0.8), np.float32(-np.inf))),
                        )
                    )
                    filter_seconds += elapsed
                    for row in hits:
                        i = int(row["mol_id"]) - 1
                        fp = DataStructs.CreateFromBinaryText(packed[i])
                        if (q & fp).GetNumOnBits() * 5 > (q | fp).GetNumOnBits() * 4:
                            rejected.add(ids[i])
                del engine
            else:
                for qid, qbytes in queries:
                    q = np.frombuffer(qbytes, dtype=np.uint8)
                    start = time.perf_counter()
                    inter = POP[q & matrix].sum(1)
                    union = int(POP[q].sum()) + POP[matrix].sum(1) - inter
                    sim = inter / union
                    boundary = np.partition(sim, max(0, len(sim) - args.k))[
                        max(0, len(sim) - args.k)
                    ]
                    candidates = np.flatnonzero(sim >= boundary)
                    chosen = sorted(candidates, key=lambda i: (-sim[i], ids[i]))[
                        : args.k
                    ]
                    matches[qid].extend(
                        (ids[i], int(inter[i]), int(union[i])) for i in chosen
                    )
                    search_seconds += time.perf_counter() - start
                    start = time.perf_counter()
                    rejected.update(
                        ids[i] for i in np.flatnonzero(inter * 5 > union * 4)
                    )
                    filter_seconds += time.perf_counter() - start
            del matrix, table, packed, ids
        for qid, values in matches.items():
            values.sort(key=lambda r: (-Fraction(r[1], r[2]), r[0]))
            matches[qid] = values[: args.k]
        _write_validation(out, args.phase, matches, rejected)
        timings.update(
            load_seconds=load_seconds,
            search_seconds=search_seconds,
            filter_seconds=filter_seconds,
            rejected_count=len(rejected),
        )
        timings["timing_note"] = (
            "Serial native initial top-k/threshold kernels only (n_workers=1, matching package kernels); exact boundary expansion, identity/count verification and output are outside kernel timings."
            if args.phase == "direct"
            else "Historical per-query NumPy byte-popcount scan; filter mask reuses computed similarities. Stable identity ties included in search time."
        )
    timings.update(peak_rss_mib=rss(), phase=args.phase)
    (out / f"{args.phase}.json").write_text(json.dumps(timings, indent=2) + "\n")


def real(args):
    out = Path(args.out).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(
            "Use a fresh output directory; resumed timings would be misleading"
        )
    out.mkdir(parents=True, exist_ok=True)
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
        "direct_kernel_workers": 1,
        "benchmark_sha256": checksum(Path(__file__)),
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
        "phase_order": ["build", "package", "direct"]
        + ([] if args.skip_numpy else ["numpy"]),
        "query_selection": "First unique molecules in stable-ID order; varies with subset size",
        "source_checksums": {
            p.name: checksum(p)
            for p in sorted((ROOT / "src/molecular_similarity_toolkit").glob("*.py"))
        },
        "code_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
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
        process = subprocess.Popen(command)
        peak_tree = 0
        while process.poll() is None:
            try:
                parent = psutil.Process(process.pid)
                total = 0
                for child in [parent, *parent.children(recursive=True)]:
                    try:
                        total += child.memory_info().rss
                    except psutil.NoSuchProcess:
                        pass
                peak_tree = max(peak_tree, total)
            except psutil.NoSuchProcess:
                pass
            time.sleep(0.05)
        if process.returncode:
            raise subprocess.CalledProcessError(process.returncode, command)
        phase_path = out / f"{phase}.json"
        measurement = json.loads(phase_path.read_text())
        measurement["peak_process_tree_rss_mib"] = peak_tree / 2**20
        measurement["memory_note"] = (
            "Process-tree RSS sampled every 50 ms; includes preparation subprocesses. peak_rss_mib is the parent OS high-water mark."
        )
        phase_path.write_text(json.dumps(measurement, indent=2) + "\n")
        print(f"Completed {phase}: {out}", flush=True)
    report = {
        **metadata,
        "measurements": [
            json.loads((out / f"{p}.json").read_text()) for p in metadata["phase_order"]
        ],
    }
    import molecular_similarity_toolkit as mst

    index = mst.MoleculeIndex.open(out / "index")
    report["checks"] = compare_results(
        out, metadata["phase_order"][1:], min(args.k, index.size)
    )
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
        "dependency_versions": {
            name: version(name)
            for name in ("molecular-similarity-toolkit", "rdkit", "FPSim2", "numpy")
        },
        "source_checksums": {
            p.name: checksum(p)
            for p in sorted((ROOT / "src/molecular_similarity_toolkit").glob("*.py"))
        },
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
        "benchmark_sha256": checksum(Path(__file__)),
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
            p.add_argument(
                "--skip-numpy",
                action="store_true",
                help="Large-query run: validate every identity against direct FPSim2; omit the slower NumPy baseline",
            )
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
