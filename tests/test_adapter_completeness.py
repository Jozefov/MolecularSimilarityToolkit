from fractions import Fraction
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import DataStructs

from molecular_similarity_toolkit.adapter import ShardEngine, write_shard
from molecular_similarity_toolkit.config import FingerprintConfig


def test_random_and_constructed_fingerprint_candidates(tmp_path):
    rng = np.random.default_rng(319)
    checks = 0
    with TemporaryDirectory(dir=tmp_path) as temp:
        for bits in [64, 128, 2048, 4096]:
            fps = []
            for density in [0.01, 0.1, 0.5, 0.9, 1]:
                for _ in range(30):
                    b = rng.random(bits) < density
                    b[0] = True
                    fps.append(np.packbits(b, bitorder="little").tobytes())
            # An exact 0.8 pair and numerous exact ties.
            for n in [20, 25, 30, 40, 50]:
                b = np.zeros(bits, dtype=np.uint8)
                b[:n] = 1
                fps.append(np.packbits(b, bitorder="little").tobytes())
            path = Path(temp) / f"{bits}.h5"
            write_shard(
                path,
                ({"local_id": i + 1, "fingerprint": f} for i, f in enumerate(fps)),
                FingerprintConfig(n_bits=bits),
            )
            engine = ShardEngine(path)
            for query in fps[::11]:
                exact = {
                    i + 1: Fraction(
                        (
                            DataStructs.CreateFromBinaryText(query)
                            & DataStructs.CreateFromBinaryText(f)
                        ).GetNumOnBits(),
                        (
                            DataStructs.CreateFromBinaryText(query)
                            | DataStructs.CreateFromBinaryText(f)
                        ).GetNumOnBits(),
                    )
                    for i, f in enumerate(fps)
                }
                for k in [1, 5, 20, len(fps)]:
                    candidates = set(map(int, engine.top_k_candidates(query, k)))
                    boundary = sorted(exact.values(), reverse=True)[k - 1]
                    needed = {i for i, v in exact.items() if v >= boundary}
                    assert needed <= candidates, (bits, k, needed - candidates)
                    checks += 1
                for t in [
                    0,
                    0.00001,
                    0.1,
                    0.5,
                    0.799999999999,
                    0.8,
                    0.800000000001,
                    0.99,
                    1,
                ]:
                    candidates = set(map(int, engine.threshold_candidates(query, t)))
                    needed = {i for i, v in exact.items() if v >= Fraction(str(t))}
                    assert needed <= candidates, (bits, t, needed - candidates)
                    checks += 1
    assert checks == 780
