"""The only dependency boundary to pinned FPSim2 0.7.4 (including its file writer)."""

from pathlib import Path

import numpy as np
import tables as tb
from FPSim2 import FPSim2Engine, __version__
from FPSim2.io.backends.pytables import create_schema, sort_db_file
from FPSim2.io.chem import process_fp
from rdkit import DataStructs

from .chemistry import RDKIT_VERSION
from .config import FingerprintConfig


def write_shard(path: Path, rows, config: FingerprintConfig):
    """Create FPSim2 HDF5 from staged fingerprints, without recomputation."""
    params = {
        "radius": config.radius,
        "fpSize": config.n_bits,
        "includeChirality": config.use_chirality,
    }
    with tb.open_file(str(path), "w") as handle:
        table = handle.create_table(
            handle.root,
            "fps",
            create_schema(config.n_bits),
            filters=tb.Filters(complib="blosc2", complevel=5),
        )
        settings = handle.create_vlarray(handle.root, "config", atom=tb.ObjectAtom())
        for value in ["Morgan", params, RDKIT_VERSION, __version__]:
            settings.append(value)
        pending = []
        for row in rows:
            local_id = int(row["local_id"])
            if not 0 < local_id <= 1_000_000:
                raise ValueError("FPSim2 IDs must be shard-local integers <= 1,000,000")
            pending.append(
                process_fp(
                    DataStructs.CreateFromBinaryText(row["fingerprint"]), local_id
                )
            )
            if len(pending) == 8192:
                table.append(pending)
                pending = []
        if pending:
            table.append(pending)
        table.cols.popcnt.create_index(kind="full")
    sort_db_file(str(path))


def conservative_threshold(value: float) -> float:
    """Lower by two float32 ULPs to include every rounded boundary candidate."""
    number = np.float32(value)
    return max(
        0.0,
        float(
            np.nextafter(np.nextafter(number, np.float32(-np.inf)), np.float32(-np.inf))
        ),
    )


class ShardEngine:
    """One resident read-only engine; query-level threads share its fingerprint array."""

    def __init__(self, path: Path):
        self.engine = FPSim2Engine(str(path), in_memory_fps=True)

    def top_k_candidates(self, fp: bytes, k: int):
        query = DataStructs.CreateFromBinaryText(fp)
        initial = self.engine.top_k(query, k=k, threshold=0.0, n_workers=1)
        # FPSim2 returns float32 scores and may choose arbitrary boundary ties.
        lower = (
            conservative_threshold(float(initial["coeff"].min()))
            if len(initial)
            else 0.0
        )
        return self.engine.similarity(query, threshold=lower, n_workers=1)["mol_id"]

    def threshold_candidates(self, fp: bytes, threshold: float):
        query = DataStructs.CreateFromBinaryText(fp)
        return self.engine.similarity(
            query, threshold=conservative_threshold(threshold), n_workers=1
        )["mol_id"]
