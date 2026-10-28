"""Resource caps for one outer pool and serial FPSim2 kernels."""

import os


def effective_workers(
    workers: int = 4, memory_gb: float = 4, *, per_worker_gb: float = 0.25
) -> int:
    """Cap requested workers by CPU affinity, common scheduler limits and memory."""
    if workers < 1 or memory_gb <= 0:
        raise ValueError("workers and memory_gb must be positive")
    cpu = (
        len(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else (os.cpu_count() or 1)
    )
    for key in ("SLURM_CPUS_PER_TASK", "PBS_NP", "NSLOTS"):
        if os.environ.get(key):
            try:
                cpu = min(cpu, int(os.environ[key]))
            except ValueError:
                pass
    return max(1, min(workers, cpu, int(memory_gb / per_worker_gb)))
