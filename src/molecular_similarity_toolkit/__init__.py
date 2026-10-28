"""Exact, disk-backed molecular similarity and protected-set filtering."""

from .config import ChemistryConfig, FingerprintConfig, Source
from .index import MoleculeIndex, build_index, plan_build
from .inspection import status, verify
from .operations import (
    execute,
    filter_pool,
    merge,
    plan_filter,
    plan_search,
    run_task,
    search,
)
from .results import Result

__version__ = "0.1.0"
__all__ = [
    "ChemistryConfig",
    "FingerprintConfig",
    "MoleculeIndex",
    "Result",
    "Source",
    "build_index",
    "execute",
    "filter_pool",
    "merge",
    "plan_build",
    "plan_filter",
    "plan_search",
    "run_task",
    "search",
    "status",
    "verify",
]
