"""Exact, disk-backed molecular similarity and protected-set filtering."""

from .config import ChemistryConfig, FingerprintConfig, Source
from .index import MoleculeIndex, build_index, export_index, plan_build
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
from .versions import __version__

__all__ = [
    "ChemistryConfig",
    "FingerprintConfig",
    "MoleculeIndex",
    "Result",
    "Source",
    "__version__",
    "build_index",
    "execute",
    "export_index",
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
