"""Immutable chemistry, fingerprint and source descriptions."""

from collections.abc import Iterable, Mapping
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class ChemistryConfig:
    """Structure identity policy. Defaults retain salts, charges, tautomers and isotopes."""

    remove_stereochemistry: bool = True
    remove_isotopes: bool = False
    normalized_parent: bool = False
    canonical_tautomer: bool = False

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict) -> "ChemistryConfig":
        return cls(**value)


@dataclass(frozen=True)
class FingerprintConfig:
    """Explicit Morgan binary fingerprint settings (unlike FPSim2 defaults)."""

    radius: int = 2
    n_bits: int = 2048
    use_chirality: bool = False

    def __post_init__(self):
        if self.radius < 0 or self.n_bits < 64 or self.n_bits % 64:
            raise ValueError(
                "radius must be nonnegative; n_bits must be a positive multiple of 64"
            )

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict) -> "FingerprintConfig":
        return cls(**value)


def validate_configs(chemistry: ChemistryConfig, fingerprint: FingerprintConfig):
    if chemistry.remove_stereochemistry and fingerprint.use_chirality:
        raise ValueError("chirality fingerprinting requires retained stereochemistry")


@dataclass(frozen=True)
class Source:
    """One input with explicit column mapping and optional split selection.

    Paths accept CSV, TSV, Parquet, SMI and SDF, including gzipped text/SDF.
    Records are snapshotted by ``from_records``; selected metadata is JSON encoded
    in alias tables. Source rows are one-based data rows (excluding CSV headers).
    """

    path: str | Path | None = None
    name: str | None = None
    id_column: str = "id"
    smiles_column: str = "smiles"
    split_column: str | None = None
    splits: tuple[str, ...] | None = None
    metadata_columns: tuple[str, ...] = ()
    _records: tuple[dict[str, Any], ...] | None = None

    def __post_init__(self):
        if (self.path is None) == (self._records is None):
            raise ValueError("provide either path or records")
        if self.splits is not None and self.split_column is None:
            raise ValueError("splits requires split_column")
        object.__setattr__(
            self, "splits", None if self.splits is None else tuple(self.splits)
        )
        object.__setattr__(self, "metadata_columns", tuple(self.metadata_columns))

    @classmethod
    def from_records(
        cls, records: Iterable[Mapping[str, Any]], *, name: str = "records", **kwargs
    ) -> "Source":
        """Snapshot an in-memory iterable. Prefer file sources for large datasets."""
        return cls(name=name, _records=tuple(dict(row) for row in records), **kwargs)

    def to_dict(self) -> dict:
        value = asdict(self)
        if self.path is not None:
            value["path"] = str(Path(self.path).resolve())
        return value

    @classmethod
    def from_dict(cls, value: dict) -> "Source":
        return cls(**value)
