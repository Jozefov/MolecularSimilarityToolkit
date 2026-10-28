"""Scientific compatibility versions; bump when preparation or results change."""

__version__ = "0.2.0"
FORMAT_VERSION = 2
PREPARATION_VERSION = 2
OPERATION_VERSION = 2


def check_preparation(value):
    if (
        value.get("format_version") != FORMAT_VERSION
        or value.get("preparation_version") != PREPARATION_VERSION
    ):
        raise ValueError(
            "unsupported preparation/format version; rebuild legacy or incompatible indexes in a new directory"
        )
