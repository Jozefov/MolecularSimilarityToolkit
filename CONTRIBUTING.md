# Contributing

Use Python 3.11–3.13. Install the locked development environment with
`uv sync --locked --group dev`, then run `uv run ruff check .`,
`uv run pytest`, and `uv run mkdocs build --strict`.

Keep changes focused. New search or filtering behavior needs an independent
RDKit oracle test. Changes to storage or scheduling need an interruption and
resume test. Scientific results must agree across worker counts and shard sizes.

Do not commit molecular databases, benchmark inputs, credentials, or local
indexes. Small hand-written test molecules are welcome. State the source and
license of any third-party fixture. Benchmark reports must identify inputs,
chemistry, hardware, memory, and cache conditions. Never describe synthetic
fingerprints as real compounds.

Please discuss storage-format changes before implementing them. Existing
indexes are immutable; unsupported formats must fail clearly.

The repository starts private. Public visibility, hosted documentation, and
package publication require a separate release review by the owner.
