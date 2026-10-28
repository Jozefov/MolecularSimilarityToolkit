import importlib.util
from pathlib import Path

from readme_renderer.markdown import render
from typer.testing import CliRunner

import molecular_similarity_toolkit as mst
from molecular_similarity_toolkit.cli import app


def test_package_readme_is_generated_and_renders_without_relative_assets():
    root = Path(__file__).parents[1]
    spec = importlib.util.spec_from_file_location(
        "package_readme", root / "scripts/build_package_readme.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    text = (root / "PACKAGE_README.md").read_text()
    assert text == module.build((root / "README.md").read_text())
    html = render(text)
    assert html and "<img" in html and "<source" not in html
    assert 'href="docs/' not in html and 'src="docs/' not in html


def test_cli_export_and_compact(tmp_path):
    pool = mst.build_index(
        mst.Source.from_records([{"smiles": "CCO"}]), tmp_path / "pool"
    )
    runner = CliRunner()
    out = tmp_path / "compact"
    result = runner.invoke(app, ["compact", str(pool.path), str(out)])
    assert result.exit_code == 0, result.output
    hits = mst.search(
        out, mst.Source.from_records([{"smiles": "CCO"}]), 1, tmp_path / "search"
    )
    result = runner.invoke(app, ["export", str(hits.path), str(tmp_path / "hits.tsv")])
    assert result.exit_code == 0, result.output
    assert len((tmp_path / "hits.tsv").read_text().splitlines()) == 2
