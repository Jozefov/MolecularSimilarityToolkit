# Artwork sources

The SVG files are the editable originals; PNG files are exports for tools that
do not support SVG. Both are covered by the repository's MIT license.

`logo.svg` and `logo-dark.svg` combine molecular nodes with a ranked stack.
`workflow*.svg` use RDKit molecule drawings and measured Morgan radius-2,
2,048-bit Tanimoto examples. The small bit grids are schematic, not literal
fingerprints. A separate tall composition is selected on narrow screens.

Rebuild workflow SVGs with `uv run python scripts/generate_assets.py`. Export
PNGs with Inkscape, for example:

```bash
inkscape docs/assets/workflow.svg --export-type=png --export-filename=docs/assets/workflow.png
```

README visual review covers 1,280-pixel desktop and 390-pixel mobile viewports,
in both color schemes. The SVG title/description and Markdown alt text describe
the information without relying on color alone.
