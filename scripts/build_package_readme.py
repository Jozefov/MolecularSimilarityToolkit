"""Create a PyPI description with absolute links and static PNG illustrations."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = "https://github.com/Jozefov/MolecularSimilarityToolkit"
RAW = "https://raw.githubusercontent.com/Jozefov/MolecularSimilarityToolkit/main"


def build(text):
    def picture(match):
        image = re.search(
            r'<img[^>]*alt="([^"]*)"[^>]*src="([^"]*)"', match.group(0), re.DOTALL
        )
        alt, path = image.groups()
        return f"![{alt}]({RAW}/{path.removesuffix('.svg')}.png)"

    text = re.sub(r"<picture>.*?</picture>", picture, text, flags=re.DOTALL)
    text = re.sub(
        r"\]\(((?:docs/|CONTRIBUTING|LICENSE)[^)]*)\)",
        lambda m: f"]({REPO}/blob/main/{m.group(1)})",
        text,
    )
    return text


if __name__ == "__main__":
    (ROOT / "PACKAGE_README.md").write_text(build((ROOT / "README.md").read_text()))
