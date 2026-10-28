"""Rebuild original editable vector artwork (PNG export uses Inkscape)."""

from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D

OUT = Path(__file__).resolve().parents[1] / "docs" / "assets"
MOLS = {
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "methyl": "CC(=O)Oc1ccccc1C(=O)OC",
    "ethyl": "CC(=O)Oc1ccccc1C(=O)OCC",
    "salicylic": "O=C(O)c1ccccc1O",
    "ethanol": "CCO",
    "benzene": "c1ccccc1",
}
GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
FPS = {k: GEN.GetFingerprint(Chem.MolFromSmiles(v)) for k, v in MOLS.items()}


def molecule(key, x, y, w=180, h=120, dark=False):
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    opt = drawer.drawOptions()
    opt.clearBackground = False
    opt.bondLineWidth = 2.1
    opt.minFontSize = 14
    ink = (0.86, 0.93, 0.98) if dark else (0.12, 0.23, 0.34)
    opt.updateAtomPalette(
        {6: ink, 0: ink, 8: (0.88, 0.38, 0.17), 7: (0.18, 0.52, 0.85)}
    )
    drawer.DrawMolecule(Chem.MolFromSmiles(MOLS[key]))
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    inner = svg[
        svg.index("<!-- END OF HEADER -->")
        + len("<!-- END OF HEADER -->") : svg.rindex("</svg>")
    ]
    return f'<g transform="translate({x},{y})">{inner}</g>'


def make(dark=False):
    bg = "#091827" if dark else "#f4f8fb"
    paper = "#102b40" if dark else "#ffffff"
    ink = "#e8f2fa" if dark else "#123047"
    muted = "#9bb7cb" if dark else "#557183"
    line = "#28485e" if dark else "#d4e2ea"
    teal = "#42d6b2" if dark else "#04816d"
    blue = "#63adff" if dark else "#2362bd"
    amber = "#ffbf64" if dark else "#b56b0b"

    def text(x, y, s, size=25, color=ink, weight=400):
        return f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" font-weight="{weight}">{s}</text>'

    def rect(x, y, w, h, stroke=line, fill=paper, r=18):
        return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="1.5"/>'

    def tag(x, y, label, color):
        return f'<circle cx="{x}" cy="{y - 7}" r="6" fill="{color}"/>' + text(
            x + 17, y, label, 20, muted
        )

    def bits(x, y, cols=10, rows=4, seed=0, size=8, color=teal):
        return "".join(
            f'<rect x="{x + i * (size + 5)}" y="{y + j * (size + 5)}" width="{size}" height="{size}" rx="2" fill="{color if (i * 7 + j * 11 + seed) % 5 < 2 else line}"/>'
            for j in range(rows)
            for i in range(cols)
        )

    pieces = [
        f'''<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="1040" viewBox="0 0 1600 1040" role="img" aria-labelledby="title desc">
<title id="title">Molecules in. Close matches or a protected pool out.</title><desc id="desc">Six example molecules enter a packed fingerprint index. Aspirin finds ranked matches with Tanimoto scores 1.00, 0.67 and 0.61. With aspirin protected and cutoff greater than 0.8, aspirin is removed while less similar structures are kept. Scores use Morgan radius 2, 2048 bits. Molecular drawings are examples, and the bit grids are schematic.</desc>
<defs><linearGradient id="flow" x1="0" x2="1"><stop stop-color="{teal}"/><stop offset="1" stop-color="{blue}"/></linearGradient><pattern id="dots" width="24" height="24" patternUnits="userSpaceOnUse"><circle cx="2" cy="2" r="1" fill="{line}"/></pattern><marker id="arr" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0 0 10 5 0 10Z" fill="{blue}"/></marker><filter id="shadow" x="-10%" y="-10%" width="130%" height="135%"><feDropShadow dx="0" dy="6" stdDeviation="9" flood-color="#173b55" flood-opacity=".08"/></filter></defs>
<style>text{{font-family:Arial,Helvetica,sans-serif}}.card{{filter:url(#shadow)}}</style>
<rect width="1600" height="1040" rx="26" fill="{bg}"/><rect x="20" y="115" width="1560" height="870" fill="url(#dots)" opacity=".55"/>'''
    ]
    pieces += [
        text(50, 65, "Molecules in. Useful choices out.", 40, ink, 700),
        text(50, 104, "MolecularSimilarityToolkit", 21, muted),
    ]
    # Collection samples: overlapping cards give the collection physical depth.
    pieces += [text(50, 183, "YOUR COLLECTIONS", 20, muted, 700)]
    for i, (x, y, key) in enumerate(
        [
            (50, 216, "aspirin"),
            (269, 240, "methyl"),
            (65, 421, "ethyl"),
            (284, 445, "salicylic"),
            (50, 626, "ethanol"),
            (269, 650, "benzene"),
        ]
    ):
        pieces += [
            rect(x + 7, y - 7, 204, 174),
            "<g>",
            rect(x, y, 204, 174),
            molecule(key, x + 12, y + 5, 180, 118, dark),
            bits(x + 17, y + 139, 12, 1, i, 8),
            "</g>",
        ]
    pieces += [
        tag(67, 866, "PubChem", teal),
        tag(239, 866, "ZINC", blue),
        text(67, 900, "…or any molecule files", 21, muted),
    ]
    # Molecules condense into packed bits and a reusable shard stack.
    pieces += [
        '<path d="M494 514C548 514 539 514 580 514" fill="none" stroke="url(#flow)" stroke-width="5" marker-end="url(#arr)"/>'
    ]
    pieces += [text(583, 282, "PACK ONCE", 20, muted, 700)]
    for i in (2, 1, 0):
        pieces += [rect(587 + i * 12, 334 + i * 21, 260, 288, teal)]
    pieces += [
        text(610, 376, "Fingerprint index", 25, ink, 700),
        bits(612, 402, 13, 8, 2, 10),
        text(611, 584, "2,048 bits / molecule", 21, muted),
    ]
    pieces += [tag(609, 716, "Reuse", teal), tag(733, 716, "Resume", blue)]
    # Curved search and filtering paths.
    pieces += [
        f'<path d="M874 480C940 480 906 270 984 270" fill="none" stroke="{blue}" stroke-width="4" marker-end="url(#arr)"/>',
        f'<path d="M874 548C952 548 911 740 984 740" fill="none" stroke="{blue}" stroke-width="4" marker-end="url(#arr)"/>',
    ]
    pieces += [text(999, 183, "01  CLOSEST CANDIDATES", 20, blue, 700)]
    pieces += [
        rect(998, 206, 548, 122),
        molecule("aspirin", 1007, 207, 166, 113, dark),
        text(1190, 254, "Query: aspirin", 25, ink, 700),
        text(1190, 293, "Closest first ↓", 22, muted),
    ]
    for i, key in enumerate(["aspirin", "methyl", "ethyl"]):
        x = 998 + i * 185
        score = DataStructs.TanimotoSimilarity(FPS["aspirin"], FPS[key])
        pieces += [
            "<g>",
            rect(x, 347, 177, 192),
            text(x + 14, 379, str(i + 1), 23, blue, 700),
            molecule(key, x + 3, 383, 169, 105, dark),
            text(x + 18, 518, f"{score:.2f}", 25, blue, 700),
            f'<rect x="{x + 81}" y="503" width="{78 * score}" height="12" rx="6" fill="{blue}"/>',
            "</g>",
        ]
    pieces += [text(998, 613, "02  PROTECTED TRAINING POOL", 20, teal, 700)]
    pieces += [
        rect(998, 638, 548, 119, amber),
        molecule("aspirin", 1007, 640, 166, 113, dark),
        text(1188, 680, "Protected: aspirin", 25, ink, 700),
        text(1188, 720, "Remove matches > 0.8", 22, amber),
    ]
    pieces += [
        rect(998, 778, 177, 191, amber),
        molecule("aspirin", 1001, 791, 169, 110, dark),
        f'<path d="m1015 807 22 22m0-22-22 22" stroke="{amber}" stroke-width="4"/>',
        text(1016, 942, "Remove  1.00", 22, amber, 700),
    ]
    pieces += [
        rect(1190, 778, 356, 191, teal),
        molecule("methyl", 1198, 791, 169, 110, dark),
        molecule("ethyl", 1367, 791, 169, 110, dark),
        f'<path d="m1205 819 8 8 18-23" stroke="{teal}" stroke-width="4" fill="none"/>',
        text(1209, 942, "Keep  0.67", 22, teal, 700),
        text(1384, 942, "Keep  0.61", 22, teal, 700),
    ]
    pieces += [
        text(
            50,
            1009,
            "Real example scores · Morgan radius 2 / 2,048 bits · Bit grids are schematic",
            18,
            muted,
        ),
        "</svg>",
    ]
    return "".join(pieces)


for dark in (False, True):
    (OUT / ("workflow" + ("-dark" if dark else "") + ".svg")).write_text(make(dark))

# A taller composition keeps molecule examples readable on a phone.
for dark in (False, True):
    bg, ink, muted, paper, teal, blue = (
        ("#091827", "#e8f2fa", "#9bb7cb", "#102b40", "#42d6b2", "#63adff")
        if dark
        else ("#f4f8fb", "#123047", "#557183", "#ffffff", "#04816d", "#2362bd")
    )
    desktop = make(dark)
    inner = desktop[desktop.index("<defs>") : desktop.rindex("</svg>")]
    mobile = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="700" height="1510" viewBox="0 0 700 1510" role="img" aria-labelledby="mobile-title"><title id="mobile-title">Real molecule examples: closest matches and a protected pool</title><rect width="700" height="1510" rx="24" fill="{bg}"/><g font-family="Arial,Helvetica,sans-serif" fill="{ink}"><text x="30" y="60" font-size="34" font-weight="700">Molecules in. Useful choices out.</text><text x="30" y="107" font-size="24" fill="{muted}">PubChem, ZINC, or your own files</text></g>'
    ]
    for i, key in enumerate(["aspirin", "methyl", "ethyl", "salicylic"]):
        x = 28 + i * 164
        mobile += [
            f'<rect x="{x}" y="143" width="151" height="151" rx="16" fill="{paper}" stroke="{muted}" stroke-opacity=".3"/>',
            molecule(key, x + 4, 158, 143, 122, dark),
        ]
    mobile += [
        f'<path d="M350 309v30m-10-10 10 10 10-10" stroke="{blue}" stroke-width="4" fill="none"/><rect x="70" y="358" width="560" height="146" rx="18" fill="{paper}" stroke="{teal}" stroke-width="2"/><text x="101" y="404" font-family="Arial,Helvetica,sans-serif" font-size="29" font-weight="700" fill="{ink}">Reusable fingerprint index</text>'
    ]
    for j in range(2):
        for i in range(29):
            mobile += [
                f'<rect x="{104 + i * 17}" y="{428 + j * 23}" width="12" height="12" rx="2" fill="{teal if (i * 7 + j * 3) % 5 < 2 else muted}" opacity="{1 if (i * 7 + j * 3) % 5 < 2 else 0.25}"/>'
            ]
    mobile += [
        f'<svg x="30" y="540" width="640" height="434" viewBox="990 165 560 380" overflow="hidden">{inner}</svg>',
        f'<svg x="30" y="1002" width="640" height="434" viewBox="990 600 560 380" overflow="hidden">{inner}</svg>',
        f'<text x="30" y="1480" font-family="Arial,Helvetica,sans-serif" font-size="20" fill="{muted}">Morgan radius 2 / 2,048 bits · Actual example scores</text></svg>',
    ]
    (OUT / ("workflow-mobile" + ("-dark" if dark else "") + ".svg")).write_text(
        "".join(mobile)
    )
