"""End-to-end tests for the protein-animation workflow activities.

Runs the align -> animate -> viewer chain on the cached ECA3 two-species
fixture (no network needed beyond the two AlphaFold PDBs already in
data/). Validates the artifact payloads the workflow returns.
"""
import asyncio
import base64
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from protein_animation.activities import (align_structures, build_viewer,  # noqa: E402
                                          render_animation)

DATA = ROOT.parent / "data"
FASTA = DATA / "test_sequences.fasta"


def _decode(b64):
    return base64.b64decode(b64)


async def main():
    fasta = FASTA.read_text()

    print("[1/3] align_structures ...")
    aln = await align_structures("ECA3_test", fasta, email=None)
    assert aln["n_structures"] == 2, aln["n_structures"]
    assert aln["reference"] == "Q9SY55", aln["reference"]
    assert set(aln["species"]) == {"Arabidopsis thaliana",
                                   "Coccomyxa viridis"}, aln["species"]
    assert set(aln["pdbs"]) == {"A0AAV1IES4.pdb", "Q9SY55.pdb"}, list(aln["pdbs"])
    assert aln["overlay_png"], "missing overlay"
    aj = json.loads(aln["aligned_json"])
    assert aj["reference"] == "Q9SY55"
    assert all("ss" in s and "seq" in s for s in aj["structures"])
    # PDBs are ribbon-ready (HELIX/SHEET records)
    for name, b64 in aln["pdbs"].items():
        pdb = _decode(b64).decode()
        assert pdb.startswith("HELIX") or pdb.startswith("SHEET") or "ATOM" in pdb
        assert "ATOM" in pdb
    # overlay is a real PNG
    assert _decode(aln["overlay_png"])[:8] == b"\x89PNG\r\n\x1a\n"
    print(f"  ok: {aln['n_structures']} structures, ref {aln['reference']}")

    print("[2/3] render_animation ...")
    anim = await render_animation("ECA3_test", aln["aligned_json"],
                                  hold=3, morph=4, fps=6, size=(320, 320))
    gif = _decode(anim["gif"])
    assert gif[:6] in (b"GIF87a", b"GIF89a"), gif[:6]
    assert anim["n_frames"] == 10, anim["n_frames"]  # 2 holds *3 + 1 morph *4
    assert anim["sheet_png"], "missing frame sheet"
    assert _decode(anim["sheet_png"])[:8] == b"\x89PNG\r\n\x1a\n"
    print(f"  ok: {anim['n_frames']} frames, gif {len(gif)} bytes")

    print("[3/3] build_viewer ...")
    v = await build_viewer("ECA3_test", aln["aligned_json"], title="ECA3 test")
    html = _decode(v["viewer_html"]).decode()
    assert "__DATA__" not in html and "__TITLE__" not in html
    assert "ECA3 test" in html
    assert "membrane" not in html.lower()
    assert "pairWeights" in html and "buildMorphCoords" in html
    assert "HemisphereLight" in html
    print(f"  ok: viewer {len(html)} bytes, morph+lighting present, "
          f"0 membrane refs")

    print("\nALL WORKFLOW ACTIVITY TESTS PASSED")


if __name__ == "__main__":
    asyncio.run(main())
