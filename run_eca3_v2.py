#!/usr/bin/env python3
"""Re-run the ECA3 pipeline with the v3 aligner + v3 animator + v10 viewer.

Uses the cached sequences.fasta and the cached AlphaFold PDBs (pre-seeded
into the structure fetcher's raw cache) so only UniProt resolution and
taxonomy lookups hit the network.
"""
import asyncio, json, os, shutil, sys, time

sys.path.insert(0, "/workspace/workflow")
sys.path.insert(0, "/workspace/workflow/src")

from protein_animation.activities import (align_structures, build_viewer,
                                          render_animation)

OUT = "/workspace/outputs/ECA3_run_v2"
os.makedirs(OUT, exist_ok=True)

EMAIL = "nils.woerz@student.uni-halle.de"

async def main():
    t0 = time.time()
    fasta = open(f"{OUT}/sequences.fasta").read()

    print("[1/3] align_structures (v3 aligner, threaded seq fix) ...")
    t1 = time.time()
    aln = await align_structures("ECA3", fasta, email=EMAIL)
    print(f"  done in {time.time()-t1:.0f}s: {aln['n_structures']} structures, "
          f"ref {aln['reference']}")
    aj_path = f"{OUT}/aligned.json"
    with open(aj_path, "w") as f:
        f.write(aln["aligned_json"])
    with open(f"{OUT}/summary.tsv", "w") as f:
        f.write(aln["summary"])
    os.makedirs(f"{OUT}/structures/aligned", exist_ok=True)
    for name, b64 in aln["pdbs"].items():
        import base64
        with open(f"{OUT}/structures/aligned/{name}", "wb") as f:
            f.write(base64.b64decode(b64))
    if aln.get("overlay_png"):
        import base64
        with open(f"{OUT}/overlay.png", "wb") as f:
            f.write(base64.b64decode(aln["overlay_png"]))

    print("[2/3] render_animation (v3 animator, ref-anchored + N/C labels) ...")
    t1 = time.time()
    anim = await render_animation("ECA3", aln["aligned_json"],
                                  hold=20, morph=12, fps=12, size=(720, 720))
    import base64
    with open(f"{OUT}/evolution.gif", "wb") as f:
        f.write(base64.b64decode(anim["gif"]))
    if anim.get("sheet_png"):
        with open(f"{OUT}/evolution_sheet.png", "wb") as f:
            f.write(base64.b64decode(anim["sheet_png"]))
    print(f"  done in {time.time()-t1:.0f}s: {anim['n_frames']} frames")

    print("[3/3] build_viewer (v10 viewer, smooth morph + N/C markers) ...")
    t1 = time.time()
    v = await build_viewer("ECA3", aln["aligned_json"],
                           title="ECA3 protein structure across lineages")
    with open(f"{OUT}/viewer.html", "wb") as f:
        f.write(base64.b64decode(v["viewer_html"]))
    print(f"  done in {time.time()-t1:.0f}s: viewer.html")

    print(f"\nPIPELINE COMPLETE in {time.time()-t0:.0f}s -> {OUT}")

asyncio.run(main())
