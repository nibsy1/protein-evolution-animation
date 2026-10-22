#!/usr/bin/env python3
"""One-command protein evolution animation pipeline.

Chains the three skills:
    1. sequence-scraper        -> sequences FASTA + TSV
    2. structure-aligner       -> aligned structures + aligned.json
    3. morph-animation         -> evolution GIF + frame sheet

Usage:
    python run_evolution_pipeline.py ECA3
    python run_evolution_pipeline.py AT1G10130 --email you@example.com
    python run_evolution_pipeline.py "ribulose bisphosphate carboxylase large chain"

Outputs are written to ./<QUERY>/ :
    <QUERY>/sequences.fasta     scraped homologs (one per species)
    <QUERY>/sequences.tsv       quality table
    <QUERY>/structures/        aligned PDBs, aligned.json, overlay.png
    <QUERY>/evolution.gif      the final animation
    <QUERY>/evolution_sheet.png one still per species
"""
import argparse
import os
import re
import subprocess
import sys
import time

SCRIPTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)))


def find_script(name, search_paths):
    for base in search_paths:
        for root, _, files in os.walk(base):
            if name in files:
                return os.path.join(root, name)
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("query",
                   help="gene name (ECA3), AGI (AT1G10130), or protein name")
    p.add_argument("-o", "--outdir", default=None,
                   help="output directory (default: ./<QUERY>)")
    p.add_argument("--email", default=None,
                   help="contact email for NCBI (recommended)")
    p.add_argument("--hold", type=int, default=20,
                   help="animation hold frames per structure")
    p.add_argument("--morph", type=int, default=12,
                   help="animation morph frames per transition")
    args = p.parse_args()

    safe = re.sub(r"[^A-Za-z0-9_-]+", "_", args.query)[:40] or "query"
    outdir = args.outdir or safe
    os.makedirs(outdir, exist_ok=True)
    struct_dir = os.path.join(outdir, "structures")
    fasta = os.path.join(outdir, "sequences.fasta")
    table = os.path.join(outdir, "sequences.tsv")
    gif = os.path.join(outdir, "evolution.gif")

    search = [SCRIPTS_DIR,
              os.path.join(SCRIPTS_DIR, "scripts"),
              os.path.join(SCRIPTS_DIR, "..", "scripts"),
              os.getcwd()]
    scraper = find_script("scrape_sequences.py", search)
    aligner = find_script("fetch_and_align_structures.py", search)
    animator = find_script("morph_animation.py", search)
    missing = [n for n, s in
               [("sequence-scraper", scraper),
                ("structure-aligner", aligner),
                ("morph-animation", animator)] if s is None]
    if missing:
        sys.exit(f"Could not find scripts for skill(s): {', '.join(missing)}. "
                 "Download the three skill scripts into ./scripts/ first.")

    def run(step, cmd):
        print(f"\n{'=' * 60}\n[{step}] {' '.join(cmd)}\n{'=' * 60}",
              flush=True)
        t0 = time.time()
        r = subprocess.run(cmd)
        if r.returncode != 0:
            sys.exit(f"[{step}] failed with exit code {r.returncode}")
        print(f"[{step}] done in {time.time() - t0:.0f}s", flush=True)

    email_args = ["--email", args.email] if args.email else []

    run("1/3 sequence-scraper",
        [sys.executable, scraper, args.query, "-o", fasta,
         "--table", table] + email_args)

    if not os.path.exists(fasta) or os.path.getsize(fasta) < 100:
        sys.exit("Scraper produced no sequences; aborting.")

    run("2/3 structure-aligner",
        [sys.executable, aligner, fasta, "-o", struct_dir] + email_args)

    aligned_json = os.path.join(struct_dir, "aligned.json")
    if not os.path.exists(aligned_json):
        sys.exit("Aligner produced no aligned.json; aborting.")

    run("3/3 morph-animation",
        [sys.executable, animator, aligned_json, "-o", gif,
         "--hold", str(args.hold), "--morph", str(args.morph)])

    print(f"\n{'=' * 60}\nPipeline complete.\n"
          f"  animation : {gif}\n"
          f"  structures: {struct_dir}\n"
          f"  sequences : {fasta}\n{'=' * 60}")


if __name__ == "__main__":
    main()
