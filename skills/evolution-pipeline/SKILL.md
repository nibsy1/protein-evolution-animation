---
name: evolution-pipeline
description: One-command pipeline chaining sequence-scraper, structure-aligner,
  and morph-animation: from a protein query (gene name, AGI, or protein name) to an
  evolutionary morph animation GIF with all intermediate artifacts. v2 chain uses
  the ribbon-style, stability-masked morph renderer. Validated end-to-end on ECA3
  and rbcL.
version: 2
---
# Evolution Animation Pipeline (all-in-one)

One command from protein query to evolutionary morph animation.
Chains the three workspace skills: sequence-scraper → structure-aligner
→ morph-animation.

## Usage

```bash
python scripts/run_evolution_pipeline.py ECA3 --email you@example.com
python scripts/run_evolution_pipeline.py AT1G10130
python scripts/run_evolution_pipeline.py "ribulose bisphosphate carboxylase large chain"
```

## What it does

1. **sequence-scraper** — resolves the query (gene name / AGI / protein
   name) and scrapes one verified homolog per model species (16 species,
   green algae → angiosperms) → `sequences.fasta` + `sequences.tsv`
2. **structure-aligner** (v2) — fetches AlphaFold DB models (sequence-
   verified) or threads missing ones, superposes all into one 3D frame,
   writes ribbon-ready PDBs + `aligned.json` with seq/ss → `structures/`
3. **morph-animation** (v2) — renders the ribbon-style, stability-masked
   morphing GIF in evolutionary order → `evolution.gif` + `evolution_sheet.png`

Everything lands in `./<QUERY>/`.

## Requirements

The three skill scripts must be findable: put them (and this script)
in one directory — e.g. the skill assets folder `scripts/` — or in the
working directory. Dependencies: Python 3.8+, requests, numpy,
biopython, matplotlib, pillow, plus `src/proteinlib.py` from the workspace.

## Validated

- **ECA3** (16 species, 500-frame GIF): ~8 min total (v1 chain)
- **rbcL** (16 species, 500-frame GIF): 188 s end-to-end (v1 chain)
- v2 renderer validated on the ECA3 algae↔Arabidopsis pair (ribbon style,
  fixed medium lighting, 19% of core residues frozen during the morph)

## Notes

- Pass `--hold`/`--morph` to tune animation pacing; the morph skill also
  accepts `--disp-max` to control how much change is needed before a
  residue animates.
- The wrapper aborts with a clear message if any stage fails or
  produces no output.
- For a fully automated agent (scheduled runs, retries, delivery),
  wrap this script in a Mistral workflow next.
