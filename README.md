# Plant Protein Evolution → Structure Animation

Animates the evolutionary development of a plant protein structure
(default target: **ECA3**, the ER-type Ca²⁺-transporting ATPase) across
green algae → mosses → ferns → conifers → flowering plants.

Give it a gene name (e.g. `ECA3`) or an Arabidopsis AGI locus and it fetches
the homologous sequences, predicts/fetches their 3D structures, superposes
them into one common frame, and renders a smooth ribbon-style morph animation
(GIF) plus an interactive HTML viewer — showing only the *real* structural
changes between species, with conserved helices/strands held frozen.

## What's in this repo

```
skills/                 the four Mistral skills (SKILL.md + scripts)
  sequence-scraper/     fetch homologous protein sequences (NCBI/UniProt)
  structure-aligner/    fetch 3D structures (AlphaFold DB + threading fallback),
                        superpose all into a common frame            (v3)
  morph-animation/      render the evolutionary morph GIF, ribbon style,
                        stability-masked smooth transitions, N/C labels (v3)
  interactive-viewer/   single-file interactive HTML viewer (Three.js),
                        same smooth morph + N/C markers as the GIF    (v10)
  evolution-pipeline/   end-to-end runner tying the skills together
src/                    shared library + multi-agent pipeline
  proteinlib.py         SS assignment, BLOSUM62 correspondence, Kabsch
                        superposition, StabilityMask, morph_frames,
                        reference-anchored correspondence, orientation check
  agents/               taxonomy / search / structure / animation agents
  state/, config.py, main.py, run_pipeline.py
workflow/               Mistral "protein animation" workflow (activities + workflow def)
workflow_for_user/      user-facing copy of the workflow (+ DEPLOY.md)
tests/                  unit tests for proteinlib
outputs/ECA3_run_v2/    the canonical ECA3 demo run (GIF, viewer, PDBs, ...)
workflow_demo/          small 2-species end-to-end demo outputs
data/                   test sequences + test structures
```

## Quick start

```bash
# setup
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python -r requirements.txt

# run the full pipeline for a gene (needs network: NCBI, UniProt, AlphaFold DB)
.venv/bin/python skills/evolution-pipeline/scripts/run_evolution_pipeline.py ECA3 --email you@example.com

# or step by step
.venv/bin/python skills/sequence-scraper/scripts/scrape_sequences.py ECA3 -o ECA3.fasta --table ECA3.tsv
.venv/bin/python skills/structure-aligner/scripts/fetch_and_align_structures.py ECA3.fasta -o structures/ --email you@example.com
.venv/bin/python skills/morph-animation/scripts/morph_animation.py structures/aligned.json -o evolution.gif --hold 20 --morph 12
.venv/bin/python skills/interactive-viewer/scripts/interactive_viewer.py structures/aligned.json -o viewer.html

# tests
.venv/bin/python tests/test_proteinlib.py
```

## The smooth, change-only transition (the core idea)

For each pair of consecutive structures (already superposed in a common
frame by the aligner):

1. **Reference-anchored correspondence** — residue pairs are computed by
   aligning *both* structures to the alignment reference and intersecting
   (`pl.reference_pairs`), not by a direct pairwise alignment. A direct
   pairwise alignment between distant homologs can slip over long terminal
   extensions, which made the whole protein appear to reshuffle.
2. **Stability mask** — per corresponding residue: Cα displacement + whether
   the H/E/C secondary-structure class changed → a morph weight in [0, 1]
   (smoothstep-ed).
3. **Weight ≈ 0 → frozen**: a helix that is a helix in both species holds
   its position for the whole transition (settles onto the target only at
   the very end, a sub-Å correction).
4. **Weight = 1 → animated**: genuinely changed regions travel the full
   source→target distance with smoothstep easing.
5. **Insertions** (residues present only in the target) grow in from the
   source centroid instead of popping in.

Measured on the ECA3 run (16 species, 15 transitions): **16–94% of
corresponding residues stay frozen per transition** (median ≈ 64%). The
transitions that used to reshuffle the whole protein (0.0–0.3% frozen in
the old version) now freeze 49–72%.

Tune with `--disp-max` (default 4.0 Å): lower values freeze more residues,
so only the biggest changes move.

## Ribbon style & lighting

- **GIF**: helices as wide flat ribbons (linewidth 8.5) with darker edge
  passes, strands as wide ribbons (6.5) with arrow heads, loops as a thin
  tube (2.2). Fixed view-space light (0.45 ambient + 0.55 diffuse lambert)
  → shading is identical on every frame and while the camera rotates.
- **Viewer**: Three.js ribbon boxes (helices), flat arrows (strands), smooth
  tube (loops); camera-attached hemisphere + directional light at constant
  medium intensity.
- **No membrane visualization** anywhere.
- **N-/C-termini labeled** on every GIF frame and in the viewer (blue "N"
  sphere, orange "C" sphere) so chain direction is always visible.
- Aligned PDBs carry HELIX/SHEET records, so they open as ribbons in
  PyMOL / ChimeraX / Mol* too.

## Orientation check

`pl.chain_orientation_ok` verifies every structure's chain runs in the same
direction as the reference (superpose via the sequence correspondence, then
check which reference terminus each chain end is near). All 16 ECA3
structures pass — no reversed N-/C-termini.

## Demo output

`outputs/ECA3_run_v2/` contains the full ECA3 run:

| File | What it is |
|---|---|
| `evolution.gif` | 500-frame morph animation across 16 plant species |
| `viewer.html` | interactive viewer with the same smooth morph |
| `evolution_sheet.png` | one still per species |
| `aligned.json` | all structures, superposed, with seq/SS/coords |
| `sequences.fasta` / `summary.tsv` | input sequences / per-species summary |
| `structures/aligned/*.pdb` | 16 ribbon-ready aligned PDBs |
| `manifest.json` | file sizes + sha256 checksums |

## Developing with Mistral

Your skills also live in the Mistral skills registry and your outputs in the
document library — both reachable from a Mistral agent session via the
`mistral_mcp_*` and `document_library_*` tools. This repo is the canonical,
version-controlled copy. Workflow:

1. Edit and test here with the venv above.
2. When a skill is ready, publish it back with `mistral_mcp_write_skill`.
3. Upload new outputs (GIF/HTML) to the library with `mistral_mcp_create_document`.

## License

MIT
