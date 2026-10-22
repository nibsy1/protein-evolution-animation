---
name: morph-animation
description: Renders an evolutionary morph animation GIF from the aligned.json
  produced by the structure-aligner skill v2. v2 renders a clean ribbon-style
  cartoon (helices as wide ribbons, strands as arrows, loops as a thin tube) with
  fixed medium lighting, and morphs with a stability mask so only genuinely
  changing regions move. No membrane visualization. Steps through species in
  evolutionary order (green algae to angiosperms). Outputs GIF + frame sheet PNG.
  v3 anchors residue correspondence through the alignment reference and labels
  the N-/C-termini.
version: 3
---
# Morph Animation (v3)

Renders the evolutionary-development animation: structures morph from
green algae to flowering plants in a smooth rotating 3D GIF.

## Usage

```bash
python scripts/morph_animation.py structures/aligned.json -o evolution.gif --hold 20 --morph 12
```

## What's new in v2

- **Ribbon style** — helices render as wide flat ribbons (linewidth 8.5) with
  darker edge passes on both sides, strands as wide ribbons (linewidth 6.5)
  with arrow heads, loops as a thin tube (linewidth 2.2). SS-class colors
  (purple helix / yellow strand / grey loop) are tinted per segment by mean pLDDT,
  so confidence is still readable.
- **Fixed medium lighting** — a constant view-space light direction with
  0.45 ambient + 0.55 diffuse lambert shading. Shading is identical on every
  frame and while the camera rotates; the structure never appears to rotate
  under a static light.
- **Stability-masked morphing** (the big one). For each pair of consecutive
  structures, every corresponding residue gets a morph weight from its Cα
  displacement and whether its secondary-structure class changed:
  - stable + same SSE (a helix that is a helix in both species) → **frozen**,
    it holds its position for the whole transition;
  - large displacement or SSE change (new/lost helix or strand, moved loop) →
    **animated**, travelling the full A→B distance with smoothstep easing;
  - insertions (residues present only in the target) **grow in** from the
    source centroid instead of popping in.
  This removes the old "whole protein reshuffles with an abrupt end" behaviour:
  only the biggest changes are visualized.
- **No membrane planes** anywhere.
- **Reference-anchored correspondence (v3)** — residue pairs for the morph are
  computed by aligning both structures to the alignment reference and
  intersecting (`pl.reference_pairs`), not by a direct pairwise alignment.
  A pairwise alignment between two distant homologs can slip over long
  terminal extensions, which made the whole protein appear to reshuffle;
  the reference-anchored path keeps every species on the same residue
  numbering, so conserved regions stay frozen.
- **N-/C-terminus labels (v3)** — every frame and sheet panel marks the
  N-terminus (blue "N") and C-terminus (orange "C") with a small sphere +
  text label, so chain direction is always visible.
- The morph math lives in `src/proteinlib.py` (`StabilityMask`, `morph_frames`,
  `reference_pairs`, `chain_orientation_ok`) and is shared with the interactive
  viewer and the library agents, so GIF and web viewer transitions behave
  identically.

## How it works

1. **Ordering** — structures sorted by lineage (Green Algae →
   Bryophytes → Ferns → Gymnosperms → Angiosperms); within a lineage
   by core identity (highest-confidence first).
2. **Stability mask** — per transition, local BLOSUM62 correspondences +
   per-residue displacement + SS-class comparison → morph weights.
3. **Masked morph** — stable residues frozen, changing residues eased
   (smoothstep) from A to B, insertions fade in. pLDDT colors cross-fade.
4. **Rendering** — each structure holds for N frames (default 25), then
   morphs over M frames (default 15), with slow continuous rotation, ribbon
   cartoon, fixed medium lighting, dark theme, title = lineage + species,
   subtitle = accession + source + core identity, and an
   evolutionary-progression percentage bar.

## Outputs

- `<name>.gif` — the animation (default evolution.gif)
- `<name>_sheet.png` — one still per species for quick inspection

## Tuning

- `--disp-max` (default 4.0 Å): Cα displacement at which a residue animates
  fully. Lower = more residues frozen (only big changes move).
- `--hold` / `--morph` / `--fps` / `--size`: pacing and resolution.

## Validation (ECA3 pair, algae ↔ Arabidopsis)

- 992 corresponding residues, core identity 0.62, core-refined fit 2.47 Å
- 188 residues (19%) frozen: mean displacement 1.26 Å (max 1.73 Å) — these hold
  position through the whole morph; animated residues average 3.63 Å of travel
- 49 secondary-structure class changes, all fully animated
- Rendered frames verified: ribbon cartoon, consistent lighting, no membrane

## Dependencies

Python 3.8+, `numpy`, `biopython`, `matplotlib`, `pillow`, plus
`src/proteinlib.py` from the workspace.
