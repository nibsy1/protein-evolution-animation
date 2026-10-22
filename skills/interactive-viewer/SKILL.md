---
name: interactive-viewer
description: Generates an interactive HTML viewer: clickable phylogenetic tree +
  rotatable 3D ribbon-style structure (Three.js, works from CA-only traces) with
  consistent medium camera-attached lighting, hydropathy-tinted SS coloring, and a
  stability-masked smooth morph on species selection. No membrane visualization.
  Single self-contained HTML from aligned.json (structure-aligner v2).
  v10: smooth ~1s eased morph with color/SS cross-fade (matches the GIF),
  N-/C-terminus markers.
version: 10
---
# Interactive Evolution Viewer (v10)

Generates a single self-contained HTML file: an interactive
phylogenetic tree next to a rotatable 3D structure viewer.
Click a species in the tree — the structure morphs to it.

## Usage

```bash
python scripts/interactive_viewer.py structures/aligned.json -o viewer.html --title "ECA3 protein structure across lineages"
```

Open `viewer.html` in any browser (needs internet for the Three.js CDN).

## What's new in v9

- **Ribbon style** — helices render as wide flat ribbon boxes, strands as flat
  arrows (shaft + cone head), loops as the thin smooth tube. SS-class base colors
  (purple helix / yellow strand / grey loop) are tinted per segment by smoothed
  Kyte-Doolittle hydropathy (blue hydrophilic → white → orange hydrophobic).
- **Consistent medium lighting** — hemisphere (0.55) + directional (0.45) lights
  attached to the camera at fixed medium intensity. Illumination stays fixed
  relative to the viewer while orbiting, and shading never blows out or goes dark.
- **Stability-masked morphing** — switching species computes per-residue morph
  weights (Cα displacement + secondary-structure class change, same rule as the
  GIF renderer): stable residues that keep their SSE class hold their position
  (a helix that is a helix in both species does not reshuffle), changing residues
  travel with smoothstep easing, and insertions grow in from the source centroid.
  No more whole-protein reshuffle with an abrupt end.
- **Membrane visualization removed** — the predicted membrane slab planes are gone.

## Tree

SVG phylogenetic tree from live NCBI taxonomy; leaves at uniform
depth; species labels colored #e6edf3 on the dark theme; internal
labels only at branching nodes (≥2 children) in #58a6ff italic, with
wider row spacing. Click handlers on the whole leaf row (circle, text,
and hit area).

## Morphing

O(n) with precomputed residue→reference maps; per-pair stability weights
computed once at morph start; the target tube geometry is built once and
updated **in place** (vertex positions refreshed per frame, no geometry
rebuild) — smooth 36-frame morphs even for 1000+ residue proteins.
The morph math mirrors `src/proteinlib.py` (StabilityMask / morph_frames)
so the web viewer and the GIF renderer behave identically.

## Input

`aligned.json` from the structure-aligner skill v2 (structures should carry
`seq` and `ss`; if absent they are derived from the CA trace). Both assets
(`scripts/interactive_viewer.py` + `scripts/viewer_template.html`)
must sit side by side.

## Validation

- Node logic tests with THREE stubs: weights are 0 for stable same-SSE
  residues, 1 for displaced / SS-changed / inserted residues; morph frames
  start exactly at the source and converge exactly to the target; frozen
  residues hold the source position at mid-morph; insertions interpolate
  from the source centroid.
- End-to-end HTML generation on the ECA3 test pair; payload verified
  (0 membrane references, morph helpers present, SS/hydro/pairs consistent).

## Dependencies

Python 3.8+: requests, numpy, biopython, plus `src/proteinlib.py`.
Browser: WebGL.
