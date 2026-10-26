---
name: interactive-viewer
description: Generates an interactive HTML viewer: clickable phylogenetic tree +
  rotatable 3D ribbon-style structure (Three.js, works from CA-only traces) with
  consistent medium camera-attached lighting, pure hydropathy coloring (blue
  hydrophilic → white → red hydrophobic), and a stability-masked smooth morph
  on species selection. No membrane visualization. Single self-contained HTML
  from aligned.json (structure-aligner v2).
  v17: all sequence rows always drawn (no off-screen culling — every row is
  always active and visible when scrolled into view), draggable left panel
  resize handle (180–500px, ew-resize cursor at right edge of tree panel).
  v16: native scrolling, draggable bottom panel resize handle (60–600px),
  click individual amino acid → 3D gold sphere+ring highlight that persists
  across species switches (auto-removed if residue doesn't exist in next structure).
  v15: sequence alignment panel — MSA with BLOSUM62 similarity coloring
  (red=dissimilar, white=neutral, green=identical), selected species row
  highlighted with smooth transition, click-to-select, collapsible panel.
  v14: union-space morph handles both insertions AND deletions. Insertions
  fade in (radius scales 0→1), deletions collapse and fade out (1→0).
  DoubleSide material fixes "see inside tube" issue. No more clumped
  coincident points distorting the CatmullRom spline near conserved helices.
  v13: standard cartoon — round tube with variable thickness (helix thick,
  strand medium, loop thin), smooth transitions at SS boundaries. No more
  elliptical cross-sections. Clean PyMOL/VMD-style cartoon look.
  v12: parallel-transport frames fix oval/twisted helix ribbons; insertions
  grow from nearest paired neighbor (not global centroid) — no more
  "expanding from inner center" on distant species transitions.
  v11: variable-radius ribbon tube (helix wide/flat, strand wide, loop thin)
  replaces separate helix boxes + strand arrows — no more broken morph rendering.
  No geometry rebuild at morph end — no more orientation jump. Pure hydropathy
  coloring. SS legend. N-/C-terminus markers.
version: 17
---
# Interactive Evolution Viewer (v17)

Generates a single self-contained HTML file: an interactive
phylogenetic tree next to a rotatable 3D structure viewer.
Click a species in the tree — the structure morphs to it.

## Usage

```bash
python scripts/interactive_viewer.py structures/aligned.json -o viewer.html --title "ECA3 protein structure across lineages"
```

Open `viewer.html` in any browser (needs internet for the Three.js CDN).

## What's new in v17

- **All sequence rows always drawn** — removed the off-screen row culling in
  `drawSeqPanel()`. Previously, rows outside the visible panel area were
  skipped during canvas redraw, which meant they stayed invisible (stale
  pixels) until you clicked them. Now every row is always drawn on every
  redraw, so all sequences are always active and visible when scrolled
  into view. The canvas is sized to the full content height and the
  container scrolls natively with `overflow:auto`.
- **Draggable left panel border** — a resize handle (`#leftresize`, 6px
  vertical bar at the right edge of the phylogenetic tree panel,
  `ew-resize` cursor) lets the user adjust the left panel width between
  180px and 500px. The 3D renderer resizes in real time via `resizeRenderer()`.
  The left panel now has `flex-shrink:0` so it doesn't get compressed by
  the center panel.

## What's new in v16

- **Native scrolling** — the sequence panel container uses native browser
  `overflow:auto` scrolling (`scrollTop`/`scrollLeft`) instead of a manual
  `scrollY` offset. All sequence rows are always active and visible when
  scrolled into view — no more invisible rows that only appear on click.
- **Adjustable panel height** — a drag handle (`#seqresize`, 6px bar at the
  top of the sequence panel, `n-resize` cursor) lets the user resize the
  panel between 60px and 600px. Dragging up shrinks, dragging down grows.
  The 3D renderer and sequence canvas resize in real time.
- **Residue click → 3D highlight** — clicking an individual amino acid in
  the sequence panel places a pulsing gold sphere + camera-facing torus ring
  at that residue's Cα position in the 3D structure. The highlight is added
  to the scene root (not the structure group) so it persists during morphs.
  When switching species, the highlight automatically follows to the
  corresponding residue (mapped via `ref_pairs`); if the residue doesn't
  exist in the new structure (insertion/deletion), the highlight is removed.
  Clicking the same residue again deselects it. A gold selection box is
  drawn behind the residue in every sequence row, with a vertical connector
  line across all rows.

## What's new in v15

- **Sequence alignment panel** — a third panel at the bottom shows a multiple
  sequence alignment (MSA) of all species, anchored to the reference sequence.
  The reference row has no gaps; insertions from other species appear as
  extra columns between reference columns.
- **BLOSUM62 similarity coloring** — each residue is colored by its BLOSUM62
  substitution score relative to the reference residue at the same column,
  normalized by the reference residue's self-score. Color scale:
  red (dissimilar, e.g. charged vs hydrophobic) → white (neutral) →
  green (identical/similar, e.g. D↔E both acidic, K↔R both basic).
- **Selected species highlight** — the selected species' row is highlighted
  with a blue background and left accent border. Switching species triggers
  a smooth 400ms cross-fade transition between the old and new highlights.
- **Identity percentage** — each row shows the overall BLOSUM62-normalized
  identity to the reference, making it easy to spot outliers (e.g. Picea
  abies at ~60% vs Hordeum vulgare at ~98%).
- **Click-to-select** — clicking a sequence row selects that species
  (triggers the 3D morph). Mouse wheel scrolls horizontally through the
  alignment. Panel is collapsible.

## What's new in v14

- **Union-space morph** — the morph now works in a merged index space that
  includes paired residues, insertions (target-only), AND deletions
  (source-only). Previously, source-only C/N-terminal tails simply vanished
  and all target insertions clumped at a single anchor point, creating
  CatmullRom spline artifacts that distorted nearby conserved helices.
- **Insertions fade in** — new residues start with `radiusScale=0` (invisible)
  and grow to full thickness. Positions are spread along the anchor→target
  direction to prevent spline degeneracy from coincident points.
- **Deletions collapse and fade out** — source-only residues (e.g., a long
  C-terminal tail in Oryza that Phaseolus lacks) stay at their source position
  at morph start, then collapse toward the nearest paired neighbor while
  fading to zero thickness. At morph end they are invisible at the anchor.
- **DoubleSide material** — the tube mesh now renders both sides, fixing the
  "see inside the tubing" issue where back-facing surfaces were invisible.

## What's new in v13

- **Standard cartoon: round variable-thickness tube** — replaced the elliptical
  cross-sections with a clean round tube whose thickness varies by SS class:
  helices thick (2.2 Å radius), strands medium (1.6 Å), loops thin (0.7 Å).
  This is the standard cartoon representation used by PyMOL, VMD, and ChimeraX.
- **Smooth thickness transitions** — a 3-point moving average blends the radius
  at SS boundaries (helix→loop, strand→coil) so there are no hard steps.
- **Simplified SS legend** — shows three dots of decreasing size (thick=helix,
  medium=strand, thin=loop) instead of shape indicators.

## What's new in v12

- **Parallel-transport frames** — the ribbon tube's cross-section frame now
  follows the curve using parallel transport (projecting the previous frame's
  normal onto the plane perpendicular to the new tangent). This eliminates the
  oval/twisted appearance of helix ribbons that occurred when the binormal was
  locked to the XZ plane (`B = (N.z, 0, -N.x)`), which broke down for helices
  running along the Y axis.
- **Insertions grow from nearest paired neighbor** — when a species has residues
  not present in the source (insertions, up to 34% for distant pairs like
  Chlamydomonas ↔ Arabidopsis), they now grow from the morphing position of the
  nearest paired residue along the chain, not from the global centroid. This
  eliminates the "expanding from inner center outwards" effect on transitions
  between dissimilar structures. Same fix applied to `proteinlib.morph_frames`
  so the GIF renderer matches.

## What's new in v11

- **Variable-radius ribbon tube** — a single TubeGeometry with per-ring elliptical
  cross-sections replaces the old combination of uniform tube + separate helix
  ribbon boxes + strand arrows. Helices render wide and flat (3.0×1.0 Å),
  strands wide (2.4×0.9 Å), loops thin (0.9×0.9 Å). This eliminates the sync
  problem where boxes stayed at source positions while the tube morphed
  (the "broken" look during animation).
- **No geometry rebuild at morph end** — `morphTick()` at f=1 snaps the tube
  to exact target coords/colors in place instead of calling `clearGroup()` +
  `showStructure()`. The camera orientation is preserved — no more jump to
  a new orientation when the animation finishes.
- **Pure hydropathy coloring** — residue colors are mapped directly from
  smoothed Kyte-Doolittle hydropathy (blue hydrophilic → white → red hydrophobic)
  with no SS-class tint. The ribbon width/shape still conveys secondary structure.
- **SS legend** — small legend at bottom-left shows the three ribbon shapes
  (wide flat bar = helix, medium bar = strand, thin dot = loop).
- **Consistent medium lighting** — hemisphere (0.55) + directional (0.45) lights
  attached to the camera at fixed medium intensity. Illumination stays fixed
  relative to the viewer while orbiting, and shading never blows out or goes dark.
- **Stability-masked morphing** — switching species computes per-residue morph
  weights (Cα displacement + secondary-structure class change, same rule as the
  GIF renderer): stable residues that keep their SSE class hold their position
  (a helix that is a helix in both species does not reshuffle), changing residues
  travel with smoothstep easing, and insertions grow in from the source centroid.
  The tube cross-section (SS radius) also cross-fades during the morph.
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
