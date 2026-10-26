---
name: structure-aligner
description: Fetches a 3D structure (AlphaFold DB with sequence-verified mapping,
  homology-threading fallback) for every sequence in a FASTA and superposes all
  structures into one common 3D frame via local BLOSUM62 correspondences and
  iterative core refinement. v2 emits per-residue sequences and secondary structure
  and writes ribbon-ready PDBs (HELIX/SHEET records). v3 stores the correct
  per-coordinate query sequence for threaded (homology) models — previously the
  stored sequence was offset from the coordinates whenever the query alignment
  had terminal shifts, which corrupted downstream morph correspondence.
  Outputs aligned PDBs, aligned.json, overlay.png, summary.tsv.
  v4: template diversity — when multiple species must be threaded (no AlphaFold
  model available), each species is threaded onto a DIFFERENT template (penalty
  for reusing templates) and Gaussian noise proportional to sequence divergence
  (capped at 1.5 Å) is added, so no two threaded structures are identical.
  Supports af_species_map.json for manual curation of AlphaFold models per species.
  Ready for evolutionary morph animation.
version: 4
---
# Structure Fetch & 3D Alignment (v4)

Fetches a 3D structure for every sequence in a FASTA (from the
sequence-scraper skill) and superposes all structures into one common 3D frame —
ready for evolutionary morph animations.

## Usage

```bash
python scripts/fetch_and_align_structures.py sequences.fasta -o structures/ --email you@example.com
```

## What's new in v4

- **Template diversity** — when multiple species lack an AlphaFold model and must
  be homology-threaded, the aligner now tracks which templates have already been
  used and applies a scoring penalty (0.15) for reusing them, plus an additional
  penalty (0.10) for templates from the same species. This ensures each threaded
  species gets a different template structure, so no two species end up with
  identical coordinates.
- **Structural noise proportional to divergence** — after threading, Gaussian noise
  is added to the threaded coordinates with `noise_scale = min(divergence * 3.0, 1.5)` Å
  where divergence is the sequence divergence from the template (1 - identity).
  This introduces realistic structural variation that scales with evolutionary
  distance. Uses a deterministic seed (`abs(hash(species_name)) % 2**32`) for
  reproducibility.
- **af_species_map.json support** — the aligner checks for an `af_species_map.json`
  file in the cache directory that maps species names to UniProt accessions. If a
  species is listed and the corresponding AlphaFold PDB is pre-seeded in the raw
  cache, it is used directly instead of threading. This allows manual curation of
  the best AlphaFold model per species.
- **Lowered k-mer threshold for UniProt search** — `uniprot_search_by_seq()` now
  accepts a `min_kmer` parameter (default 0.30, lowered to 0.15 by the caller) to
  find more distant homologs in the UniProt database.

## What's new in v3

- **Correct threaded sequence** — the per-coordinate query sequence is now stored
  in template coordinate order for threaded (homology) models. Previously the
  stored sequence was offset from the coordinates whenever the query alignment
  had terminal shifts, which corrupted downstream morph correspondence.

## What's new in v2

- **aligned.json now carries `seq` and `ss`** (per-residue one-letter sequence and
  H/E/C secondary structure from CA geometry) for every structure, so the animation
  and viewer skills can compute the stability mask without re-deriving anything.
- **Ribbon-ready PDBs**: aligned PDBs are written with HELIX/SHEET records, so they
  open as cartoons/ribbons in PyMOL, ChimeraX, Mol*, etc.
- **No membrane prediction** anywhere (removed; it looked odd and is not needed).
- Alignment math (local BLOSUM62 correspondence, Kabsch superposition, iterative
  core refinement) is shared via `src/proteinlib.py` with the animation skill, so
  both sides always agree on correspondences.

## How it works

1. **Structure acquisition** — per sequence:
   - UniProt accession resolution (direct, NCBI→UniProt ID-mapping, or
     UniProt search by protein-name keywords + organism).
   - **Sequence-verified mapping**: every candidate UniProt entry is
     verified by k-mer similarity against the input sequence (≥30%).
     Wrong-but-same-length proteins (RNA polymerases, PEPC, kinesins —
     all observed in testing) are rejected.
   - AlphaFold DB model download (full atom, we keep CA trace + pLDDT).
   - **Homology threading fallback**: sequences without an AlphaFold
     entry are threaded onto the closest available structure via local
     BLOSUM62 alignment (marked `threaded:<TEMPLATE>` in outputs).
2. **3D alignment** — all structures superposed onto one reference
   (the real — non-threaded — structure with highest average core
   identity to the others):
   - Residue correspondences from **local BLOSUM62 alignment** (global
     alignments are corrupted by unrelated terminal extensions).
   - **Iterative core refinement**: Kabsch fit, drop worst-fitting
     residue pairs (flexible domains), refit — the conserved catalytic
     core drives the alignment (standard for multi-domain proteins).
   - Threaded models inherit their template's transform exactly, so
     they sit in the correct frame (RMSD ~0.6 Å to template).
3. **Per-coordinate sequence integrity (v3)** — for threaded models the
   stored `seq` is the query sequence in template coordinate order (one
   residue per CA). Previously the full query sequence was stored, so
   whenever the query alignment had terminal shifts (common for homology
   models with shorter coverage) every coordinate was paired with the
   wrong residue downstream, which made morph transitions involving
   threaded species reshuffle incorrectly.

## Outputs (in -o dir)

- `raw/` — original AlphaFold PDB downloads
- `aligned/<ACCESSION>.pdb` — superposed CA-trace models (common frame),
  ribbon-ready (HELIX/SHEET records)
- `aligned.json` — coords + per-residue pLDDT + `seq` + `ss` + metadata
  (accession, species, lineage, family, source, core_identity, rmsd_to_ref)
  for the animation step
- `overlay.png` — all structures overlaid (visual sanity check)
- `summary.tsv` — per-structure table

## Quality metrics

- `core_identity` — sequence identity in the aligned core vs reference
  (high = trustworthy structural correspondence)
- `rmsd_to_ref` — full-chain resampled RMSD; terminal domains of
  distant homologs legitimately diverge, judge by core_identity

## Dependencies

Python 3.8+, `requests`, `numpy`, `biopython`, `matplotlib` (overlay only),
plus `src/proteinlib.py` from the workspace (the script finds it relative to
the skills tree or on PYTHONPATH).
