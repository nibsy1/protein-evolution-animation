#!/usr/bin/env python3
"""Generate an interactive evolutionary structure viewer (single HTML). (v17)

v17 changes:
    - All sequence rows always drawn (no off-screen culling) — every row is
      active and visible when scrolled into view via native overflow:auto.
    - Draggable left panel resize handle (180–500px, ew-resize cursor).

Usage:
    python interactive_viewer.py structures/aligned.json -o viewer.html

Reads aligned.json (structure-aligner skill v2) and embeds:
  - an interactive phylogenetic tree (SVG) of the species, colored by
    lineage; clicking a species morphs the structure
  - a Three.js cartoon-style 3D viewer (works from CA-only traces):
    round tube with variable thickness (helix thick, strand medium,
    loop thin), smooth transitions, rotatable/zoomable
  - per-residue PURE HYDRopathy coloring (blue hydrophilic -> white ->
    red hydrophobic), no SS-class tint
  - consistent MEDIUM lighting: hemisphere + directional lights attached
    to the camera at fixed medium intensity, so shading is identical while
    orbiting
  - STABILITY-MASKED morphing: when switching species, residues that are
    spatially stable AND keep their secondary-structure class stay put
    (a helix that is a helix in both species does not reshuffle); only
    genuinely changing regions move, insertions fade in, and timing is
    smoothstep-eased
  - v16 NATIVE SCROLL + RESIZE + RESIDUE HIGHLIGHT: sequence panel uses
    native browser scrolling (all rows always visible), draggable resize
    handle (60-600px), click individual amino acid -> 3D gold sphere+ring
    highlight that persists across species switches (auto-removed if the
    residue doesn't exist in the next structure).
  - v15 SEQUENCE ALIGNMENT PANEL: MSA with BLOSUM62 similarity coloring
    (red=dissimilar -> white -> green=identical), selected species row
    highlighted with smooth cross-fade transition, click-to-select,
    collapsible panel at bottom of viewer.
  - v14 UNION-SPACE MORPH + DoubleSide: handles both insertions AND
    deletions. Insertions fade in (radius 0->1), deletions collapse and
    fade out (1->0). DoubleSide material fixes "see inside tube" issue.
    No more spline artifacts from clumped coincident points near
    conserved helices when C/N-terminal tails differ in length.
  - v13 STANDARD CARTOON: round tube with variable thickness (helix=2.2,
    strand=1.6, loop=0.7 radius), smooth 3-point-averaged transitions at
    SS boundaries. No elliptical cross-sections — clean PyMOL/VMD look.
  - v12 PARALLEL-TRANSPORT FRAMES + LOCAL INSERTION GROWTH: the tube's
    cross-section frame follows the curve via parallel transport (fixes
    oval/twisted helix ribbons), and insertions grow from the nearest
    paired neighbor's morphing position instead of the global centroid
    (fixes "expanding from inner center" on distant-species transitions).
  - v11 VARIABLE-RADIUS TUBE + NO REBUILD AT MORPH END: the tube's
    elliptical cross-section (width/height) cross-fades with SS changes
    during the morph, and at f=1 the tube snaps in place (no clearGroup /
    showStructure) so the camera orientation is preserved — no more
    broken-looking morph or orientation jump at the end.
  - N-/C-termini are marked with labeled spheres on every structure.
  - NO membrane visualization (removed in v7; template asset restored in v8)

The output is one self-contained HTML file (Three.js from CDN).
"""
import argparse
import json
import os
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np
import requests

_HERE = os.path.dirname(os.path.abspath(__file__))
for _cand in (os.path.join(_HERE, "..", "..", "..", "src"),
              os.path.join(_HERE, "..", "..", "src"),
              os.path.join(_HERE, "src")):
    if os.path.isdir(_cand) and _cand not in sys.path:
        sys.path.insert(0, os.path.abspath(_cand))
try:
    import proteinlib as pl
except ImportError:
    pl = None

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

LINEAGE_ORDER = ["Green Algae", "Bryophytes", "Ferns",
                 "Gymnosperms", "Angiosperms"]
LINEAGE_COLORS = {
    "Green Algae": "#2ca02c", "Bryophytes": "#8c564b",
    "Ferns": "#ff7f0e", "Gymnosperms": "#9467bd",
    "Angiosperms": "#d62728",
}

KD = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9,
      'A': 1.8, 'G': -0.4, 'T': -0.7, 'S': -0.8, 'W': -0.9,
      'P': -1.6, 'Y': -1.3, 'H': -3.2, 'E': -3.5, 'Q': -3.5,
      'D': -3.5, 'N': -3.5, 'K': -3.9, 'R': -4.5, 'X': 0.0}

AA1TO3 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
          'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
          'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
          'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
          'X': 'GLY'}

# BLOSUM62 substitution matrix as a flat dict for sequence similarity coloring.
# Scores range from -4 (very dissimilar) to +11 (identical, e.g. W-W).
# Used to color aligned sequence residues by similarity to the reference.
BLOSUM62 = {
    'A': {'A':4,'R':-1,'N':-2,'D':-2,'C':0,'Q':-1,'E':-1,'G':0,'H':-2,'I':-1,'L':-1,'K':-1,'M':-1,'F':-2,'P':-1,'S':1,'T':0,'W':-3,'Y':-2,'V':0},
    'R': {'A':-1,'R':5,'N':0,'D':-2,'C':-3,'Q':1,'E':0,'G':-2,'H':0,'I':-3,'L':-2,'K':2,'M':-1,'F':-3,'P':-2,'S':-1,'T':-1,'W':-3,'Y':-2,'V':-3},
    'N': {'A':-2,'R':0,'N':6,'D':1,'C':-3,'Q':0,'E':0,'G':0,'H':1,'I':-3,'L':-3,'K':0,'M':-2,'F':-3,'P':-2,'S':1,'T':0,'W':-4,'Y':-2,'V':-3},
    'D': {'A':-2,'R':-2,'N':1,'D':6,'C':-3,'Q':0,'E':2,'G':-1,'H':-1,'I':-3,'L':-3,'K':-1,'M':-3,'F':-3,'P':-1,'S':0,'T':-1,'W':-4,'Y':-3,'V':-3},
    'C': {'A':0,'R':-3,'N':-3,'D':-3,'C':9,'Q':-3,'E':-4,'G':-3,'H':-3,'I':-1,'L':-1,'K':-3,'M':-1,'F':-2,'P':-3,'S':-1,'T':-1,'W':-2,'Y':-2,'V':-1},
    'Q': {'A':-1,'R':1,'N':0,'D':0,'C':-3,'Q':5,'E':2,'G':-2,'H':0,'I':-3,'L':-2,'K':1,'M':0,'F':-3,'P':-1,'S':0,'T':-1,'W':-2,'Y':-1,'V':-2},
    'E': {'A':-1,'R':0,'N':0,'D':2,'C':-4,'Q':2,'E':5,'G':-2,'H':0,'I':-3,'L':-3,'K':1,'M':-2,'F':-3,'P':-1,'S':0,'T':-1,'W':-3,'Y':-2,'V':-2},
    'G': {'A':0,'R':-2,'N':0,'D':-1,'C':-3,'Q':-2,'E':-2,'G':6,'H':-2,'I':-4,'L':-4,'K':-2,'M':-3,'F':-3,'P':-2,'S':0,'T':-2,'W':-2,'Y':-3,'V':-3},
    'H': {'A':-2,'R':0,'N':1,'D':-1,'C':-3,'Q':0,'E':0,'G':-2,'H':8,'I':-3,'L':-2,'K':-1,'M':-2,'F':-1,'P':-2,'S':-1,'T':-2,'W':-2,'Y':2,'V':-3},
    'I': {'A':-1,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-3,'E':-3,'G':-4,'H':-3,'I':4,'L':2,'K':-3,'M':1,'F':0,'P':-3,'S':-2,'T':-1,'W':-3,'Y':-1,'V':3},
    'L': {'A':-1,'R':-2,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-3,'G':-4,'H':-2,'I':2,'L':4,'K':-2,'M':2,'F':0,'P':-3,'S':-2,'T':-1,'W':-2,'Y':-1,'V':1},
    'K': {'A':-1,'R':2,'N':0,'D':-1,'C':-3,'Q':1,'E':1,'G':-2,'H':-1,'I':-3,'L':-2,'K':5,'M':-1,'F':-3,'P':-1,'S':0,'T':-1,'W':-3,'Y':-2,'V':-2},
    'M': {'A':-1,'R':-1,'N':-2,'D':-3,'C':-1,'Q':0,'E':-2,'G':-3,'H':-2,'I':1,'L':2,'K':-1,'M':5,'F':0,'P':-2,'S':-1,'T':-1,'W':-1,'Y':-2,'V':1},
    'F': {'A':-2,'R':-3,'N':-3,'D':-3,'C':-2,'Q':-3,'E':-3,'G':-3,'H':-1,'I':0,'L':0,'K':-3,'M':0,'F':6,'P':-4,'S':-2,'T':-1,'W':1,'Y':3,'V':-1},
    'P': {'A':-1,'R':-2,'N':-2,'D':-1,'C':-3,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-3,'L':-3,'K':-1,'M':-2,'F':-4,'P':7,'S':-1,'T':-2,'W':-4,'Y':-3,'V':-2},
    'S': {'A':1,'R':-1,'N':1,'D':0,'C':-1,'Q':0,'E':0,'G':0,'H':-1,'I':-2,'L':-2,'K':0,'M':-1,'F':-2,'P':-1,'S':4,'T':1,'W':-3,'Y':-2,'V':-2},
    'T': {'A':0,'R':-1,'N':0,'D':-1,'C':-1,'Q':-1,'E':-1,'G':-2,'H':-2,'I':-1,'L':-1,'K':-1,'M':-1,'F':-1,'P':-2,'S':1,'T':5,'W':-2,'Y':-2,'V':0},
    'W': {'A':-3,'R':-3,'N':-4,'D':-4,'C':-2,'Q':-2,'E':-3,'G':-2,'H':-2,'I':-3,'L':-2,'K':-3,'M':-1,'F':1,'P':-4,'S':-3,'T':-2,'W':11,'Y':2,'V':-3},
    'Y': {'A':-2,'R':-2,'N':-2,'D':-3,'C':-2,'Q':-1,'E':-2,'G':-3,'H':2,'I':-1,'L':-1,'K':-2,'M':-2,'F':3,'P':-3,'S':-2,'T':-2,'W':2,'Y':7,'V':-1},
    'V': {'A':0,'R':-3,'N':-3,'D':-3,'C':-1,'Q':-2,'E':-2,'G':-3,'H':-3,'I':3,'L':1,'K':-2,'M':1,'F':-1,'P':-2,'S':-2,'T':0,'W':-3,'Y':-1,'V':4},
}


def blosum_score(a, b):
    """BLOSUM62 substitution score between two amino acids.
    Returns 0 for unknown residues (X, gaps, etc.)."""
    if a not in BLOSUM62 or b not in BLOSUM62:
        return 0
    return BLOSUM62[a][b]


def build_msa(structs, ref_seq, ref_acc):
    """Build a multiple sequence alignment (MSA) anchored to the reference.

    Uses a reference-anchored star alignment: the reference sequence defines
    the column backbone (one column per reference residue, no gaps in the
    reference row). Insertions from other sequences are placed in additional
    columns between reference columns. The reference alignment row has NO
    gaps, making it easy to map structure residue indices to MSA columns.

    Returns:
      - msa_seqs: list of (species, aligned_seq) in MSA column space
      - msa_cols: number of columns in the MSA
      - ref_col_map: list mapping reference residue index -> MSA column
        (identity mapping since reference defines the backbone, but kept
        for generality; ref residue i is at column ref_col_map[i])
    """
    n_ref = len(ref_seq)

    # For each structure, get alignment pairs to reference
    all_pairs = {}
    for s in structs:
        seq = s.get("seq") or ""
        if not seq:
            all_pairs[s["species"]] = []
            continue
        if s["accession"] == ref_acc:
            all_pairs[s["species"]] = [(i, i) for i in range(min(len(seq), n_ref))]
        else:
            pairs = pl.sequence_pairs(seq, ref_seq)
            all_pairs[s["species"]] = pairs

    # Build seq_to_ref mapping
    seq_to_ref = {}
    for sp, pairs in all_pairs.items():
        seq_to_ref[sp] = {si: ri for si, ri in pairs}

    # For each non-reference sequence, find unpaired residues (insertions)
    # and group them into blocks, noting which reference position they
    # insert after.
    seq_insertions = {}  # species -> list of (seq_start, seq_end, after_ref_idx)
    for s in structs:
        sp = s["species"]
        if sp == ref_acc or sp not in seq_to_ref:
            continue
        seq = s.get("seq") or ""
        if not seq:
            continue
        paired = set(seq_to_ref[sp].keys())
        unpaired = [i for i in range(len(seq)) if i not in paired]
        insertions = []
        i = 0
        while i < len(unpaired):
            start = unpaired[i]
            end = start
            while i + 1 < len(unpaired) and unpaired[i + 1] == end + 1:
                i += 1
                end = unpaired[i]
            # Find the reference position this block inserts after
            before = -1
            for j in range(start - 1, -1, -1):
                if j in seq_to_ref[sp]:
                    before = seq_to_ref[sp][j]
                    break
            insertions.append((start, end, before))
            i += 1
        seq_insertions[sp] = insertions

    # Collect all insertion points: after_ref_idx -> list of (species, seq_start, seq_end)
    # Each insertion block from each species gets its own set of columns.
    # Multiple species inserting at the same point get separate column groups.
    ins_points = {}  # after_ref_idx -> list of (species, seq_start, seq_end)
    for sp, inss in seq_insertions.items():
        for seq_start, seq_end, after_ref in inss:
            if after_ref not in ins_points:
                ins_points[after_ref] = []
            ins_points[after_ref].append((sp, seq_start, seq_end))

    # Build column layout:
    # - Reference columns: one per reference residue (ref row has no gaps)
    # - Insertion columns: between reference columns, one per inserted residue
    #   from each species (each species gets its own columns)
    # Column structure: {type: 'ref', ref_idx: int} or {type: 'ins', species: str, seq_idx: int}
    columns = []
    # Pre-insertion columns (insertions before first ref residue, after_ref = -1)
    if -1 in ins_points:
        for sp, ss, se in sorted(ins_points[-1]):
            for k in range(ss, se + 1):
                columns.append({"type": "ins", "species": sp, "seq_idx": k})
    # Reference columns with interleaved insertion columns
    for ri in range(n_ref):
        columns.append({"type": "ref", "ref_idx": ri})
        if ri in ins_points:
            for sp, ss, se in sorted(ins_points[ri]):
                for k in range(ss, se + 1):
                    columns.append({"type": "ins", "species": sp, "seq_idx": k})
    # Post-insertion columns (insertions after last ref residue)
    if n_ref in ins_points:
        for sp, ss, se in sorted(ins_points[n_ref]):
            for k in range(ss, se + 1):
                columns.append({"type": "ins", "species": sp, "seq_idx": k})

    n_cols = len(columns)

    # Build aligned sequences
    msa_seqs = []
    for s in structs:
        sp = s["species"]
        seq = s.get("seq") or ""
        if not seq:
            msa_seqs.append((sp, "-" * n_cols))
            continue
        # Build column -> residue mapping
        col_to_res = {}
        if sp in seq_to_ref:
            # Map reference columns
            for ci, col in enumerate(columns):
                if col["type"] == "ref":
                    ri = col["ref_idx"]
                    # Find seq_idx that maps to this ref_idx
                    for si, rj in seq_to_ref[sp].items():
                        if rj == ri:
                            col_to_res[ci] = seq[si]
                            break
        # Map insertion columns
        for ci, col in enumerate(columns):
            if col["type"] == "ins" and col["species"] == sp:
                col_to_res[ci] = seq[col["seq_idx"]]
        # Build the aligned string
        aligned = []
        for ci in range(n_cols):
            aligned.append(col_to_res.get(ci, "-"))
        msa_seqs.append((sp, "".join(aligned)))

    # ref_col_map: reference residue i is at the column where type='ref' and ref_idx==i
    ref_col_map = []
    for ri in range(n_ref):
        for ci, col in enumerate(columns):
            if col["type"] == "ref" and col["ref_idx"] == ri:
                ref_col_map.append(ci)
                break

    return msa_seqs, n_cols, ref_col_map


def compute_similarity_scores(msa_seqs, ref_species):
    """Compute per-column BLOSUM62 similarity scores for each sequence.

    For each sequence, for each MSA column, compute the BLOSUM62 score
    between the sequence's residue and the reference residue at that column.
    Returns: list of (species, [score_per_column]) where scores are normalized
    to [0, 1] (0 = dissimilar, 1 = identical).
    """
    ref_aligned = None
    for sp, aln in msa_seqs:
        if sp == ref_species:
            ref_aligned = aln
            break
    if ref_aligned is None:
        return [(sp, [0.5] * len(aln)) for sp, aln in msa_seqs]

    n_cols = len(ref_aligned)
    results = []
    for sp, aln in msa_seqs:
        scores = []
        for ci in range(n_cols):
            ref_aa = ref_aligned[ci]
            seq_aa = aln[ci]
            if ref_aa == "-" or seq_aa == "-":
                scores.append(None)
            else:
                raw = blosum_score(seq_aa, ref_aa)
                self_score = blosum_score(ref_aa, ref_aa)
                if self_score > 0:
                    norm = (raw + 4) / (self_score + 4)
                    norm = max(0, min(1, norm))
                else:
                    norm = 0.5
                scores.append(round(norm, 3))
        results.append((sp, scores))
    return results


def hydrophobicity(seq, window=9):
    """Smoothed Kyte-Doolittle hydropathy per residue."""
    raw = np.array([KD.get(a, 0.0) for a in seq], dtype=float)
    if len(raw) == 0:
        return raw
    half = window // 2
    sm = np.zeros_like(raw)
    for i in range(len(raw)):
        lo = max(0, i - half)
        hi = min(len(raw), i + half + 1)
        sm[i] = raw[lo:hi].mean()
    return sm


def to_pdb_text(struct):
    """CA-trace PDB with column-exact ATOM and HELIX/SHEET records."""
    coords = np.asarray(struct["coords"], dtype=float)
    seq = struct.get("seq") or "A" * len(coords)
    ss = struct.get("ss") or (pl.assign_secondary(coords)
                              if pl else ['C'] * len(coords))
    ss = "".join(ss)
    lines = []
    nhelix = 0
    nsheet = 0
    i = 0
    while i < len(ss):
        if ss[i] == 'H':
            j = i
            while j < len(ss) and ss[j] == 'H':
                j += 1
            if j - i >= 4:
                nhelix += 1
                r1 = AA1TO3.get(seq[i], 'GLY')
                r2 = AA1TO3.get(seq[j - 1], 'GLY')
                lines.append(f"HELIX  {nhelix:3d} {nhelix:3d} "
                             f"{r1:>3s} A {i + 1:4d}  {r2:>3s} "
                             f"A {j:4d}  1")
            i = j
        else:
            i += 1
    i = 0
    while i < len(ss):
        if ss[i] == 'E':
            j = i
            while j < len(ss) and ss[j] == 'E':
                j += 1
            if j - i >= 3:
                nsheet += 1
                r1 = AA1TO3.get(seq[i], 'GLY')
                r2 = AA1TO3.get(seq[j - 1], 'GLY')
                lines.append(f"SHEET  {nsheet:3d} {nsheet:3d} 1 "
                             f"{r1:>3s} A{i + 1:4d}  {r2:>3s} "
                             f"A{j:4d}  0")
            i = j
        else:
            i += 1
    for k, xyz in enumerate(coords, start=1):
        x, y, z = xyz
        res = AA1TO3.get(seq[k - 1] if k - 1 < len(seq) else 'A', 'GLY')
        lines.append(f"ATOM  {k:5d}  CA  {res:>3s} A{k:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00")
    lines.append("TER")
    lines.append("END")
    return chr(10).join(lines)


def fetch_lineage(species, session, cache):
    if species in cache:
        return cache[species]
    for _ in range(3):
        try:
            r = session.get(f"{BASE}/esearch.fcgi",
                            params={"db": "taxonomy", "term": species,
                                    "retmax": 1}, timeout=30)
            root = ET.fromstring(r.text)
            ids = [e.text for e in root.findall(".//IdList/Id")]
            if not ids:
                break
            r2 = session.get(f"{BASE}/efetch.fcgi",
                             params={"db": "taxonomy", "id": ids[0]},
                             timeout=30)
            root2 = ET.fromstring(r2.text)
            ranks = [(el.findtext("Rank"), el.findtext("ScientificName"))
                     for el in root2.findall(".//Taxon/LineageEx/Taxon")]
            ranks = [name for rank, name in ranks
                     if rank in ("superkingdom", "kingdom", "phylum",
                                 "subphylum", "class", "subclass", "order",
                                 "suborder", "family", "subfamily", "genus")]
            cache[species] = ranks
            time.sleep(0.35)
            return ranks
        except Exception:
            time.sleep(1)
    cache[species] = []
    return []


def build_tree(species_lineages):
    root = {"name": "root", "children": {}, "leaf": None}
    for sp, ranks in species_lineages.items():
        node = root
        for name in ranks:
            if name not in node["children"]:
                node["children"][name] = {"name": name, "children": {},
                                          "leaf": None}
            node = node["children"][name]
        node["children"][sp] = {"name": sp, "children": {},
                                "leaf": sp}
    return root


def layout_tree(root, order):
    leaves = []

    def collect(node):
        if node["leaf"]:
            leaves.append(node["leaf"])
            return
        for ch in node["children"].values():
            collect(ch)
    collect(root)

    rank0 = {sp: i for i, sp in enumerate(leaves)}
    def leaf_key(sp):
        lin = SPECIES_LINEAGE.get(sp, "")
        r = order.index(lin) if lin in order else 99
        return (r, rank0[sp])
    leaves.sort(key=leaf_key)

    pos_y = {sp: i for i, sp in enumerate(leaves)}
    nodes = []
    edges = []

    def maxdepth(node, d=0):
        if node["leaf"] or not node["children"]:
            return d
        return max(maxdepth(ch, d + 1)
                   for ch in node["children"].values())
    leaf_depth = maxdepth(root)

    def place(node, depth):
        if node["leaf"]:
            y = pos_y[node["leaf"]]
            nodes.append({"x": leaf_depth, "y": y,
                          "label": node["leaf"], "leaf": True})
            return y, leaf_depth
        ys = []
        child_ys = []
        for ch in node["children"].values():
            cy, cmax = place(ch, depth + 1)
            child_ys.append((cy, cmax))
            ys.append(cy)
        my = float(np.mean(ys))
        nodes.append({"x": depth, "y": my, "label": node["name"],
                      "leaf": False})
        for cy, cmax in child_ys:
            edges.append([depth, my, cmax, cy])
        return my, depth
    place(root, 0)
    return nodes, edges, leaves


SPECIES_LINEAGE = {}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("aligned_json")
    p.add_argument("-o", "--output", default="viewer.html")
    p.add_argument("--title", default="Protein evolution viewer")
    args = p.parse_args()

    if pl is None:
        sys.exit("proteinlib.py not found - put it next to the skills tree "
                 "(src/proteinlib.py) or on PYTHONPATH")

    with open(args.aligned_json) as f:
        data = json.load(f)
    structs = data["structures"]
    if not structs:
        sys.exit("No structures in aligned.json")

    for s in structs:
        SPECIES_LINEAGE[s["species"]] = s.get("lineage", "")

    print(f"[tree] fetching taxonomy for {len(structs)} species ...",
          file=sys.stderr)
    session = requests.Session()
    cache = {}
    lineages = {}
    for s in structs:
        lineages[s["species"]] = fetch_lineage(s["species"], session, cache)
    tree = build_tree(lineages)
    nodes, edges, leaves = layout_tree(tree, LINEAGE_ORDER)

    print("[structures] computing SS / hydropathy ...", file=sys.stderr)
    out_structs = []
    ref_acc = data.get("reference")
    ref = next((s for s in structs if s["accession"] == ref_acc),
               structs[0])
    ref_coords = np.asarray(ref["coords"], dtype=float)
    ref_seq = ref.get("seq") or "A" * len(ref_coords)

    for s in structs:
        coords = np.asarray(s["coords"], dtype=float)
        seq = s.get("seq") or "A" * len(coords)
        n_coord = len(coords)
        coord_seq = seq[:n_coord] if len(seq) >= n_coord else seq + "X" * (n_coord - len(seq))
        h = hydrophobicity(coord_seq)
        ss = s.get("ss") or "".join(pl.assign_secondary(coords))
        pairs = pl.sequence_pairs(seq, ref_seq)
        pairs = [[a, b] for a, b in pairs
                 if a < len(coords) and b < len(ref_coords)]
        out_structs.append({
            "accession": s["accession"], "species": s["species"],
            "lineage": s.get("lineage", ""),
            "source": s.get("source", ""),
            "core_identity": s.get("core_identity"),
            "n_residues": len(coords),
            "ss": ss,
            "pdb": to_pdb_text({"coords": coords.tolist(), "seq": seq,
                                "ss": ss}),
            "hydro": [round(float(x), 2) for x in h],
            "ref_pairs": pairs,
            "ref_coords": None,
        })
    ref_pts = [[round(float(c), 3) for c in xyz] for xyz in ref_coords]
    for s in out_structs:
        s["ref_coords"] = ref_pts
        if s["accession"] == ref["accession"]:
            s["ref_pairs"] = [[i, i] for i in range(len(ref_coords))]

    # v15: Build MSA and BLOSUM62 similarity scores for sequence panel
    print("[sequences] building MSA and similarity scores ...", file=sys.stderr)
    msa_seqs, n_msa_cols, ref_col_map = build_msa(structs, ref_seq, ref["accession"])
    sim_scores = compute_similarity_scores(msa_seqs, ref["species"])
    # Pack as compact arrays: {species: {aln: "...", sim: [0.5, null, ...]}}
    seq_data = {}
    for sp, aln in msa_seqs:
        sim = next((s for ssp, s in sim_scores if ssp == sp), [0.5] * len(aln))
        seq_data[sp] = {"aln": aln, "sim": sim}
    # Overall identity per species (mean of non-null similarity scores)
    identity = {}
    for sp, aln in msa_seqs:
        sim = seq_data[sp]["sim"]
        valid = [s for s in sim if s is not None]
        identity[sp] = round(sum(valid) / max(len(valid), 1), 3) if valid else 0

    payload = {
        "title": args.title,
        "reference": ref["accession"],
        "referenceSpecies": ref["species"],
        "lineageOrder": LINEAGE_ORDER,
        "lineageColors": LINEAGE_COLORS,
        "tree": {"nodes": nodes, "edges": edges, "leaves": leaves},
        "structures": out_structs,
        "sequences": {
            "msaCols": n_msa_cols,
            "refColMap": ref_col_map,
            "data": seq_data,
            "identity": identity,
        },
    }

    template = open(os.path.join(os.path.dirname(
        os.path.abspath(__file__)), "viewer_template.html")).read()
    html = (template.replace("__DATA__", json.dumps(payload))
                .replace("__TITLE__", args.title))
    with open(args.output, "w") as f:
        f.write(html)
    print(f"[done] viewer -> {args.output} "
          f"({os.path.getsize(args.output) / 1e6:.1f} MB)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
