#!/usr/bin/env python3
"""Generate an interactive evolutionary structure viewer (single HTML). (v13)

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

    payload = {
        "title": args.title,
        "reference": ref["accession"],
        "lineageOrder": LINEAGE_ORDER,
        "lineageColors": LINEAGE_COLORS,
        "tree": {"nodes": nodes, "edges": edges, "leaves": leaves},
        "structures": out_structs,
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
