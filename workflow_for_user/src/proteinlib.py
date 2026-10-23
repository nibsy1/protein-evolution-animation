#!/usr/bin/env python3
"""Shared library for the protein evolution animation project.

Structure handling (CA-trace PDB parsing/writing), secondary-structure
assignment, local sequence correspondence (BLOSUM62), Kabsch superposition,
and the stability mask that drives the smooth morph transitions.

The stability mask is the core of the new animation behaviour: for every
pair of consecutive structures we compute, per corresponding residue,
  - the C-alpha displacement after superposition, and
  - whether the secondary-structure element (helix/strand/loop) is conserved.
Residues that are both spatially stable and structurally conserved are
FROZEN during the morph (helix stays a helix, in place); only residues
with real structural change are interpolated, and residues present in only
one structure fade in/out. This removes the old "whole protein reshuffles
with an abrupt end" transitions.
"""
import numpy as np

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}
AA1TO3 = {v: k for k, v in AA3TO1.items()}
AA1TO3['X'] = 'GLY'


# ---------------------------------------------------------------- PDB I/O

def parse_pdb_ca(path):
    """Parse CA atoms from a PDB file -> dict(coords, plddts, residues, seq)."""
    coords, plddts, residues = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                residues.append(line[17:20].strip())
                coords.append((float(line[30:38]), float(line[38:46]),
                               float(line[46:54])))
                try:
                    plddts.append(float(line[60:66]))
                except ValueError:
                    plddts.append(50.0)
    seq = "".join(AA3TO1.get(r, 'X') for r in residues)
    return {"coords": np.array(coords, dtype=float),
            "plddts": np.array(plddts, dtype=float),
            "residues": residues, "seq": seq}


def write_ca_pdb(coords, plddts, path, name, ss=None, seq=None):
    """Write a CA-trace PDB, optionally with HELIX/SHEET records.

    HELIX/SHEET records make the file render as ribbons in PyMOL, ChimeraX,
    Mol*, etc. - the ribbon style the animation and viewer use.
    """
    coords = np.asarray(coords, dtype=float)
    lines = [f"REMARK   aligned CA-trace model: {name}"]
    if ss:
        n = 0
        i = 0
        while i < len(ss):
            if ss[i] == 'H':
                j = i
                while j < len(ss) and ss[j] == 'H':
                    j += 1
                if j - i >= 4:
                    n += 1
                    r1 = AA1TO3.get(seq[i] if seq else 'A', 'GLY')
                    r2 = AA1TO3.get(seq[j - 1] if seq else 'A', 'GLY')
                    lines.append(f"HELIX  {n:3d} {n:3d} {r1:>3s} A {i + 1:4d}  "
                                 f"{r2:>3s} A {j:4d}  1")
                i = j
            else:
                i += 1
        n = 0
        i = 0
        while i < len(ss):
            if ss[i] == 'E':
                j = i
                while j < len(ss) and ss[j] == 'E':
                    j += 1
                if j - i >= 3:
                    n += 1
                    r1 = AA1TO3.get(seq[i] if seq else 'A', 'GLY')
                    r2 = AA1TO3.get(seq[j - 1] if seq else 'A', 'GLY')
                    lines.append(f"SHEET  {n:3d} {n:3d} 1 {r1:>3s} A{i + 1:4d}  "
                                 f"{r2:>3s} A{j:4d}  0")
                i = j
            else:
                i += 1
    for k, (xyz, p) in enumerate(zip(coords, plddts), start=1):
        x, y, z = xyz
        res = AA1TO3.get(seq[k - 1] if seq and k - 1 < len(seq) else 'A', 'GLY')
        lines.append(f"ATOM  {k:5d}  CA  {res:>3s} A{k:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 {p:6.2f}")
    lines.append("TER")
    lines.append("END")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


# ------------------------------------------------------- secondary structure

def assign_secondary(coords):
    """Assign H/E/C per residue from CA geometry.

    Helix: CA(i)..CA(i+3) distance ~4.5-5.8 A. Strand: CA(i)..CA(i+5)
    distance ~13-15 A (parallel/antiparallel sheet repeat). Short runs are
    demoted to coil. Pure CA heuristics - good enough to detect conserved
    vs changed secondary structure between homologs.
    """
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    ss = ['C'] * n
    if n < 6:
        return ss
    d3 = np.linalg.norm(coords[3:] - coords[:-3], axis=1)
    d5 = np.linalg.norm(coords[5:] - coords[:-5], axis=1)
    for i in range(n):
        if i < len(d3) and 4.5 <= d3[i] <= 5.8:
            ss[i] = 'H'
        elif i < len(d5) and 13.0 <= d5[i] <= 15.0:
            ss[i] = 'E'
    out = list(ss)
    for sym, minlen in (('H', 4), ('E', 3)):
        i = 0
        while i < n:
            if ss[i] == sym:
                j = i
                while j < n and ss[j] == sym:
                    j += 1
                if j - i < minlen:
                    for k in range(i, j):
                        out[k] = 'C'
                i = j
            else:
                i += 1
    # smooth: a single H inside a strand run etc. stays; keep simple.
    return out


def ss_segments(ss):
    """Yield (start, end, type) runs of H/E (half-open, min lengths applied)."""
    segs = []
    i = 0
    while i < len(ss):
        if ss[i] in 'HE':
            j = i
            while j < len(ss) and ss[j] == ss[i]:
                j += 1
            segs.append((i, j, ss[i]))
            i = j
        else:
            i += 1
    return segs


# ---------------------------------------------------------- sequence align

def sequence_pairs(seq_a, seq_b):
    """Local BLOSUM62 correspondence between two residue sequences.

    Returns a list of (i, j) residue index pairs covering the conserved
    core. Local alignment is essential: distant homologs carry unrelated
    terminal extensions that corrupt global alignments.
    """
    from Bio.Align import PairwiseAligner
    from Bio.Align import substitution_matrices
    aligner = PairwiseAligner()
    aligner.mode = "local"
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:
        pass
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    try:
        alignment = aligner.align(seq_a, seq_b)[0]
    except Exception:
        n = min(len(seq_a), len(seq_b))
        return list(zip(range(n), range(n)))
    pairs = []
    i_a, i_b = 0, 0
    for a_range, b_range in zip(alignment.aligned[0], alignment.aligned[1]):
        a0, a1 = a_range
        b0, b1 = b_range
        while i_a < a0:
            i_a += 1
        while i_b < b0:
            i_b += 1
        while i_a < a1 and i_b < b1:
            pairs.append((i_a, i_b))
            i_a += 1
            i_b += 1
    if len(pairs) < 50:
        n = min(len(seq_a), len(seq_b))
        pairs = list(zip(range(n), range(n)))
    return pairs


def core_identity(seq_a, seq_b, pairs):
    if not pairs:
        return 0.0
    same = sum(1 for i, j in pairs if seq_a[i] == seq_b[j])
    return same / len(pairs)


# ---------------------------------------------------------------- Kabsch

def kabsch(mobile, target):
    """Superpose mobile onto target; returns (rotation, translation).

    Convention: transformed = mobile @ rotation + translation.
    """
    mob = np.asarray(mobile, dtype=float)
    tgt = np.asarray(target, dtype=float)
    mob_c = mob - mob.mean(axis=0)
    tgt_c = tgt - tgt.mean(axis=0)
    cov = mob_c.T @ tgt_c
    u, _, vt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(u @ vt))
    corr = np.diag([1.0, 1.0, d])
    r = u @ corr @ vt
    shift = tgt.mean(axis=0) - mob.mean(axis=0) @ r
    return r, shift


def rmsd(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def superpose(mobile, reference, core_iterations=4):
    """Superpose mobile onto reference via sequence-matched CA pairs with
    iterative core refinement (conserved core drives the alignment)."""
    pairs = sequence_pairs(mobile["seq"], reference["seq"])
    if pairs is None:
        pairs = [(i, i) for i in range(min(len(mobile["coords"]),
                                           len(reference["coords"])))]
    m_idx = np.array([p[0] for p in pairs])
    r_idx = np.array([p[1] for p in pairs])
    m_idx = m_idx[m_idx < len(mobile["coords"])]
    r_idx = r_idx[r_idx < len(reference["coords"])]
    n = min(len(m_idx), len(r_idx))
    m_idx, r_idx = m_idx[:n], r_idx[:n]
    r, shift = kabsch(mobile["coords"][m_idx], reference["coords"][r_idx])
    keep = np.ones(len(m_idx), dtype=bool)
    for _ in range(core_iterations):
        dist = np.sqrt(np.sum(
            (mobile["coords"][m_idx] @ r + shift
             - reference["coords"][r_idx]) ** 2, axis=1))
        median = np.median(dist[keep])
        new_keep = dist <= max(2.0 * median, median + 4.0)
        if new_keep.sum() < 50 or (new_keep == keep).all():
            break
        keep = new_keep
        r, shift = kabsch(mobile["coords"][m_idx[keep]],
                          reference["coords"][r_idx[keep]])
    fit = rmsd(mobile["coords"][m_idx[keep]] @ r + shift,
               reference["coords"][r_idx[keep]])
    new_coords = mobile["coords"] @ r + shift
    return new_coords, fit, r, shift


# -------------------------------------------------------- stability mask

class StabilityMask:
    """Per-residue morph weights between two structures.

    For each corresponding residue pair (i in A, j in B) we compute:
      displacement  - CA distance between A and B after superposition
                      (both structures already share a common frame from
                      the aligner; this measures local divergence)
      ss_change     - True if the secondary-structure class differs

    weight = 0.0  -> frozen (stable + same SSE): residue does not move
    weight = 1.0  -> fully animated (large displacement or SSE change)

    Residues present in only one structure get fade weights handled by
    the caller (insertions fade in, deletions fade out).
    """

    def __init__(self, disp_max=4.0, disp_min=1.0):
        # displacement (A) below which a structurally-conserved residue is
        # fully frozen, and at/above which it is fully animated
        self.disp_max = disp_max
        self.disp_min = disp_min

    def weights(self, coords_a, coords_b, pairs, ss_a=None, ss_b=None):
        coords_a = np.asarray(coords_a, dtype=float)
        coords_b = np.asarray(coords_b, dtype=float)
        pairs = [(a, b) for a, b in pairs
                 if a < len(coords_a) and b < len(coords_b)]
        a_idx = np.array([p[0] for p in pairs], dtype=int)
        b_idx = np.array([p[1] for p in pairs], dtype=int)
        disp = np.linalg.norm(coords_a[a_idx] - coords_b[b_idx], axis=1)
        # conserved secondary structure + small displacement -> frozen;
        # ramps up to full animation at disp_max
        w_disp = np.clip((disp - self.disp_min)
                         / max(self.disp_max - self.disp_min, 1e-6),
                         0.0, 1.0)
        if ss_a is not None and ss_b is not None:
            ss_change = np.array([1.0 if ss_a[a] != ss_b[b] else 0.0
                                  for a, b in zip(a_idx, b_idx)])
        else:
            ss_change = np.zeros(len(pairs))
        # an SSE change always animates fully (new/lost helix or strand);
        # otherwise the displacement ramp decides
        w = np.maximum(w_disp, ss_change)
        # smoothstep for a softer transition edge
        w = w * w * (3.0 - 2.0 * w)
        return w, a_idx, b_idx


def ease(t):
    """Smoothstep easing for morph parameter t in [0, 1]."""
    t = min(max(t, 0.0), 1.0)
    return t * t * (3.0 - 2.0 * t)


# ------------------------------------------- reference-anchored morphing

def reference_pairs(struct_a, struct_b, reference):
    """Correspondence between two structures anchored through a common
    reference structure.

    Direct pairwise local alignment between two distant homologs can slip
    (especially over long terminal extensions), which makes the morph look
    like the whole protein reshuffles. Anchoring both structures through
    the alignment reference (the structure all coordinates were superposed
    onto) yields a consistent residue correspondence for the morph.

    struct_a/struct_b: dicts with "seq" (per-coordinate residue string)
    reference: dict with "seq" and "coords"
    Returns a list of (i, j) index pairs into a/b's coordinate arrays.
    """
    ref_seq = reference.get("seq") or "A" * len(reference["coords"])
    pa = sequence_pairs(struct_a.get("seq", ""), ref_seq)
    pb = sequence_pairs(struct_b.get("seq", ""), ref_seq)
    a_of_ref = {}
    for a, r in pa:
        if a < len(struct_a.get("coords", ())):
            a_of_ref[r] = a
    b_of_ref = {}
    for b, r in pb:
        if b < len(struct_b.get("coords", ())):
            b_of_ref[r] = b
    common = sorted(set(a_of_ref) & set(b_of_ref))
    return [(a_of_ref[r], b_of_ref[r]) for r in common]


def chain_orientation_ok(struct, reference):
    """Sanity check: the structure's coordinate chain must run in the same
    direction as the reference (N->C both, or C->N both).

    The check compares, in the superposed common frame, how far the
    structure's coordinate ends sit from the reference's N- and C-terminal
    regions. A correctly oriented chain has both ends near their
    corresponding reference ends; a reversed chain has the coordinate
    head near the reference C-terminus and vice versa. Returns True if the
    orientation matches (or the structure is too short to tell).
    """
    coords = np.asarray(struct.get("coords", ()), dtype=float)
    ref_coords = np.asarray(reference.get("coords", ()), dtype=float)
    if len(coords) < 24 or len(ref_coords) < 24:
        return True
    pairs = reference_pairs(struct, reference, reference)
    if not pairs:
        return True
    a_idx = np.array([p[0] for p in pairs])
    r_idx = np.array([p[1] for p in pairs])
    n = min(len(a_idx), len(r_idx))
    a_idx, r_idx = a_idx[:n], r_idx[:n]
    # superpose struct onto reference using the sequence correspondence
    r, shift = kabsch(coords[a_idx], ref_coords[r_idx])
    moved = coords @ r + shift
    k = min(12, len(ref_coords) // 4)
    head_ref = ref_coords[:k]
    tail_ref = ref_coords[-k:]
    # distance from each structure end to the NEAREST residue of each
    # reference end region (robust to local wiggle at the termini)
    def region_dist(pt, region):
        return float(np.min(np.linalg.norm(region - pt, axis=1)))
    d_head = region_dist(moved[0], head_ref)
    d_tail = region_dist(moved[-1], tail_ref)
    d_head_rev = region_dist(moved[0], tail_ref)
    d_tail_rev = region_dist(moved[-1], head_ref)
    return (d_head + d_tail) <= (d_head_rev + d_tail_rev)


# ------------------------------------------------------------ morph frame

def morph_frames(coords_a, coords_b, pairs, weights, n_frames,
                 fade_in=None):
    """Generate n_frames masked-morph frames (union space: handles both
    insertions and deletions).

    Returns frames in UNION space: paired residues + insertions (B-only)
    + deletions (A-only). Each frame is a dict with:
      "coords": (n_union, 3) array
      "radius_scale": (n_union,) array — 1.0 for paired, t for insertions
        (fade in), 1-t for deletions (fade out)
      "is_insertion": (n_union,) bool array
      "is_deletion": (n_union,) bool array

    Stable pairs (weight ~0) are anchored to A's position for the whole
    morph and only snap to B at the very end, so a conserved helix stays
    put. Insertions grow from the nearest paired neighbor's morphing
    position and fade in (radius_scale = t). Deletions collapse toward
    the nearest paired neighbor and fade out (radius_scale = 1-t).

    fade_in: optional boolean array over B marking insertions explicitly.
    """
    coords_a = np.asarray(coords_a, dtype=float)
    coords_b = np.asarray(coords_b, dtype=float)
    n_a = len(coords_a)
    n_b = len(coords_b)
    a_idx = np.array([p[0] for p in pairs], dtype=int)
    b_idx = np.array([p[1] for p in pairs], dtype=int)
    w = np.asarray(weights, dtype=float)
    paired_b = np.zeros(n_b, dtype=bool)
    paired_b[b_idx] = True
    paired_a = np.zeros(n_a, dtype=bool)
    paired_a[a_idx] = True
    frozen = w < 0.15
    ins = ~paired_b
    if fade_in is not None:
        ins = ins | np.asarray(fade_in, dtype=bool)
    dele = ~paired_a

    # Build union index mapping: walk through B, insert A-only residues
    # (deletions) at the right positions. Each union entry has:
    # (b_idx or -1, a_idx or -1, type)
    union_entries = []  # list of (b_idx, a_idx, 'paired'|'ins'|'del')
    b_to_a = {b: a for a, b in zip(a_idx, b_idx)}
    a_to_b = {a: b for a, b in zip(a_idx, b_idx)}

    j = 0  # B pointer
    while j < n_b:
        if paired_b[j]:
            union_entries.append((j, b_to_a[j], 'paired'))
            j += 1
        else:
            # Insertion block in B; check for deletions in A between
            # surrounding paired residues
            ins_start = j
            while j < n_b and not paired_b[j]:
                j += 1
            ins_end = j
            # Find A range for deletions: between the A indices of the
            # paired residues surrounding this gap
            a_before = b_to_a[ins_start - 1] if ins_start > 0 else -1
            a_after = b_to_a[ins_end] if ins_end < n_b else n_a
            # Add deletions (A-only residues in this gap)
            for ai in range(a_before + 1, a_after):
                if not paired_a[ai]:
                    union_entries.append((-1, ai, 'del'))
            # Add insertions
            for bi in range(ins_start, ins_end):
                union_entries.append((bi, -1, 'ins'))
    # Trailing deletions (A-only after last paired B)
    if len(a_idx) > 0:
        last_paired_a = max(a_idx)
        for ai in range(last_paired_a + 1, n_a):
            if not paired_a[ai]:
                union_entries.append((-1, ai, 'del'))

    n_union = len(union_entries)
    # Precompute weight lookup per union entry
    w_union = np.ones(n_union)
    bpos_union = np.full(n_union, -1, dtype=int)
    apos_union = np.full(n_union, -1, dtype=int)
    for u, (bi, ai, typ) in enumerate(union_entries):
        bpos_union[u] = bi
        apos_union[u] = ai
        if typ == 'paired':
            # find weight for this pair
            w_union[u] = w[np.where(b_idx == bi)[0][0]]

    frames = []
    for k in range(n_frames):
        t = ease(k / max(n_frames - 1, 1))
        out = np.zeros((n_union, 3))
        rs = np.ones(n_union)
        is_ins = np.zeros(n_union, dtype=bool)
        is_del = np.zeros(n_union, dtype=bool)

        for u, (bi, ai, typ) in enumerate(union_entries):
            if typ == 'paired':
                a = coords_a[ai]
                b = coords_b[bi]
                wu = w_union[u]
                if wu < 0.15:
                    settle = 1.0 if t >= 1.0 else 0.0
                    out[u] = a * (1 - settle) + b * settle
                else:
                    eff = 1.0 if t >= 1.0 else wu * t
                    out[u] = a * (1 - eff) + b * eff
                rs[u] = 1.0
            elif typ == 'ins':
                # Grow from nearest paired neighbor's morphing position
                anchor_u = u - 1  # preceding union entry
                while anchor_u >= 0 and union_entries[anchor_u][2] != 'paired':
                    anchor_u -= 1
                if anchor_u < 0:
                    # Search forward
                    anchor_u = u + 1
                    while anchor_u < n_union and union_entries[anchor_u][2] != 'paired':
                        anchor_u += 1
                if anchor_u >= 0 and anchor_u < n_union:
                    anchor_pos = out[anchor_u]
                else:
                    anchor_pos = coords_a.mean(axis=0)
                b = coords_b[bi]
                out[u] = anchor_pos * (1 - t) + b * t
                rs[u] = t
                is_ins[u] = True
            elif typ == 'del':
                # Collapse toward nearest paired neighbor
                anchor_u = u - 1
                while anchor_u >= 0 and union_entries[anchor_u][2] != 'paired':
                    anchor_u -= 1
                if anchor_u < 0:
                    anchor_u = u + 1
                    while anchor_u < n_union and union_entries[anchor_u][2] != 'paired':
                        anchor_u += 1
                if anchor_u >= 0 and anchor_u < n_union:
                    anchor_pos = out[anchor_u]
                else:
                    anchor_pos = coords_a[ai]
                a = coords_a[ai]
                out[u] = a * (1 - t) + anchor_pos * t
                rs[u] = 1.0 - t
                is_del[u] = True

        frames.append({
            "coords": out,
            "radius_scale": rs,
            "is_insertion": is_ins,
            "is_deletion": is_del,
        })

    # Exact endpoints: at t=0, paired at A, insertions at anchor, deletions at A
    # At t=1, paired at B, insertions at B, deletions at anchor (rs=0)
    return frames


# ---------------------------------------------------------- ribbon data

def ribbon_segments(coords, ss):
    """Split a CA trace into ribbon segments per SSE run.

    Returns a list of dicts: {type: 'H'|'E'|'C', idx: [i, j), coords}.
    Helices and strands become flat ribbon pieces; loops become a thin
    tube. Used by the matplotlib renderer to draw a cartoon-style ribbon.
    """
    coords = np.asarray(coords, dtype=float)
    segs = []
    for (i, j, typ) in ss_segments(ss):
        segs.append({"type": typ, "idx": (i, j),
                     "coords": coords[i:j]})
    # coil = everything not covered
    covered = np.zeros(len(coords), dtype=bool)
    for s in segs:
        covered[s["idx"][0]:s["idx"][1]] = True
    i = 0
    n = len(coords)
    while i < n:
        if not covered[i]:
            j = i
            while j < n and not covered[j]:
                j += 1
            segs.append({"type": 'C', "idx": (i, j), "coords": coords[i:j]})
            i = j
        else:
            i += 1
    segs.sort(key=lambda s: s["idx"][0])
    return segs


# ------------------------------------------------------------- colors

def plddt_colors(plddts, n=None):
    """RdYlGn colors from per-residue pLDDT, padded/truncated to n."""
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap("RdYlGn")
    p = np.asarray(plddts, dtype=float)
    if n is not None:
        if len(p) < n:
            p = np.concatenate([p, np.full(n - len(p), 70.0)])
        else:
            p = p[:n]
    if len(p) == 0:
        p = np.array([70.0])
    q = np.clip(p / 100.0, 0, 1)
    return cmap(0.15 + 0.85 * q)
