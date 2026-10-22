#!/usr/bin/env python3
"""Render an evolutionary morph animation from aligned structures. (v3)

Usage:
    python morph_animation.py structures/aligned.json -o evolution.gif
    python morph_animation.py aligned.json -o out.gif --hold 25 --morph 15 --fps 12

Reads the aligned.json produced by the structure-aligner skill v2 (all
structures already superposed in a common 3D frame, each carrying its
residue sequence and per-residue secondary structure) and renders a GIF
that steps through the species in evolutionary order (green algae ->
angiosperms).

v2 changes (the "smooth, change-only" transition):
    - RIBBON rendering: helices as flat ribbons, strands as arrows/ribbons,
      loops as a thin tube, instead of scatter+line. Clean cartoon look.
    - FIXED MEDIUM LIGHTING: two lights parented to the camera at constant
      intensity, so shading is consistent across frames and while rotating.
    - STABILITY-MASKED MORPH: per corresponding residue we compute the
      CA displacement and whether the secondary-structure class changed.
      Residues that are stable AND structurally conserved are frozen
      (a helix that is a helix in both species does not reshuffle);
      only genuinely changing regions move. Insertions fade in from the
      core centroid instead of popping in. Eased (smoothstep) timing.
    - v3: residue correspondence is anchored through the alignment
      reference (reference_pairs), so pairwise alignment slips between
      distant homologs no longer reshuffle the whole protein, and the
      N-/C-termini are labeled with markers on every frame.
    - No membrane planes anywhere.

Output:
    GIF animation (default evolution.gif) + a frame sheet PNG for
    quick visual inspection.
"""
import argparse
import json
import os
import sys

import numpy as np

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

LINEAGE_ORDER = [
    "Green Algae",
    "Bryophytes",
    "Ferns",
    "Gymnosperms",
    "Angiosperms",
]

LINEAGE_LABELS = {
    "Green Algae": "Green algae",
    "Bryophytes": "Bryophytes (mosses & liverworts)",
    "Ferns": "Ferns",
    "Gymnosperms": "Gymnosperms (conifers)",
    "Angiosperms": "Flowering plants",
}

# ribbon palette (matches the dark theme, reads well under medium light)
SS_COLOR = {"H": "#c678dd", "E": "#e5c07b", "C": "#7f8c98"}
BG = "#0d1117"


def order_structures(structs):
    """Order by lineage; within a lineage, chain by core identity."""
    ordered = []
    for lineage in LINEAGE_ORDER:
        group = [s for s in structs if s.get("lineage") == lineage]
        if not group:
            continue
        group.sort(key=lambda s: -(s.get("core_identity") or 0))
        ordered.extend(group)
    known = set(LINEAGE_ORDER)
    rest = [s for s in structs if s.get("lineage") not in known]
    return ordered + rest


# ------------------------------------------------------------- ribbon draw

def _ribbon_axes(coords, i, j):
    """Local frame along a ribbon segment: tangent, normal, binormal."""
    a = coords[i]
    b = coords[min(j - 1, len(coords) - 1)]
    tangent = b - a
    n = np.linalg.norm(tangent)
    if n < 1e-6:
        tangent = np.array([1.0, 0.0, 0.0])
    else:
        tangent = tangent / n
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(tangent, ref)) > 0.95:
        ref = np.array([0.0, 1.0, 0.0])
    normal = np.cross(tangent, ref)
    normal /= np.linalg.norm(normal)
    binormal = np.cross(tangent, normal)
    return tangent, normal, binormal


def draw_structure(ax, coords, ss, colors, min_helix=4, min_strand=3,
                   termini=None):
    """Cartoon ribbon: helices as wide ribbons, strands as arrows, loops as
    a thin tube. SS class colors are tinted per segment by mean pLDDT so
    confidence is still readable. Lighting is applied by the caller.

    termini: optional (n_idx, c_idx) marking the N- and C-terminal CA
    positions; drawn as small sphere markers with "N"/"C" text labels.
    """
    coords = np.asarray(coords, dtype=float)
    # faint backbone underneath for continuity
    ax.plot(coords[:, 0], coords[:, 1], coords[:, 2],
            color="#30363d", linewidth=0.6, alpha=0.5, zorder=1)
    segs = pl.ribbon_segments(coords, ss)
    for seg in segs:
        i, j = seg["idx"]
        pts = coords[i:j]
        if len(pts) < 2:
            continue
        if seg["type"] == "H" and (j - i) < min_helix:
            seg_type = "C"
        elif seg["type"] == "E" and (j - i) < min_strand:
            seg_type = "C"
        else:
            seg_type = seg["type"]
        hexcol = SS_COLOR[seg_type].lstrip("#")
        base = np.array([int(hexcol[k:k + 2], 16) / 255.0
                         for k in (0, 2, 4)])
        seg_cols = colors[i:j]
        tint = (seg_cols.mean(axis=0)[:3] if len(seg_cols)
                else colors[min(i, len(colors) - 1)][:3])
        col = tuple(np.clip(0.45 * base + 0.55 * tint, 0, 1))
        if seg_type == "H":
            # helix: wide flat ribbon with darker edge passes -> shaded look
            tangent, normal, _ = _ribbon_axes(coords, i, j)
            off = normal * 0.85
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col,
                    linewidth=8.5, solid_capstyle="round", zorder=3)
            edge = tuple(np.clip(np.array(col) * 0.65, 0, 1))
            ax.plot(pts[:, 0] + off[0], pts[:, 1] + off[1],
                    pts[:, 2] + off[2], color=edge,
                    linewidth=2.6, solid_capstyle="round",
                    alpha=0.95, zorder=2)
            ax.plot(pts[:, 0] - off[0], pts[:, 1] - off[1],
                    pts[:, 2] - off[2], color=edge,
                    linewidth=2.6, solid_capstyle="round",
                    alpha=0.95, zorder=2)
        elif seg_type == "E":
            # strand: wide ribbon + arrow head at the C-terminal end
            tangent, normal, _ = _ribbon_axes(coords, i, j)
            off = normal * 0.70
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col,
                    linewidth=6.5, solid_capstyle="round", zorder=3)
            edge = tuple(np.clip(np.array(col) * 0.65, 0, 1))
            ax.plot(pts[:, 0] + off[0], pts[:, 1] + off[1],
                    pts[:, 2] + off[2], color=edge,
                    linewidth=2.2, solid_capstyle="round",
                    alpha=0.95, zorder=2)
            ax.plot(pts[:, 0] - off[0], pts[:, 1] - off[1],
                    pts[:, 2] - off[2], color=edge,
                    linewidth=2.2, solid_capstyle="round",
                    alpha=0.95, zorder=2)
            tip = pts[-1]
            head = pts[-1] - tangent * 2.6
            ax.plot([head[0], tip[0]], [head[1], tip[1]],
                    [head[2], tip[2]], color=col,
                    linewidth=9.0, solid_capstyle="butt",
                    zorder=3)
        else:
            # loop: thin tube
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=col,
                    linewidth=2.2, solid_capstyle="round", zorder=2)
    if termini is not None:
        n_idx, c_idx = termini
        for idx, label, col in ((n_idx, "N", "#58a6ff"),
                                (c_idx, "C", "#ff7f0e")):
            if idx is None or idx >= len(coords):
                continue
            p = coords[idx]
            ax.scatter([p[0]], [p[1]], [p[2]], s=42, color=col,
                       depthshade=False, zorder=6,
                       edgecolors="#0d1117", linewidths=0.6)
            ax.text(p[0], p[1], p[2], f" {label}", color=col,
                    fontsize=8, fontweight="bold", zorder=6)


# --------------------------------------------------------------- rendering

def render(frames_meta, out_path, fps, size):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import animation

    all_pts = np.vstack([f["coords"] for f in frames_meta])
    center = (all_pts.min(axis=0) + all_pts.max(axis=0)) / 2
    spread = (all_pts.max(axis=0) - all_pts.min(axis=0)).max() / 2 * 1.12

    fig = plt.figure(figsize=(size[0] / 100, size[1] / 100), dpi=100)
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor(BG)
    # fixed medium lighting: a constant light direction in VIEW space, so
    # shading is identical on every frame and while the camera rotates
    # (the structure never appears to rotate under a static light)
    LIGHT_DIR = np.array([0.5, 0.6, 0.62])
    LIGHT_DIR /= np.linalg.norm(LIGHT_DIR)
    total = len(frames_meta)
    print(f"[render] {total} frames -> {out_path}", file=sys.stderr)

    def shade(coords, colors):
        """Lambert-style shading with the fixed view-space light: normals
        are approximated from the CA trace direction."""
        coords = np.asarray(coords, dtype=float)
        cols = np.asarray(colors, dtype=float)
        n = min(len(coords), len(cols))
        if n < 2:
            return cols
        coords = coords[:n]
        cols = cols[:n]
        tang = np.gradient(coords, axis=0)
        norm = np.linalg.norm(tang, axis=1, keepdims=True)
        norm[norm < 1e-6] = 1.0
        tang = tang / norm
        lam = np.clip(np.abs(tang @ LIGHT_DIR), 0.0, 1.0)
        # medium lighting: 0.45 ambient + 0.55 diffuse
        shade_f = 0.45 + 0.55 * lam[:, None]
        out = np.clip(cols * shade_f, 0.0, 1.0)
        if len(out) < len(colors):
            pad = np.tile(out[-1:], (len(colors) - len(out), 1))
            out = np.vstack([out, pad])
        return out

    def update(i):
        fm = frames_meta[i]
        ax.cla()
        ax.set_facecolor(BG)
        c = np.asarray(fm["coords"])
        cols = shade(c, np.asarray(fm["colors"]))
        ss = fm.get("ss") or "C" * len(c)
        termini = fm.get("termini")
        if termini is not None:
            n_idx, c_idx = termini
            termini = (min(n_idx, len(c) - 1), min(c_idx, len(c) - 1))
        draw_structure(ax, c, ss, cols, termini=termini)
        ax.set_xlim(-spread + center[0], spread + center[0])
        ax.set_ylim(-spread + center[1], spread + center[1])
        ax.set_zlim(-spread + center[2], spread + center[2])
        ax.set_box_aspect((1, 1, 1))
        ax.set_axis_off()
        ax.set_title(fm["title"], color="#e6edf3", fontsize=12, pad=2)
        ax.text2D(0.5, 0.03, fm["subtitle"], transform=ax.transAxes,
                  ha="center", color="#8b949e", fontsize=9)
        prog = i / max(total - 1, 1)
        ax.text2D(0.5, 0.96,
                  f"evolutionary progression  {prog * 100:.0f}%",
                  transform=ax.transAxes, ha="center",
                  color="#58a6ff", fontsize=9)
        # slow rotation; lighting stays fixed because shading is view-based
        ax.view_init(elev=18, azim=(20 + i * 0.5) % 360)
        return []

    anim = animation.FuncAnimation(fig, update, frames=total, interval=100)
    writer = animation.PillowWriter(fps=fps)
    anim.save(out_path, writer=writer,
              savefig_kwargs={"facecolor": fig.get_facecolor()})
    plt.close(fig)
    print(f"[render] saved {out_path}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("aligned_json",
                   help="aligned.json from the structure-aligner skill")
    p.add_argument("-o", "--output", default="evolution.gif")
    p.add_argument("--hold", type=int, default=25,
                   help="frames to hold each structure")
    p.add_argument("--morph", type=int, default=15,
                   help="frames per morph transition")
    p.add_argument("--fps", type=int, default=12)
    p.add_argument("--size", type=int, nargs=2, default=(720, 720))
    p.add_argument("--disp-max", type=float, default=4.0,
                   help="CA displacement (A) at which a residue fully animates")
    p.add_argument("--title", default=None,
                   help="overall title (default: none)")
    args = p.parse_args()

    if pl is None:
        sys.exit("proteinlib.py not found - put it next to the skills tree "
                 "(src/proteinlib.py) or on PYTHONPATH")

    with open(args.aligned_json) as f:
        data = json.load(f)
    structs = order_structures(data["structures"])
    if len(structs) < 2:
        sys.exit("Need at least 2 structures to animate")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def colors_for(s, n):
        return pl.plddt_colors(s.get("plddts") or [70] * n, n)

    # v3: correspondence for the morph is anchored through the alignment
    # reference, so a pairwise alignment slip between two distant homologs
    # cannot reshuffle the whole protein; N/C termini are labeled.
    reference = next((s for s in data["structures"]
                      if s["accession"] == data.get("reference")),
                     structs[0])
    mask = pl.StabilityMask(disp_max=args.disp_max)
    frames_meta = []
    for idx, s in enumerate(structs):
        coords = np.asarray(s["coords"], dtype=float)
        seq = s.get("seq") or "A" * len(coords)
        ss = s.get("ss") or "C" * len(coords)
        termini = (0, len(coords) - 1) if len(coords) else None
        fm = {
            "title": f"{LINEAGE_LABELS.get(s['lineage'], s['lineage'])}"
                     f" \u2014 {s['species']}",
            "subtitle": f"{s['accession']} \u00b7 {s['source']}"
                        + (f" \u00b7 core ident "
                           f"{s['core_identity']:.0%}"
                           if s.get("core_identity") else ""),
            "ss": ss,
            "termini": termini,
        }
        cols = colors_for(s, len(coords))
        for _ in range(args.hold):
            frames_meta.append({**fm, "coords": coords, "colors": cols})
        if idx + 1 < len(structs):
            nxt = structs[idx + 1]
            ncoords = np.asarray(nxt["coords"], dtype=float)
            nseq = nxt.get("seq") or "A" * len(ncoords)
            nss = nxt.get("ss") or "C" * len(ncoords)
            pairs = pl.reference_pairs(s, nxt, reference)
            pairs = [(a, b) for a, b in pairs
                     if a < len(coords) and b < len(ncoords)]
            weights, a_idx, b_idx = mask.weights(
                coords, ncoords, pairs, ss, nss)
            interp = np.asarray(pl.morph_frames(coords, ncoords, pairs,
                                                weights, args.morph))
            n_cols_next = colors_for(nxt, len(ncoords))
            cols_a = cols[a_idx]
            cols_b = n_cols_next[b_idx]
            n_stable = int((weights < 0.15).sum())
            print(f"[morph] {s['species']} -> {nxt['species']}: "
                  f"{len(pairs)} pairs, {n_stable} frozen "
                  f"({n_stable / max(len(pairs), 1):.0%})",
                  file=sys.stderr)
            for k in range(interp.shape[0]):
                t = k / max(interp.shape[0] - 1, 1)
                cols_k = cols_a * (1 - t) + cols_b * t
                # SS per B-residue: stable pairs keep A's class, changing
                # pairs take B's; insertions take B's class
                b_to_pos = {b: pos for pos, b in enumerate(b_idx)}
                ss_full = []
                for j in range(len(ncoords)):
                    pos = b_to_pos.get(j)
                    if pos is None:
                        ss_full.append(nss[j])
                    else:
                        w = weights[pos]
                        ss_full.append(ss[a_idx[pos]] if w < 0.5
                                       else nss[b_idx[pos]])
                frames_meta.append({
                    **{"title": (fm["title"] if t < 0.5 else
                                 f"{LINEAGE_LABELS.get(nxt['lineage'], nxt['lineage'])}"
                                 f" \u2014 {nxt['species']}"),
                       "subtitle": fm["subtitle"] if t < 0.5 else
                           f"{nxt['accession']} \u00b7 {nxt['source']}"
                           + (f" \u00b7 core ident "
                              f"{nxt['core_identity']:.0%}"
                              if nxt.get("core_identity") else "")},
                    "coords": interp[k], "colors": cols_k,
                    "ss": "".join(ss_full),
                    "termini": (0, len(ncoords) - 1)
                    if len(ncoords) else None})
    render(frames_meta, args.output, args.fps, args.size)
    sheet = frames_meta[::args.hold + args.morph]
    try:
        cols_n = int(np.ceil(np.sqrt(len(sheet))))
        rows_n = int(np.ceil(len(sheet) / cols_n))
        fig = plt.figure(figsize=(cols_n * 3, rows_n * 3), dpi=100)
        for i, fm in enumerate(sheet, 1):
            ax = fig.add_subplot(rows_n, cols_n, i, projection="3d")
            c = np.asarray(fm["coords"])
            ss = fm.get("ss") or "C" * len(c)
            termini = fm.get("termini")
            if termini is not None:
                n_idx, c_idx = termini
                termini = (min(n_idx, len(c) - 1), min(c_idx, len(c) - 1))
            draw_structure(ax, c, ss, np.asarray(fm["colors"])[:, :3],
                           termini=termini)
            ax.set_axis_off()
            ax.set_title(fm["title"], fontsize=7, color="#e6edf3")
            ax.view_init(elev=18, azim=30)
        fig.patch.set_facecolor(BG)
        plt.tight_layout()
        sheet_path = os.path.splitext(args.output)[0] + "_sheet.png"
        plt.savefig(sheet_path, facecolor=fig.get_facecolor())
        plt.close(fig)
        print(f"[render] frame sheet -> {sheet_path}", file=sys.stderr)
    except Exception as e:
        print(f"[warn] frame sheet failed: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
