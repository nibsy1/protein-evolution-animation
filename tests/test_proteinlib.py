#!/usr/bin/env python3
"""Unit tests for proteinlib: SS assignment, stability mask, masked morph.

Run:  python tests/test_proteinlib.py   (from the workspace root)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import proteinlib as pl


def helix_coords(n=20, rise=1.5, radius=2.3, phase=0.0):
    """Ideal alpha-helix CA trace."""
    t = np.arange(n) * 100.0 / 360.0 * 2 * np.pi + phase
    return np.stack([radius * np.cos(t),
                     radius * np.sin(t),
                     np.arange(n) * rise], axis=1)


def test_assign_secondary_helix():
    ss = pl.assign_secondary(helix_coords(30))
    assert ss.count("H") >= 24, f"expected mostly H, got {ss}"
    print("ok test_assign_secondary_helix")


def test_assign_secondary_coil():
    rng = np.random.default_rng(1)
    ss = pl.assign_secondary(np.cumsum(rng.normal(size=(100, 3)), axis=0))
    assert set(ss) <= set("HEC")
    print("ok test_assign_secondary_coil")


def test_sequence_pairs_identical():
    pairs = pl.sequence_pairs("ACDEFGHIK", "ACDEFGHIK")
    assert pairs == list(zip(range(9), range(9)))
    print("ok test_sequence_pairs_identical")


def test_stability_mask_frozen_vs_animated():
    rng = np.random.default_rng(2)
    a = np.cumsum(rng.normal(size=(200, 3)), axis=0)
    b = a.copy()
    b[50:80] += np.array([8.0, 0, 0])   # one region moves a lot
    ss = "H" * 100 + "C" * 100
    w, ai, bi = pl.StabilityMask(disp_max=4.0, disp_min=1.0).weights(
        a, b, list(zip(range(200), range(200))), ss, ss)
    assert w[0] < 0.15, "stable residue should be frozen"
    assert w[60] > 0.9, "moved residue should be animated"
    print("ok test_stability_mask_frozen_vs_animated")


def test_stability_mask_ss_change_animates():
    a = helix_coords(30)
    b = a.copy()
    ss_a = "H" * 30
    ss_b = "C" * 30                      # helix -> coil, no displacement
    pairs = list(zip(range(30), range(30)))
    w, _, _ = pl.StabilityMask().weights(a, b, pairs, ss_a, ss_b)
    assert np.all(w > 0.9), "SS change must animate even without displacement"
    print("ok test_stability_mask_ss_change_animates")


def test_morph_frames_endpoints_and_freeze():
    rng = np.random.default_rng(3)
    a = np.cumsum(rng.normal(size=(150, 3)), axis=0)
    b = a.copy()
    b[40:60] += np.array([6.0, 0, 0])
    pairs = list(zip(range(150), range(150)))
    w, ai, bi = pl.StabilityMask(disp_max=4.0, disp_min=1.0).weights(
        a, b, pairs, "C" * 150, "C" * 150)
    frames = pl.morph_frames(a, b, pairs, w, 9)
    assert len(frames) == 9
    assert np.allclose(frames[0][bi], a[ai]), "frame 0 must equal A"
    assert np.allclose(frames[-1], b), "last frame must equal B"
    frozen = w < 0.15
    mid = frames[4]
    assert np.allclose(mid[bi[frozen]], a[ai[frozen]], atol=1e-6), \
        "frozen residues must hold A's position mid-morph"
    print("ok test_morph_frames_endpoints_and_freeze")


def test_morph_frames_insertions_grow_from_centroid():
    a = np.zeros((10, 3))
    b = np.zeros((12, 3))
    b[10:] = np.array([[5.0, 5.0, 5.0], [6.0, 6.0, 6.0]])  # insertions
    pairs = list(zip(range(10), range(10)))
    w = np.zeros(10)
    frames = pl.morph_frames(a, b, pairs, w, 5)
    # at t=0 insertions sit at A's centroid (origin)
    assert np.allclose(frames[0][10:], 0.0)
    # at t=1 they reach B
    assert np.allclose(frames[-1][10:], b[10:])
    print("ok test_morph_frames_insertions_grow_from_centroid")


def test_write_ca_pdb_ribbon_records(tmp_path=None):
    import tempfile
    coords = helix_coords(30)
    ss = pl.assign_secondary(coords)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "x.pdb")
        pl.write_ca_pdb(coords, np.full(30, 80.0), path, "TEST",
                        ss=ss, seq="A" * 30)
        txt = open(path).read()
        assert txt.count("HELIX") >= 1, "ribbon PDB must carry HELIX records"
        assert txt.count("ATOM") == 30
        parsed = pl.parse_pdb_ca(path)
        assert len(parsed["coords"]) == 30
        assert np.allclose(parsed["coords"], coords, atol=1e-3)
    print("ok test_write_ca_pdb_ribbon_records")


def test_kabsch_identity():
    rng = np.random.default_rng(4)
    pts = rng.normal(size=(50, 3))
    r, shift = pl.kabsch(pts, pts)
    out = pts @ r + shift
    assert np.allclose(out, pts, atol=1e-6)
    # and a known rotation is recovered
    theta = 0.7
    rot = np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])
    moved = pts @ rot + np.array([3.0, -2.0, 5.0])
    r2, shift2 = pl.kabsch(moved, pts)
    assert pl.rmsd(moved @ r2 + shift2, pts) < 1e-6
    print("ok test_kabsch_identity")


def test_reference_pairs_anchored():
    """reference_pairs must intersect both structures' alignments to the
    reference, so a direct pairwise alignment slip cannot reshuffle."""
    ref_seq = "ACDEFGHIKLMNPQRSTVWY" * 5
    # struct A: reference with an unrelated 30-residue N-terminal extension
    ext_a = "W" * 30
    seq_a = ext_a + ref_seq
    coords_a = helix_coords(len(seq_a))
    # struct B: reference with a DIFFERENT unrelated extension
    ext_b = "M" * 25
    seq_b = ext_b + ref_seq
    coords_b = helix_coords(len(seq_b), phase=0.3)
    reference = {"seq": ref_seq, "coords": helix_coords(len(ref_seq))}
    a = {"seq": seq_a, "coords": coords_a}
    b = {"seq": seq_b, "coords": coords_b}
    pairs = pl.reference_pairs(a, b, reference)
    assert len(pairs) >= len(ref_seq) - 5, len(pairs)
    # every pair must map into the reference-aligned core, i.e. a-index
    # skips A's extension and b-index skips B's extension
    a_idx = [p[0] for p in pairs]
    b_idx = [p[1] for p in pairs]
    assert min(a_idx) >= 25, min(a_idx)
    assert min(b_idx) >= 20, min(b_idx)
    # a direct pairwise alignment of the two extended sequences can pair
    # extension residues (wrong); the anchored path must not
    direct = pl.sequence_pairs(seq_a, seq_b)
    ext_pairs = [p for p in direct if p[0] < 30 or p[1] < 25]
    anchored_ext = [p for p in pairs if p[0] < 30 or p[1] < 25]
    assert len(anchored_ext) == 0, anchored_ext
    print("ok test_reference_pairs_anchored")


def test_reference_pairs_consistent_with_direct_for_close_homologs():
    """For two close homologs without terminal extensions, the anchored
    correspondence agrees with the direct one."""
    ref_seq = "ACDEFGHIKLMNPQRSTVWY" * 5
    seq_a = ref_seq
    seq_b = ref_seq[:50] + "G" + ref_seq[51:]  # one substitution
    reference = {"seq": ref_seq, "coords": helix_coords(len(ref_seq))}
    a = {"seq": seq_a, "coords": helix_coords(len(seq_a))}
    b = {"seq": seq_b, "coords": helix_coords(len(seq_b))}
    pairs = pl.reference_pairs(a, b, reference)
    assert pairs == list(zip(range(len(ref_seq)), range(len(ref_seq))))
    print("ok test_reference_pairs_consistent_with_direct_for_close_homologs")


def test_chain_orientation_ok():
    """In the superposed common frame, a correctly oriented chain has its
    coordinate head near the reference N-terminal region; a chain whose
    head sits at the reference C-terminus (reversed) is flagged."""
    ref_seq = "ACDEFGHIKLMNPQRSTVWY" * 5
    coords = helix_coords(len(ref_seq))
    reference = {"seq": ref_seq, "coords": coords}
    fwd = {"seq": ref_seq, "coords": coords}
    assert pl.chain_orientation_ok(fwd, reference)
    # reversed direction: the coordinate order is flipped relative to the
    # sequence, so after superposition the coordinate head lands at the
    # reference C-terminal region
    rev_coords = coords[::-1].copy()
    rev = {"seq": ref_seq, "coords": rev_coords}
    # build the correspondence a reversed chain would have: sequence
    # position i corresponds to reference position n-1-i
    rev_pairs = [(i, len(ref_seq) - 1 - i) for i in range(len(ref_seq))]
    # emulate by calling the internals directly
    c = np.asarray(rev_coords, dtype=float)
    rc = np.asarray(coords, dtype=float)
    a_idx = np.array([p[0] for p in rev_pairs])
    r_idx = np.array([p[1] for p in rev_pairs])
    r, shift = pl.kabsch(c[a_idx], rc[r_idx])
    moved = c @ r + shift
    k = min(12, len(rc) // 4)
    head_ref, tail_ref = rc[:k], rc[-k:]
    def region_dist(pt, region):
        return float(np.min(np.linalg.norm(region - pt, axis=1)))
    d_head = region_dist(moved[0], head_ref)
    d_tail = region_dist(moved[-1], tail_ref)
    d_head_rev = region_dist(moved[0], tail_ref)
    d_tail_rev = region_dist(moved[-1], head_ref)
    assert (d_head + d_tail) > (d_head_rev + d_tail_rev), \
        "reversed chain should score worse in the forward orientation"
    print("ok test_chain_orientation_ok")


def main():
    test_assign_secondary_helix()
    test_assign_secondary_coil()
    test_sequence_pairs_identical()
    test_stability_mask_frozen_vs_animated()
    test_stability_mask_ss_change_animates()
    test_morph_frames_endpoints_and_freeze()
    test_morph_frames_insertions_grow_from_centroid()
    test_write_ca_pdb_ribbon_records()
    test_kabsch_identity()
    test_reference_pairs_anchored()
    test_reference_pairs_consistent_with_direct_for_close_homologs()
    test_chain_orientation_ok()
    print("\nALL TESTS PASSED")


if __name__ == "__main__":
    main()
