#!/usr/bin/env python3
"""
analyze_reservoirs.py -- are the two electrolyte reservoirs exchanging, and
are they compositionally symmetric?

This is the check behind NOTES.md section 13: the 150 ns extension cannot fix
the face-selectivity systematic, because the sheet spans the full periodic
cross-section and Li+ essentially never change sides. Each face is therefore
sampled by a disjoint ion set in a single realization of the packing -- which
is why replicas (rep1/rep2) exist at all.

Re-run this on each new trajectory as it lands. A replica is only useful if
its reservoirs are symmetric to the same ~few percent as the parent run;
if a replica comes out badly asymmetric, say so rather than averaging it in.

Usage:
    python3 analyze_reservoirs.py <dir-with-prod.lammpstrj> [...]

<dir> must contain prod.lammpstrj (or an existing traj_cache.npz); traj_io.py
caches the parse, so the second call on a directory is fast.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traj_io as T


def report(run_dir):
    cwd = os.getcwd()
    os.chdir(run_dir)
    try:
        d = T.load()
    finally:
        os.chdir(cwd)

    pos, typ = d['pos'], d['typ']
    z = pos[:, :, 2]
    zw = z[:, typ == T.T_W].mean()
    zs, zse = z[:, typ == T.T_S].mean(), z[:, typ == T.T_SE].mean()
    s_side = '+' if zs > zw else '-'

    print(f"\n{'='*66}\n{run_dir}\n{'='*66}")
    print(f"  W midplane z={zw:.2f}  S-face z={zs:.2f}  Se-face z={zse:.2f} "
          f"(S face on the {s_side}z side)")

    sel_li = typ == T.T_LI
    zl = z[:, sel_li] - zw
    zl = (zl + T.LZ / 2) % T.LZ - T.LZ / 2     # signed, periodic
    side = np.sign(zl)
    nf, nion = side.shape
    print(f"  {nion} Li+ over {nf} frames ({nf * T.DT_FRAME_NS:.1f} ns)")

    # --- exchange across the sheet -------------------------------------
    pop = (side > 0).sum(1)
    flips = (side[1:] != side[:-1]).sum(0)
    ever_p, ever_m = (side > 0).any(0), (side < 0).any(0)
    print(f"  population on +z side: first {pop[0]}, last {pop[-1]}, "
          f"range {pop.min()}-{pop.max()}")
    print(f"  crossing events: {flips.sum()} from {(flips > 0).sum()} distinct ions; "
          f"{nion - (ever_p & ever_m).sum()} of {nion} never changed sides")
    if flips.sum():
        print(f"    -> mean time between crossings per ion "
              f"~ {nion * nf * T.DT_FRAME_NS / flips.sum():.0f} ns")
    if (flips > 0).sum() and flips.sum() / max((flips > 0).sum(), 1) > 3:
        print("    NOTE: events concentrated in very few ions -- that is "
              "flickering at the\n          dividing surface, not transport. "
              "Treat exchange as zero.")

    # --- compositional symmetry of the two reservoirs -------------------
    print("  reservoir composition (frame-averaged):")
    sl, ol = ('S', 'Se') if s_side == '+' else ('Se', 'S')
    worst = 0.0
    for lab, mask in [("Li", sel_li), ("TFSI N", typ == 5),
                      ("PEO O", typ == T.T_PEO_O), ("PEO C", typ == T.T_PEO_C)]:
        zz = z[:, mask] - zw
        zz = (zz + T.LZ / 2) % T.LZ - T.LZ / 2
        a, b = (zz > 0).sum(1).mean(), (zz < 0).sum(1).mean()
        print(f"    {lab:7s} {sl}-side {a:8.1f}   {ol}-side {b:8.1f}   "
              f"ratio {a / b:.3f}")
        worst = max(worst, abs(a / b - 1.0))
    print(f"  worst reservoir asymmetry: {100 * worst:.1f}%"
          f"{'  <-- too large to average in blindly' if worst > 0.10 else ''}")


if __name__ == '__main__':
    dirs = sys.argv[1:] or ['.']
    for dd in dirs:
        if not (os.path.exists(os.path.join(dd, T.TRAJ))
                or os.path.exists(os.path.join(dd, T.CACHE))):
            print(f"\n{dd}: no {T.TRAJ} or {T.CACHE} yet -- skipping")
            continue
        report(dd)
