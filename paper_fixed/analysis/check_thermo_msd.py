#!/usr/bin/env python3
"""
check_thermo_msd.py -- demonstrate why the live thermo MSD stream must NOT be
used to judge convergence (NOTES.md section 13).

`compute msd` writes a SINGLE-ORIGIN MSD over only 62 Li+. At late lag it is
noise: block exponents swing wildly and 3 ns slopes can come out negative.
The paper's own analysis is fine -- analyze_transport.py uses the sliding-origin
estimator traj_io.msd_sliding_xy -- this script exists so nobody is tempted to
read convergence off out_*.dat while a run is in flight.

Usage:
    python3 check_thermo_msd.py <msd_file> [more...] [--t0-step N]

Default --t0-step is 10000000: production began at step 10 M and the compute's
reference coordinates are carried in the restart file, so t = (step - 1e7) fs.
"""
import sys

import numpy as np


def load(path):
    s, m = [], []
    for ln in open(path):
        if ln.startswith('#'):
            continue
        p = ln.split()
        if len(p) < 2:
            continue
        try:
            s.append(float(p[0])); m.append(float(p[1]))
        except ValueError:
            pass
    return np.array(s), np.array(m)


def main(paths, t0_step):
    s = np.concatenate([load(p)[0] for p in paths])
    m = np.concatenate([load(p)[1] for p in paths])
    o = np.argsort(s); s, m = s[o], m[o]
    _, u = np.unique(s, return_index=True); s, m = s[u], m[u]
    t = (s - t0_step) * 1e-6                      # fs -> ns
    print(f"continuous MSD: t = {t[0]:.4f} to {t[-1]:.3f} ns, {len(t)} points")
    print(f"MSD: {m[0]:.3f} -> {m[-1]:.1f} A^2")

    print("\n2 ns blocks -- log-log exponent beta (1.0 = diffusive):")
    bad = []
    for a in range(2, int(t[-1]) - 1, 2):
        k = (t >= a) & (t < a + 2)
        if k.sum() < 100:
            continue
        b = np.polyfit(np.log(t[k]), np.log(m[k]), 1)[0]
        bad.append(b)
        print(f"  {a:3d}-{a+2:3d} ns : beta={b:6.2f}")
    if bad:
        print(f"\n  beta range over blocks: {min(bad):.2f} to {max(bad):.2f}")
        print("  -> a single-origin MSD over ~62 ions. If this range is wide, any")
        print("     'beta ~ 1' read off a late window is noise, NOT the diffusive")
        print("     regime. Use traj_io.msd_sliding_xy on the dump instead.")

    print("\n3 ns sub-blocks -- apparent 3D D (MSD = 6Dt):")
    for a in range(int(t[-1]) - 9, int(t[-1]) - 2, 3):
        k = (t >= a) & (t < a + 3)
        if k.sum() < 100:
            continue
        sl = np.polyfit(t[k], m[k], 1)[0]
        flag = "   <-- NEGATIVE: noise, not physics" if sl < 0 else ""
        print(f"  {a:3d}-{a+3:3d} ns : slope={sl:8.2f} A^2/ns  "
              f"D={sl / 6 * 1e-7:.3e} cm^2/s{flag}")


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    t0 = 10_000_000
    if '--t0-step' in sys.argv:
        t0 = int(sys.argv[sys.argv.index('--t0-step') + 1])
        args = [a for a in args if a != str(t0)]
    if not args:
        sys.exit(__doc__)
    main(args, t0)
