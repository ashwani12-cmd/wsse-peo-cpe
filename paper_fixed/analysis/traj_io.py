#!/usr/bin/env python3
"""
traj_io.py -- shared trajectory loader for the WSSe/PEO-LiTFSI analyses.

Parses prod.lammpstrj once and caches to traj_cache.npz, so the three
analysis scripts do not each re-read a 140 MB dump.

Dump format assumed:  id mol type x y z   (dump_modify sort id)

ATOM TYPES (from the data file masses and the groups in run_restart.in)
    1  C  (PEO)      2  H  (PEO)      3  O  (PEO ether)
    4  Li+
    5  N  (TFSI)     6  S  (TFSI)     7  O  (TFSI)
    8  C  (TFSI)     9  F  (TFSI)
   10  W            11  S  (WSSe)    12  Se (WSSe)

GEOMETRY NOTE -- READ BEFORE USING THESE COORDINATES
The cell is 40 x 40 x 82.6 A with the WSSe sheet spanning the FULL
periodic xy cross-section, and ~12 A of vacuum at each z boundary (the
electrolyte is a free-standing slab, not a space-filling melt). Two
consequences for every analysis built on this module:

  * Out-of-plane motion is confined, not diffusive. Diffusion
    coefficients must be taken from the IN-PLANE (xy) components only.
    get_msd_xy() enforces this.
  * The two electrolyte regions either side of the sheet do not exchange
    (31/31 -> 32/30 Li+ over 30 ns). Face-resolved quantities are
    averages over kinetically separate reservoirs.
"""

import os
import numpy as np

TRAJ = 'prod.lammpstrj'
CACHE = 'traj_cache.npz'

LX, LY, LZ = 40.0, 40.0, 82.5976
DT_FRAME_NS = 0.05                      # 50000 steps x 1 fs

MASS = {1: 12.011, 2: 1.008, 3: 15.999, 4: 6.941, 5: 14.007, 6: 32.065,
        7: 15.999, 8: 12.011, 9: 18.998, 10: 183.84, 11: 32.065, 12: 78.971}

T_PEO_C, T_PEO_H, T_PEO_O = 1, 2, 3
T_LI = 4
T_TFSI = (5, 6, 7, 8, 9)
T_W, T_S, T_SE = 10, 11, 12


def load():
    """Return dict with pos (nframes, natoms, 3), mol, typ, steps."""
    if os.path.exists(CACHE):
        z = np.load(CACHE)
        return dict(pos=z['pos'], mol=z['mol'], typ=z['typ'], steps=z['steps'])

    print(f"Parsing {TRAJ} (first run; result cached to {CACHE}) ...")
    steps, frames, mol, typ = [], [], None, None
    with open(TRAJ) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith('ITEM: TIMESTEP'):
                continue
            steps.append(int(f.readline()))
            f.readline()
            n = int(f.readline())
            f.readline()
            for _ in range(3):
                f.readline()
            f.readline()
            p = np.empty((n, 3)); m = np.empty(n, dtype=np.int32)
            t = np.empty(n, dtype=np.int32)
            for k in range(n):
                w = f.readline().split()
                m[k] = int(w[1]); t[k] = int(w[2])
                p[k] = (float(w[3]), float(w[4]), float(w[5]))
            frames.append(p)
            if mol is None:
                mol, typ = m, t
    pos = np.array(frames)
    steps = np.array(steps)
    print(f"  -> {pos.shape[0]} frames, {pos.shape[1]} atoms, "
          f"{steps[-1] - steps[0]} steps total")
    np.savez_compressed(CACHE, pos=pos, mol=mol, typ=typ, steps=steps)
    return dict(pos=pos, mol=mol, typ=typ, steps=steps)


def masses_of(typ):
    return np.array([MASS[t] for t in typ])


def unwrap_xy(pos):
    """Unwrap x and y across periodic images (z needs no unwrapping: the
    slab is bounded by vacuum, so nothing crosses the z boundary)."""
    out = pos.copy()
    box = np.array([LX, LY])
    for fi in range(1, out.shape[0]):
        d = out[fi, :, :2] - out[fi - 1, :, :2]
        d -= box * np.round(d / box)
        out[fi, :, :2] = out[fi - 1, :, :2] + d
    return out


def molecule_com(pos, mol, typ, sel):
    """Mass-weighted COM per molecule for the atoms selected by `sel`.
    Returns (nframes, nmol, 3) and the molecule ids."""
    mids = np.unique(mol[sel])
    m = masses_of(typ)
    com = np.zeros((pos.shape[0], len(mids), 3))
    for j, mid in enumerate(mids):
        idx = np.where(sel & (mol == mid))[0]
        w = m[idx]
        com[:, j, :] = (pos[:, idx, :] * w[None, :, None]).sum(1) / w.sum()
    return com, mids


def remove_com_drift(pos, weights):
    """Subtract the mass-weighted COM of a reference set at each frame.

    Use the whole electrolyte (PEO + Li + TFSI) as the reference, NOT the
    Li group alone: removing the cation-group COM would also remove the
    collective cation motion that carries charge.
    """
    w = weights / weights.sum()
    com = (pos * w[None, :, None]).sum(axis=1, keepdims=True)
    return pos - com


def msd_sliding_xy(pos_unwrapped):
    """Sliding-origin in-plane MSD. Returns (lag_ns, msd_xy)."""
    nf = pos_unwrapped.shape[0]
    maxlag = nf // 2
    msd = np.zeros(maxlag)
    for lag in range(1, maxlag):
        d = pos_unwrapped[lag:, :, :2] - pos_unwrapped[:-lag, :, :2]
        msd[lag] = np.mean((d ** 2).sum(-1))
    return np.arange(maxlag) * DT_FRAME_NS, msd


def autocorr_time(x):
    """Integrated autocorrelation time (ns) of a per-frame series."""
    y = np.asarray(x, float) - np.mean(x)
    if y.std() == 0:
        return DT_FRAME_NS
    ac = np.correlate(y, y, 'full')[len(y) - 1:]
    ac = ac / ac[0]
    below = np.where(ac < 0.05)[0]
    cut = below[0] if len(below) else min(len(ac), 50)
    return (1 + 2 * ac[1:cut].sum()) * DT_FRAME_NS


def block_sem(x):
    """SEM from blocks of >=5 autocorrelation times, >=5 blocks.
    Returns (sem, tau_ns, n_eff)."""
    x = np.asarray(x, float)
    tau = autocorr_time(x)
    nb = int(len(x) * DT_FRAME_NS / max(5.0 * tau, 2.0))
    nb = max(5, min(nb, len(x) // 4))
    means = np.array([b.mean() for b in np.array_split(x, nb)])
    return means.std(ddof=1) / np.sqrt(nb), tau, len(x) * DT_FRAME_NS / tau


PLOT_RC = {
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 16, "axes.labelsize": 16, "axes.titlesize": 16,
    "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 13,
    "axes.linewidth": 2.0, "lines.linewidth": 2.2, "lines.markersize": 7,
}


def apply_style(ax):
    ax.tick_params(axis='both', which='major', direction='inout', length=7, width=2)
    ax.tick_params(axis='both', which='minor', direction='in', length=4, width=1.4)
    ax.minorticks_on()
