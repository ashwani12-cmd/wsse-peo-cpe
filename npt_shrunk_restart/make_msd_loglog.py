#!/usr/bin/env python3
"""
==========================================================================
make_msd_loglog.py

Produces the log-log MSD diagnostic figure for the WSSe/PEO/LiTFSI CPE paper:
  - In-plane and out-of-plane Li+ MSD on log-log axes
  - Reference slope=1 line (ideal Einstein diffusion)
  - Prints the diffusion exponent beta = d(log MSD)/d(log t)

This plot reveals whether the MSD has reached the true diffusive regime
(slope -> 1) or is still sub-diffusive (slope < 1), which is the standard
diagnostic reviewers expect for polymer-electrolyte diffusivities.

USAGE:
    python3 make_msd_loglog.py

REQUIRED INPUT:
    prod.lammpstrj   - trajectory dump (id mol type x y z)

OUTPUT:
    fig1b_msd_loglog.pdf
==========================================================================
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ----------------------------- CONFIG -----------------------------
TRAJ_FILE = 'prod.lammpstrj'
OUT_PDF   = 'fig1b_msd_loglog.pdf'

LI_TYPE = 4                       # atom type for Li+
LX, LY, LZ = 40.0, 40.0, 82.5976  # box dims (Angstrom) for PBC unwrap
FIT_LO, FIT_HI = 3.0, 10.0        # window (ns) for reporting mean exponent
# ------------------------------------------------------------------

# ============================================================
# JOURNAL-QUALITY GLOBAL SETTINGS
# ============================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 18,
    "axes.labelsize": 18,
    "axes.titlesize": 18,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
    "legend.fontsize": 15,
    "axes.linewidth": 2.0,
    "lines.linewidth": 2.5,
    "lines.markersize": 8,
})


def apply_style(ax):
    ax.tick_params(axis='both', which='major',
                   direction='inout', length=8, width=2)
    ax.tick_params(axis='both', which='minor',
                   direction='in', length=5, width=1.5)
    ax.minorticks_on()


def read_li_positions(fname):
    """Parse LAMMPS dump; return list of (timestep, Li xyz array)."""
    print(f"Reading {fname} ...")
    with open(fname) as f:
        lines = f.readlines()
    timesteps, li_frames = [], []
    i = 0
    while i < len(lines):
        if 'ITEM: TIMESTEP' in lines[i]:
            ts = int(lines[i + 1])
            natoms = int(lines[i + 3])
            li = []
            for k in range(natoms):
                p = lines[i + 9 + k].split()
                if int(p[2]) == LI_TYPE:
                    li.append([float(p[3]), float(p[4]), float(p[5])])
            timesteps.append(ts)
            li_frames.append(np.array(li))
            i += 9 + natoms
        else:
            i += 1
    print(f"  -> {len(li_frames)} frames, {len(li_frames[0])} Li+ ions")
    return np.array(timesteps), li_frames


def sliding_window_msd(timesteps, li_frames):
    """Compute directional sliding-window MSD with PBC unwrap + COM removal."""
    print("Computing sliding-window MSD ...")
    box = np.array([LX, LY, LZ])
    nf = len(li_frames)
    ni = len(li_frames[0])

    pos = np.zeros((nf, ni, 3))
    for fi in range(nf):
        pos[fi] = li_frames[fi]

    # Unwrap PBC
    for j in range(ni):
        for fi in range(1, nf):
            d = pos[fi, j] - pos[fi - 1, j]
            d -= box * np.round(d / box)
            pos[fi, j] = pos[fi - 1, j] + d

    # Remove Li-group center of mass drift
    pos -= pos.mean(axis=1, keepdims=True)

    max_lag = nf // 2
    msd_xy = np.zeros(max_lag)
    msd_z = np.zeros(max_lag)
    for lag in range(1, max_lag):
        dxy, dz = [], []
        for t0 in range(nf - lag):
            d = pos[t0 + lag] - pos[t0]
            dxy.append(np.mean(d[:, 0]**2 + d[:, 1]**2))
            dz.append(np.mean(d[:, 2]**2))
        msd_xy[lag] = np.mean(dxy)
        msd_z[lag] = np.mean(dz)

    lag_ns = (timesteps[1] - timesteps[0]) / 1e6   # 1 fs timestep -> ns
    lag_t = np.arange(max_lag) * lag_ns
    return lag_t, msd_xy, msd_z


def report_exponents(lag_t, msd_xy, msd_z):
    """Print local diffusion exponent beta = d log MSD / d log t."""
    m = lag_t > 0
    lt = lag_t[m]
    mxy = msd_xy[m]
    mz = msd_z[m]
    logt = np.log(lt)
    beta_xy = np.gradient(np.log(mxy), logt)
    beta_z = np.gradient(np.log(mz), logt)

    print("\nDiffusion exponent beta = d(log MSD)/d(log t)  [want 1.0]:")
    print(f"{'lag(ns)':>8} {'beta_xy':>9} {'beta_z':>9}")
    for target in [0.5, 1, 2, 3, 5, 8, 12]:
        idx = np.argmin(np.abs(lt - target))
        print(f"{lt[idx]:>8.1f} {beta_xy[idx]:>9.3f} {beta_z[idx]:>9.3f}")

    fit = (lt > FIT_LO) & (lt < FIT_HI)
    print(f"\nMean beta in {FIT_LO}-{FIT_HI} ns window:")
    print(f"  beta_xy = {beta_xy[fit].mean():.3f}")
    print(f"  beta_z  = {beta_z[fit].mean():.3f}")
    if beta_xy[fit].mean() < 0.9:
        print("  -> sub-diffusive (beta < 1): report D as APPARENT diffusivity.")
    return beta_xy[fit].mean(), beta_z[fit].mean()


def make_loglog_figure(lag_t, msd_xy, msd_z):
    print("Building log-log figure ...")
    m = lag_t > 0
    lt = lag_t[m]
    mxy = msd_xy[m]
    mz = msd_z[m]

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(lt, mxy, '-', color='#c0392b', label='In-plane ($xy$)')
    ax.loglog(lt, mz, '-', color='#2980b9', label='Out-of-plane ($z$)')

    # Reference slope = 1 line (ideal Einstein diffusion), anchored at t=1 ns
    anchor_idx = np.argmin(np.abs(lt - 1.0))
    ref_t = np.array([1.0, lt.max()])
    ref_y = mxy[anchor_idx] * (ref_t / 1.0) ** 1.0
    ax.loglog(ref_t, ref_y, 'k--', lw=2, label='slope = 1 (Einstein)')

    ax.set_xlabel('Time lag [ns]')
    ax.set_ylabel(r'MSD [$\mathrm{\AA}^2$]')
    ax.legend(frameon=False, loc='upper left')
    apply_style(ax)
    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {OUT_PDF}")


if __name__ == '__main__':
    timesteps, li_frames = read_li_positions(TRAJ_FILE)
    lag_t, msd_xy, msd_z = sliding_window_msd(timesteps, li_frames)
    report_exponents(lag_t, msd_xy, msd_z)
    make_loglog_figure(lag_t, msd_xy, msd_z)
    print("Done.")
