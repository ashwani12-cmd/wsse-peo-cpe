#!/usr/bin/env python3
"""
==========================================================================
make_fig_msd.py -- regenerates fig1_msd_anisotropy.pdf

The original version of this figure annotated a fitted D_z and the ratio
D_xy/D_z. That is not a meaningful diffusivity for this cell: the WSSe
sheet spans the full periodic cross-section and the cell has ~12 A of
vacuum at each z boundary, so Li+ is confined to a finite slab in z and
MSD_z relaxes toward a plateau instead of growing linearly.

This version plots the same two MSD components but annotates only the
in-plane diffusivity, and marks the confinement limit MSD_z -> W^2/6 that
the out-of-plane component is approaching, so the figure states what the
data actually shows.

USAGE   python3 make_fig_msd.py
OUTPUT  fig1_msd_anisotropy.pdf
==========================================================================
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

import traj_io as T

FIT_LO, FIT_HI = 3.0, 10.0
OUT_PDF = 'fig1_msd_anisotropy.pdf'

plt.rcParams.update(T.PLOT_RC)
C_XY, C_Z = '#c0392b', '#2980b9'


def main():
    d = T.load()
    pos, typ = d['pos'], d['typ']
    sel_li = (typ == T.T_LI)
    sel_el = sel_li | np.isin(typ, T.T_TFSI) | np.isin(
        typ, (T.T_PEO_C, T.T_PEO_H, T.T_PEO_O))

    pos_u = T.unwrap_xy(pos)
    w = T.masses_of(typ).copy(); w[~sel_el] = 0.0
    pos_u = T.remove_com_drift(pos_u, w)
    li = pos_u[:, sel_li, :]

    nf = li.shape[0]; maxlag = nf // 2
    msd_xy = np.zeros(maxlag); msd_z = np.zeros(maxlag)
    for lag in range(1, maxlag):
        dd = li[lag:] - li[:-lag]
        msd_xy[lag] = np.mean((dd[:, :, :2] ** 2).sum(-1))
        msd_z[lag] = np.mean(dd[:, :, 2] ** 2)
    lag = np.arange(maxlag) * T.DT_FRAME_NS

    m = (lag > FIT_LO) & (lag < FIT_HI)
    res = stats.linregress(lag[m] * 1e-9, msd_xy[m] * 1e-16)
    D_xy = res.slope / 4.0

    # confinement limit: Li+ accessible slab width on each side of the sheet
    zc = pos[:, typ == T.T_W, 2].mean(axis=1)
    zrel = pos[:, sel_li, 2] - zc[:, None]
    side = zrel[0] > 0
    widths = []
    for mask in (side, ~side):
        v = zrel[:, mask].ravel()
        lo, hi = np.percentile(v, [1, 99])
        widths.append(abs(hi - lo))
    plateau = np.mean([wdt ** 2 / 6.0 for wdt in widths])

    beta_xy = np.gradient(np.log(msd_xy[1:]), np.log(lag[1:]))
    beta_z = np.gradient(np.log(msd_z[1:]), np.log(lag[1:]))
    wsel = (lag[1:] > FIT_LO) & (lag[1:] < FIT_HI)
    late = lag[1:] > 12.0

    print(f"D_xy (2D fit, MSD=4Dt) : {D_xy:.3e} cm^2/s")
    print(f"beta_xy {beta_xy[wsel].mean():.2f}  beta_z {beta_z[wsel].mean():.2f}"
          f"  (3-10 ns);  beta_z late (>12 ns) {beta_z[late].mean():.2f}")
    print(f"slab widths {widths[0]:.1f} / {widths[1]:.1f} A"
          f" -> confinement plateau W^2/6 ~ {plateau:.0f} A^2")
    print(f"MSD_z at {lag[-1]:.0f} ns = {msd_z[-1]:.0f} A^2"
          f" ({100*msd_z[-1]/plateau:.0f}% of the plateau)")

    fig, ax = plt.subplots(figsize=(7.4, 6))
    k = lag > 0
    ax.plot(lag[k], msd_xy[k], '-', color=C_XY, label='In-plane ($xy$)')
    ax.plot(lag[k], msd_z[k], '-', color=C_Z, label='Out-of-plane ($z$)')
    ax.plot(lag[m], 4 * D_xy * lag[m] * 1e-9 / 1e-16, 'k--', lw=1.8,
            label=r'$\mathrm{MSD}_{xy}=4D_{xy}t$ fit')
    ax.axhline(plateau, color=C_Z, ls=':', lw=2.0)
    ax.annotate('confinement limit $W^2/6$', xy=(lag[k].max() * 0.52, plateau),
                xytext=(lag[k].max() * 0.30, plateau * 1.35),
                color=C_Z, fontsize=13,
                arrowprops=dict(arrowstyle='->', color=C_Z, lw=1.5))
    ax.text(0.04, 0.95, rf'$D_{{xy}} = {D_xy*1e7:.2f}\times10^{{-7}}$ cm$^2$/s',
            transform=ax.transAxes, va='top', fontsize=14)
    ax.set_xlabel('Time lag [ns]')
    ax.set_ylabel(r'MSD [$\mathrm{\AA}^2$]')
    ax.set_xlim(0, lag[k].max())
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, loc='lower right', fontsize=13)
    T.apply_style(ax)
    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote {OUT_PDF}")


if __name__ == '__main__':
    main()
