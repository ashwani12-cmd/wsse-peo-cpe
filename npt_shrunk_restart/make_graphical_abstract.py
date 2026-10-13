#!/usr/bin/env python3
"""
==========================================================================
make_graphical_abstract.py -- Elsevier graphical abstract

Single wide panel telling the paper's story in three steps:
  (a) the system             -- Janus WSSe monolayer in PEO/LiTFSI
  (b) cross-plane transport  -- blocked by the fully-spanning sheet
  (c) in-plane transport     -- obstructed vs. the filler-free control

Panel (c) previously showed the Se-face Li+ enrichment ("1.5x,
3.4 sigma"). That result did not survive the null/causal/replica
control battery (fig8_control_battery.pdf) and main.tex now reports it
as a tested-and-falsified hypothesis, not a finding. Panel (c) is
replaced with the result that IS attributable to the filler: the
24%/18% drop in in-plane D(Li+)/D(TFSI-) against the filler-free
control, taken from ../control_nofiller/control_comparison.txt.

Elsevier asks for a minimum of 531 x 1328 px (h x w). This renders at
1600 x 640 px, comfortably above that, and also writes a PDF.

USAGE   python3 make_graphical_abstract.py
OUTPUT  graphical_abstract.png, graphical_abstract.pdf
==========================================================================
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import traj_io as T

SNAPSHOT = '../paper_fixed/figures/PEO-system.png'
OUT = 'graphical_abstract'

# composite vs filler-free control, in-plane diffusivities
# [1e-7 cm^2/s], from ../control_nofiller/control_comparison.txt
D_LI = {'composite': 4.40, 'control': 5.79}      # -24.0%
D_TFSI = {'composite': 11.67, 'control': 14.19}  # -17.8%

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 15, "axes.labelsize": 15,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 12,
    "axes.linewidth": 1.8, "lines.linewidth": 2.4,
})
C_S, C_SE, C_XY, C_Z = '#b8860b', '#1f618d', '#c0392b', '#2980b9'


def style(ax):
    ax.tick_params(axis='both', which='major', direction='in', length=5, width=1.6)
    ax.minorticks_on()
    ax.tick_params(axis='both', which='minor', direction='in', length=3, width=1.1)


def main():
    d = T.load()
    pos, typ = d['pos'], d['typ']
    sel_li = typ == T.T_LI
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

    zc = pos[:, typ == T.T_W, 2].mean(axis=1)
    zrel = pos[:, sel_li, 2] - zc[:, None]
    side = zrel[0] > 0
    widths = [abs(np.subtract(*np.percentile(zrel[:, m].ravel(), [99, 1])))
              for m in (side, ~side)]
    plateau = np.mean([x ** 2 / 6 for x in widths])

    fig = plt.figure(figsize=(16, 6.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1, 1], wspace=0.28,
                          left=0.015, right=0.985, top=0.90, bottom=0.14)

    # (a) system snapshot
    ax0 = fig.add_subplot(gs[0])
    try:
        ax0.imshow(mpimg.imread(SNAPSHOT))
    except Exception as e:
        ax0.text(0.5, 0.5, 'snapshot\nnot found', ha='center', va='center',
                 transform=ax0.transAxes)
        print(f"  (snapshot missing: {e})")
    ax0.axis('off')
    ax0.set_title('Janus WSSe in PEO/LiTFSI', fontsize=16, pad=8)

    # (b) cross-plane blocking
    ax1 = fig.add_subplot(gs[1])
    k = lag > 0
    ax1.plot(lag[k], msd_xy[k], '-', color=C_XY, label='in-plane')
    ax1.plot(lag[k], msd_z[k], '-', color=C_Z, label='cross-plane')
    ax1.axhline(plateau, color=C_Z, ls=':', lw=2)
    ax1.text(lag[k].max() * 0.46, plateau * 1.08, 'confinement limit',
             color=C_Z, fontsize=12)
    ax1.set_xlabel('time lag [ns]')
    ax1.set_ylabel(r'Li$^+$ MSD [$\mathrm{\AA}^2$]')
    ax1.set_xlim(0, lag[k].max()); ax1.set_ylim(bottom=0)
    ax1.legend(frameon=False, loc='upper left')
    ax1.set_title(r'Cross-plane transport blocked', fontsize=16, pad=8)
    ax1.text(0.97, 0.06, r'$\sim$1 crossing in 30 ns', transform=ax1.transAxes,
             ha='right', fontsize=13, color=C_Z)
    style(ax1)

    # (c) in-plane transport obstructed by the filler (vs. filler-free
    # control) -- the result that survives, replacing the retracted
    # Se-face-enrichment panel
    ax2 = fig.add_subplot(gs[2])
    species = [r'Li$^+$', r'TFSI$^-$']
    x = np.arange(len(species))
    width = 0.32
    composite_vals = [D_LI['composite'], D_TFSI['composite']]
    control_vals = [D_LI['control'], D_TFSI['control']]
    pct = [100 * (c - f) / c for f, c in zip(composite_vals, control_vals)]

    ax2.bar(x - width / 2, control_vals, width, color='#95a5a6',
            label='filler-free control')
    ax2.bar(x + width / 2, composite_vals, width, color=C_XY,
            label='composite (with sheet)')
    for xi, (f, p) in enumerate(zip(composite_vals, pct)):
        ax2.text(xi + width / 2, f + 0.35, f'-{p:.0f}%', ha='center',
                 va='bottom', fontsize=13, color=C_XY, fontweight='bold')
    ax2.set_xticks(x); ax2.set_xticklabels(species, fontsize=14)
    ax2.set_ylabel(r'in-plane $D$ [$10^{-7}\,\mathrm{cm^2/s}$]')
    ax2.set_ylim(0, max(control_vals) * 1.22)
    ax2.legend(frameon=False, loc='upper right', fontsize=11)
    ax2.set_title('Filler obstructs in-plane transport', fontsize=16, pad=8)
    style(ax2)

    for f, dpi in ((OUT + '.png', 100), (OUT + '.pdf', 100)):
        fig.savefig(f, dpi=dpi, bbox_inches=None)
        print(f"Wrote {f}")
    plt.close(fig)
    im = mpimg.imread(OUT + '.png')
    print(f"  size {im.shape[1]} x {im.shape[0]} px "
          f"(Elsevier minimum 1328 x 531)")


if __name__ == '__main__':
    main()
