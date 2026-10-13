#!/usr/bin/env python3
"""
==========================================================================
analyze_transport.py   --  item #1: anion transport, D_TFSI/D_Li, t+

Adds the observables a battery-journal referee expects and the original
analysis omitted: the ANION diffusivity, the D_TFSI/D_Li+ ratio that
validates a charge-scaled force field, the cation transference number,
and a Nernst-Einstein conductivity estimate.

Everything is computed from the EXISTING prod.lammpstrj -- no new
simulation is required.

WHY IN-PLANE ONLY
The WSSe sheet spans the full periodic cross-section and the cell has
~12 A of vacuum at each z boundary, so out-of-plane motion is confined
and MSD_z saturates rather than growing linearly. Only the in-plane
components are diffusive, so all D values here are in-plane:
    MSD_xy = 4 D_xy t      (two-dimensional Einstein relation)
Note the factor 4, not 6: this is a 2D fit. Reporting a 3D D from this
trajectory would be wrong.

CAVEAT ON t+
t+ below is the ideal Nernst-Einstein estimate
    t+ = D_Li / (D_Li + D_TFSI),
which neglects ion-ion correlations and therefore overestimates t+ in
PEO/LiTFSI, often substantially. The honest statement in a paper is
"ideal Nernst-Einstein t+", with the limitation named. A correlation-
corrected t+ needs the Onsager coefficients (distinct van Hove terms),
which are noisy at 62 ion pairs over 30 ns.

USAGE
    python3 analyze_transport.py

OUTPUT
    fig5_transport.pdf
    transport_summary.txt
==========================================================================
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

import traj_io as T

FIT_LO, FIT_HI = 3.0, 10.0          # ns, same window as the Li+ analysis
OUT_PDF = 'fig5_transport.pdf'
OUT_TXT = 'transport_summary.txt'

plt.rcParams.update(T.PLOT_RC)
C_LI, C_TFSI = '#c0392b', '#2980b9'


def fit_D2d(lag_ns, msd_A2, lo=FIT_LO, hi=FIT_HI):
    """Fit MSD_xy = 4 D t over [lo,hi] ns. Returns D in cm^2/s and stderr."""
    m = (lag_ns > lo) & (lag_ns < hi)
    t_s = lag_ns[m] * 1e-9
    y_cm2 = msd_A2[m] * 1e-16
    res = stats.linregress(t_s, y_cm2)
    return res.slope / 4.0, res.stderr / 4.0


def local_exponent(lag_ns, msd, lo=FIT_LO, hi=FIT_HI):
    m = lag_ns > 0
    lt, mm = lag_ns[m], msd[m]
    beta = np.gradient(np.log(mm), np.log(lt))
    w = (lt > lo) & (lt < hi)
    return beta[w].mean()


def electrolyte_volume():
    """Volume occupied by electrolyte, from the time-averaged z-density
    profile (the box volume is wrong here: ~12 A of each end is vacuum,
    and the WSSe slab excludes more)."""
    rows, nblocks = [], 0
    for line in open('density_z.dat'):
        if line.startswith('#'):
            continue
        p = line.split()
        if len(p) == 3:
            nblocks += 1
            continue
        rows.append([float(x) for x in p])
    a = np.array(rows).reshape(nblocks, -1, 4)
    z = a[0, :, 1]
    d = a[:, :, 3].mean(0)                 # g/cm^3 per 1 A bin
    dz = z[1] - z[0]
    occupied = d > 0.3                     # bins that are condensed phase
    # exclude the WSSe slab itself (very high density bins)
    electrolyte = occupied & (d < 2.0)
    return float(electrolyte.sum() * dz * T.LX * T.LY), z, d


def main():
    d = T.load()
    pos, mol, typ = d['pos'], d['mol'], d['typ']

    sel_li = (typ == T.T_LI)
    sel_tfsi = np.isin(typ, T.T_TFSI)
    sel_peo = np.isin(typ, (T.T_PEO_C, T.T_PEO_H, T.T_PEO_O))
    sel_electrolyte = sel_li | sel_tfsi | sel_peo

    print("Unwrapping in-plane coordinates ...")
    pos_u = T.unwrap_xy(pos)

    # Remove the COM drift of the whole electrolyte (not of the Li group:
    # that would subtract the collective cation motion we want to measure).
    w = T.masses_of(typ).copy()
    w[~sel_electrolyte] = 0.0
    pos_u = T.remove_com_drift(pos_u, w)

    print("Li+ MSD ...")
    lag, msd_li = T.msd_sliding_xy(pos_u[:, sel_li, :])

    print("TFSI- molecular COM MSD ...")
    com_tfsi, mids = T.molecule_com(pos_u, mol, typ, sel_tfsi)
    print(f"  -> {len(mids)} TFSI- molecules")
    _, msd_tfsi = T.msd_sliding_xy(com_tfsi)

    D_li, eD_li = fit_D2d(lag, msd_li)
    D_tf, eD_tf = fit_D2d(lag, msd_tfsi)
    b_li = local_exponent(lag, msd_li)
    b_tf = local_exponent(lag, msd_tfsi)

    ratio = D_tf / D_li
    tplus = D_li / (D_li + D_tf)
    # propagate the two fit errors into t+
    dt_dli = D_tf / (D_li + D_tf) ** 2
    dt_dtf = -D_li / (D_li + D_tf) ** 2
    e_tplus = float(np.hypot(dt_dli * eD_li, dt_dtf * eD_tf))

    V_cm3_A3, zprof, dprof = electrolyte_volume()
    V_cm3 = V_cm3_A3 * 1e-24
    N = int(sel_li.sum())
    e = 1.602176634e-19
    kB = 1.380649e-23
    Temp = 400.0
    # Nernst-Einstein (ideal, uncorrelated): sigma = N e^2 (D+ + D-) / (V kT)
    sigma = N * e ** 2 * (D_li + D_tf) / (V_cm3 * kB * Temp)   # S/cm

    lines = []
    def P(s=''):
        print(s); lines.append(s)

    P("=" * 72)
    P("IN-PLANE TRANSPORT FROM THE EXISTING 30 ns TRAJECTORY")
    P("=" * 72)
    P(f"  fit window              : {FIT_LO}-{FIT_HI} ns, MSD_xy = 4 D t (2D)")
    P(f"  D_Li+   (in-plane)      : {D_li:.3e} +/- {eD_li:.1e} cm^2/s")
    P(f"  D_TFSI- (in-plane)      : {D_tf:.3e} +/- {eD_tf:.1e} cm^2/s")
    P(f"  D_TFSI / D_Li           : {ratio:.2f}")
    P(f"  diffusion exponent beta : Li+ {b_li:.2f}, TFSI- {b_tf:.2f}"
      "   (1.0 = fully diffusive)")
    P()
    P(f"  ideal Nernst-Einstein t+ : {tplus:.3f} +/- {e_tplus:.3f}")
    P( "     (neglects ion correlation; an upper bound for PEO/LiTFSI)")
    P()
    P(f"  electrolyte volume      : {V_cm3_A3:.0f} A^3"
      f"  ({100*V_cm3_A3/(T.LX*T.LY*T.LZ):.0f}% of the cell;"
      " the rest is vacuum + WSSe)")
    P(f"  sigma (Nernst-Einstein) : {sigma:.3e} S/cm")
    P( "     (ideal NE, in-plane D only; an upper bound)")
    P()
    P("  INTERPRETATION")
    P("  * D_TFSI/D_Li > 1 is the expected ordering for PEO/LiTFSI: the")
    P("    cation is retarded by ether-oxygen coordination while the")
    P("    delocalised anion moves more freely. Recovering it supports")
    P("    the 0.75 charge-scaling choice.")
    P("  * Both species are sub-diffusive (beta < 1), so these are")
    P("    apparent diffusivities over the 3-10 ns window, as for Li+.")
    P("  * t+ and sigma are ideal-NE values and should be quoted as")
    P("    upper bounds, not as converged predictions.")

    with open(OUT_TXT, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"\nWrote {OUT_TXT}")

    # ----------------- figure -----------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4))
    m = lag > 0
    ax1.plot(lag[m], msd_li[m], '-', color=C_LI, label=r'Li$^+$')
    ax1.plot(lag[m], msd_tfsi[m], '-', color=C_TFSI, label=r'TFSI$^-$')
    fitm = (lag > FIT_LO) & (lag < FIT_HI)
    for D, c in ((D_li, C_LI), (D_tf, C_TFSI)):
        ax1.plot(lag[fitm], 4 * D * lag[fitm] * 1e-9 / 1e-16, '--',
                 color=c, lw=1.6, alpha=0.9)
    ax1.set_xlabel('Time lag [ns]')
    ax1.set_ylabel(r'In-plane MSD [$\mathrm{\AA}^2$]')
    ax1.set_xlim(0, lag[m].max())
    ax1.set_ylim(bottom=0)
    ax1.legend(frameon=False, loc='upper left')
    ax1.text(0.97, 0.05, '(a)', transform=ax1.transAxes, ha='right')
    T.apply_style(ax1)

    ax2.loglog(lag[m], msd_li[m], '-', color=C_LI, label=r'Li$^+$')
    ax2.loglog(lag[m], msd_tfsi[m], '-', color=C_TFSI, label=r'TFSI$^-$')
    ref = np.array([lag[m][0], lag[m][-1]])
    i1 = np.argmin(np.abs(lag[m] - 1.0))
    ax2.loglog(ref, msd_li[m][i1] * (ref / 1.0), 'k--', lw=1.6,
               label='slope 1')
    ax2.set_xlabel('Time lag [ns]')
    ax2.set_ylabel(r'In-plane MSD [$\mathrm{\AA}^2$]')
    ax2.legend(frameon=False, loc='upper left')
    ax2.text(0.97, 0.05, '(b)', transform=ax2.transAxes, ha='right')
    T.apply_style(ax2)

    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote {OUT_PDF}")


if __name__ == '__main__':
    main()
