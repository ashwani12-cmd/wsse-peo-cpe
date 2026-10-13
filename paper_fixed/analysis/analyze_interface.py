#!/usr/bin/env python3
"""
==========================================================================
analyze_interface.py  --  item #2: interfacial polymer structure

Tests, rather than asserts, the "2D filler templates the polymer" claim
made in the introduction. Computed entirely from the existing
prod.lammpstrj -- no new simulation needed.

Observables, all resolved against signed z relative to the WSSe midplane
(negative = Se-face side, positive = S-face side):

  (a) component number-density profiles: PEO backbone, ether O, Li+, TFSI-
  (b) PEO chain orientational order parameter
          P2 = (3<cos^2 theta> - 1)/2
      of the vector joining successive ether oxygens along a chain, with
      theta measured from the sheet normal.
          P2 = +1  chain segments normal to the sheet
          P2 =  0  isotropic
          P2 = -0.5 chain segments lying in the sheet plane
  (c) Li+ ether-oxygen coordination number vs distance from the sheet,
      which shows whether the interface perturbs the solvation shell.

NOTE ON THE SHEET: it is corrugated (~1.9 A rms at 400 K), so a single
midplane is an approximation. Profiles are therefore smoothed over 1 A
bins and should be read for trends, not for sub-angstrom features.

USAGE
    python3 analyze_interface.py

OUTPUT
    fig6_interface.pdf
    interface_summary.txt
==========================================================================
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import traj_io as T

OUT_PDF = 'fig6_interface.pdf'
OUT_TXT = 'interface_summary.txt'

ZMAX = 34.0
DZ = 1.0
EO_CUT = 3.5          # Li-O first shell cutoff, as in the paper

plt.rcParams.update(T.PLOT_RC)
C_PEO, C_EO, C_LI, C_TFSI = '#6c7a89', '#c0392b', '#7d3c98', '#2980b9'


def zbins():
    edges = np.arange(-ZMAX, ZMAX + DZ, DZ)
    return edges, 0.5 * (edges[:-1] + edges[1:])


def density_profiles(pos, typ, zc_frames, edges):
    """Per-frame counts per bin for each component."""
    comps = {
        'PEO backbone': np.isin(typ, (T.T_PEO_C, T.T_PEO_O)),
        'ether O':      (typ == T.T_PEO_O),
        'Li$^+$':       (typ == T.T_LI),
        'TFSI$^-$':     np.isin(typ, T.T_TFSI),
    }
    out = {}
    nf = pos.shape[0]
    for name, sel in comps.items():
        h = np.zeros(len(edges) - 1)
        for fi in range(nf):
            z = pos[fi, sel, 2] - zc_frames[fi]
            h += np.histogram(z, bins=edges)[0]
        out[name] = h / nf
    return out


def chain_order_parameter(pos, mol, typ, zc_frames, edges):
    """P2 of successive ether-oxygen vectors vs z."""
    oid = np.where(typ == T.T_PEO_O)[0]
    omol = mol[oid]
    pairs = []
    for m in np.unique(omol):
        idx = oid[omol == m]
        idx = np.sort(idx)
        pairs += [(idx[k], idx[k + 1]) for k in range(len(idx) - 1)]
    pairs = np.array(pairs)
    print(f"  {len(pairs)} intra-chain O-O vectors from {len(np.unique(omol))} chains")

    nb = len(edges) - 1
    s_sum = np.zeros(nb); s_cnt = np.zeros(nb)
    box = np.array([T.LX, T.LY])
    for fi in range(pos.shape[0]):
        a = pos[fi, pairs[:, 0], :]
        b = pos[fi, pairs[:, 1], :]
        d = b - a
        d[:, :2] -= box * np.round(d[:, :2] / box)     # nearest image in-plane
        n = np.linalg.norm(d, axis=1)
        ok = n > 1e-6
        cos2 = (d[ok, 2] / n[ok]) ** 2
        zmid = 0.5 * (a[ok, 2] + b[ok, 2]) - zc_frames[fi]
        bi = np.digitize(zmid, edges) - 1
        good = (bi >= 0) & (bi < nb)
        np.add.at(s_sum, bi[good], cos2[good])
        np.add.at(s_cnt, bi[good], 1.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean_cos2 = np.where(s_cnt > 0, s_sum / np.maximum(s_cnt, 1), np.nan)
    return 0.5 * (3 * mean_cos2 - 1), s_cnt


def li_coordination_vs_z(pos, typ, zc_frames, edges):
    """Li+ ether-oxygen coordination number resolved by z."""
    li = np.where(typ == T.T_LI)[0]
    eo = np.where(typ == T.T_PEO_O)[0]
    nb = len(edges) - 1
    c_sum = np.zeros(nb); c_cnt = np.zeros(nb)
    box = np.array([T.LX, T.LY])
    for fi in range(pos.shape[0]):
        L = pos[fi, li, :]; O = pos[fi, eo, :]
        d = L[:, None, :] - O[None, :, :]
        d[:, :, 0] -= box[0] * np.round(d[:, :, 0] / box[0])
        d[:, :, 1] -= box[1] * np.round(d[:, :, 1] / box[1])
        r = np.sqrt((d ** 2).sum(-1))
        n = (r < EO_CUT).sum(1).astype(float)
        z = L[:, 2] - zc_frames[fi]
        bi = np.digitize(z, edges) - 1
        good = (bi >= 0) & (bi < nb)
        np.add.at(c_sum, bi[good], n[good])
        np.add.at(c_cnt, bi[good], 1.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(c_cnt > 0, c_sum / np.maximum(c_cnt, 1), np.nan), c_cnt


def main():
    d = T.load()
    pos, mol, typ = d['pos'], d['mol'], d['typ']
    zc = pos[:, typ == T.T_W, 2].mean(axis=1)       # sheet midplane per frame
    sface = pos[:, typ == T.T_S, 2].mean(axis=1) - zc
    seface = pos[:, typ == T.T_SE, 2].mean(axis=1) - zc
    edges, cen = zbins()

    print("Density profiles ...")
    dens = density_profiles(pos, typ, zc, edges)
    print("Chain orientation ...")
    P2, p2n = chain_order_parameter(pos, mol, typ, zc, edges)
    print("Li+ coordination vs z ...")
    cn, cnn = li_coordination_vs_z(pos, typ, zc, edges)

    lines = []
    def P(s=''):
        print(s); lines.append(s)

    P("=" * 72)
    P("INTERFACIAL STRUCTURE (signed z from WSSe midplane)")
    P(f"  S-face plane at z = {sface.mean():+.2f} A, "
      f"Se-face plane at z = {seface.mean():+.2f} A")
    P("=" * 72)

    bulk = (np.abs(cen) > 18) & (np.abs(cen) < 28)
    P()
    P("PEO chain orientation P2 (negative = chains lie parallel to sheet):")
    P(f"  {'z (A)':>8} {'P2':>8} {'n vectors':>11}")
    for i, z in enumerate(cen):
        if 2 < abs(z) <= 14 and p2n[i] > 200:
            P(f"  {z:>8.1f} {P2[i]:>8.3f} {p2n[i]:>11.0f}")
    P(f"  bulk reference (18<|z|<28): P2 = {np.nanmean(P2[bulk]):+.3f}")

    P()
    P("Li+ ether-oxygen coordination vs z:")
    P(f"  {'z (A)':>8} {'n_coord':>9} {'samples':>9}")
    for i, z in enumerate(cen):
        if abs(z) <= 14 and cnn[i] > 50:
            P(f"  {z:>8.1f} {cn[i]:>9.2f} {cnn[i]:>9.0f}")
    P(f"  bulk reference (18<|z|<28): n_coord = {np.nansum(cn[bulk]*cnn[bulk])/np.nansum(cnn[bulk]):.2f}")

    P()
    P("INTERPRETATION GUIDE")
    P("  * P2 markedly below the bulk value within ~10 A of the sheet")
    P("    means the filler templates chain alignment parallel to its")
    P("    plane -- the 'orientational template' asserted in the intro.")
    P("  * n_coord staying near its bulk value up to the interface means")
    P("    the filler perturbs the solvation shell only weakly, which is")
    P("    the structural counterpart of the Se-face preference acting on")
    P("    fully solvated cations.")

    with open(OUT_TXT, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f"\nWrote {OUT_TXT}")

    # ----------------- figure -----------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.2))
    ax = axes[0]
    for (name, h), c in zip(dens.items(), (C_PEO, C_EO, C_LI, C_TFSI)):
        scale = h / (T.LX * T.LY * DZ)
        ax.plot(cen, scale, '-', color=c, label=name)
    ax.axvspan(seface.mean(), sface.mean(), color='gold', alpha=0.35, lw=0)
    ax.set_xlabel(r'$z$ from WSSe midplane [$\mathrm{\AA}$]')
    ax.set_ylabel(r'number density [$\mathrm{\AA}^{-3}$]')
    ax.set_xlim(-ZMAX, ZMAX); ax.set_yscale('log')
    ax.legend(frameon=False, fontsize=11, loc='lower center')
    ax.text(0.97, 0.05, '(a)', transform=ax.transAxes, ha='right')
    T.apply_style(ax)

    ax = axes[1]
    good = p2n > 200
    ax.axhline(0, color='k', lw=1.2, ls=':')
    ax.axhline(np.nanmean(P2[bulk]), color='gray', lw=1.6, ls='--',
               label='bulk')
    ax.plot(cen[good], P2[good], '-', color='#117864', label=r'$P_2$')
    ax.axvspan(seface.mean(), sface.mean(), color='gold', alpha=0.35, lw=0)
    ax.set_xlabel(r'$z$ from WSSe midplane [$\mathrm{\AA}$]')
    ax.set_ylabel(r'PEO chain order parameter $P_2$')
    ax.set_xlim(-ZMAX, ZMAX)
    ax.legend(frameon=False, loc='upper right')
    ax.text(0.97, 0.05, '(b)', transform=ax.transAxes, ha='right')
    T.apply_style(ax)

    ax = axes[2]
    good = cnn > 50
    ax.plot(cen[good], cn[good], '-', color=C_LI)
    ax.axhline(np.nansum(cn[bulk]*cnn[bulk])/np.nansum(cnn[bulk]),
               color='gray', lw=1.6, ls='--', label='bulk')
    ax.axvspan(seface.mean(), sface.mean(), color='gold', alpha=0.35, lw=0)
    ax.set_xlabel(r'$z$ from WSSe midplane [$\mathrm{\AA}$]')
    ax.set_ylabel(r'Li$^+$ ether-O coordination')
    ax.set_xlim(-ZMAX, ZMAX)
    ax.legend(frameon=False, loc='lower center')
    ax.text(0.97, 0.05, '(c)', transform=ax.transAxes, ha='right')
    T.apply_style(ax)

    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close(fig)
    print(f"Wrote {OUT_PDF}")


if __name__ == '__main__':
    main()
