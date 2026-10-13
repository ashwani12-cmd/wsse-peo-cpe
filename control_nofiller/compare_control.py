#!/usr/bin/env python3
"""
==========================================================================
compare_control.py -- composite vs filler-free reference

Run this AFTER the control job finishes. It computes the same observables
on both systems and emits the side-by-side comparison plus a LaTeX table,
which is what converts statements of the form "the composite has X" into
"the filler causes X".

USAGE
    python3 compare_control.py                      # uses default paths
    python3 compare_control.py <composite> <control>

Both directories must contain prod.lammpstrj and rdf_LiO.dat.

OUTPUT
    control_comparison.txt
    table2_control.tex
==========================================================================
"""

import os
import sys
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
COMPOSITE = sys.argv[1] if len(sys.argv) > 2 else os.path.join(HERE, '..', 'npt_shrunk_restart')
CONTROL = sys.argv[2] if len(sys.argv) > 2 else HERE
ANALYSIS = os.path.join(HERE, '..', 'npt_shrunk_restart')
sys.path.insert(0, ANALYSIS)

FIT_LO, FIT_HI = 3.0, 10.0


def fit_D2d(lag, msd):
    m = (lag > FIT_LO) & (lag < FIT_HI)
    r = stats.linregress(lag[m] * 1e-9, msd[m] * 1e-16)
    return r.slope / 4.0


def parse_rdf_coord(path, cut=3.5):
    blocks, ls, i = [], open(path).readlines(), 0
    while i < len(ls):
        if ls[i].startswith('#'):
            i += 1; continue
        p = ls[i].split()
        if len(p) == 2:
            n = int(p[1])
            blocks.append(np.array([[float(x) for x in ls[i + 1 + k].split()]
                                    for k in range(n)]))
            i += 1 + n
        else:
            i += 1
    B = np.array(blocks)
    a = B.mean(0); r = a[:, 1]
    j = int(np.argmin(np.abs(r - cut)))
    pk = int(np.argmax(a[:, 2]))
    per_block = B[:, j, 3]
    return r[pk], a[j, 3], per_block.std(ddof=1) / np.sqrt(len(per_block))


def analyse(d, label):
    """Return dict of observables for one system directory."""
    print(f"\n### {label}: {d}")
    cwd = os.getcwd()
    os.chdir(d)
    for m in ('traj_io',):
        if m in sys.modules:
            del sys.modules[m]
    import traj_io as T

    data = T.load()
    pos, mol, typ = data['pos'], data['mol'], data['typ']
    sel_li = (typ == T.T_LI)
    sel_tf = np.isin(typ, T.T_TFSI)
    sel_el = sel_li | sel_tf | np.isin(typ, (T.T_PEO_C, T.T_PEO_H, T.T_PEO_O))

    pos_u = T.unwrap_xy(pos)
    w = T.masses_of(typ).copy(); w[~sel_el] = 0.0
    pos_u = T.remove_com_drift(pos_u, w)

    lag, msd_li = T.msd_sliding_xy(pos_u[:, sel_li, :])
    com, _ = T.molecule_com(pos_u, mol, typ, sel_tf)
    _, msd_tf = T.msd_sliding_xy(com)
    D_li, D_tf = fit_D2d(lag, msd_li), fit_D2d(lag, msd_tf)

    pk, ncoord, encoord = parse_rdf_coord('rdf_LiO.dat')

    # bulk PEO chain order parameter (far from any interface)
    import analyze_interface as AI
    has_sheet = (typ == T.T_W).sum() > 0
    zc = (pos[:, typ == T.T_W, 2].mean(axis=1) if has_sheet
          else np.full(pos.shape[0], 75.9))
    edges, cen = AI.zbins()
    P2, p2n = AI.chain_order_parameter(pos, mol, typ, zc, edges)
    bulk = (np.abs(cen) > 18) & (np.abs(cen) < 28) & (p2n > 200)
    near = (np.abs(cen) > 3) & (np.abs(cen) < 8) & (p2n > 200)

    os.chdir(cwd)
    return dict(D_li=D_li, D_tf=D_tf, ratio=D_tf / D_li,
                tplus=D_li / (D_li + D_tf), peak=pk,
                ncoord=ncoord, encoord=encoord,
                P2_bulk=np.nanmean(P2[bulk]), P2_near=np.nanmean(P2[near]),
                has_sheet=has_sheet)


def main():
    for d in (COMPOSITE, CONTROL):
        for f in ('prod.lammpstrj', 'rdf_LiO.dat'):
            if not os.path.exists(os.path.join(d, f)):
                sys.exit(f"missing {f} in {d} -- has the control run finished?")

    C = analyse(COMPOSITE, 'composite (with WSSe)')
    R = analyse(CONTROL, 'control (filler-free)')

    lines = []
    def P(s=''):
        print(s); lines.append(s)

    P()
    P("=" * 76)
    P("COMPOSITE vs FILLER-FREE CONTROL")
    P("=" * 76)
    P(f"{'observable':<34}{'composite':>14}{'control':>14}{'change':>13}")
    P("-" * 76)

    def row(name, a, b, fmt='{:.3f}', pct=True):
        ch = (f"{100*(a-b)/b:+.1f}%" if (pct and b) else f"{a-b:+.3f}")
        P(f"{name:<34}{fmt.format(a):>14}{fmt.format(b):>14}{ch:>13}")

    row('D_xy(Li+)  [1e-7 cm2/s]', C['D_li']*1e7, R['D_li']*1e7, '{:.2f}')
    row('D_xy(TFSI-) [1e-7 cm2/s]', C['D_tf']*1e7, R['D_tf']*1e7, '{:.2f}')
    row('D_TFSI / D_Li', C['ratio'], R['ratio'], '{:.2f}')
    row('ideal NE t+', C['tplus'], R['tplus'], '{:.3f}')
    row('Li-O peak position [A]', C['peak'], R['peak'], '{:.3f}')
    row('n_coord (r<3.5 A)', C['ncoord'], R['ncoord'], '{:.2f}')
    row('PEO P2, bulk', C['P2_bulk'], R['P2_bulk'], '{:+.3f}', pct=False)
    row('PEO P2, 3-8 A from sheet', C['P2_near'], R['P2_near'], '{:+.3f}', pct=False)
    P("-" * 76)
    P()
    P("HOW TO READ THIS")
    P("  * D_xy(Li+): a difference here is the filler's effect on in-plane")
    P("    mobility. The composite is expected to be SLOWER if the filler")
    P("    obstructs, or similar if it only blocks the cross-plane channel.")
    P("  * n_coord and the Li-O peak should be UNCHANGED. If they are, the")
    P("    claim that the filler is a structural rather than chemical")
    P("    perturbation is established rather than asserted.")
    P("  * P2 near the sheet vs the control's bulk value is the clean test")
    P("    of templating: in the control there is no sheet, so its 'near'")
    P("    and 'bulk' values should agree with each other.")
    P("  * D_TFSI/D_Li in the control is the force-field validation against")
    P("    neat PEO/LiTFSI literature, free of any filler contribution.")

    with open(os.path.join(HERE, 'control_comparison.txt'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')

    with open(os.path.join(HERE, 'table2_control.tex'), 'w') as fh:
        fh.write(r"""% Generated by compare_control.py
\begin{table}[htbp]
\caption{Composite versus filler-free reference under an identical force
field, thermostat and sampling protocol. Diffusivities are in-plane
apparent values from the $3$--$10$~ns window. $P_2$ is the PEO chain
orientational order parameter; in the filler-free system the two spatial
ranges are equivalent by construction and their agreement is a
consistency check.}
\label{tab:control}
\centering
\begin{tabular}{lrr}
\hline
Observable & Composite & Filler-free \\
\hline
""")
        fh.write(f"$\\Dxy(\\mathrm{{Li}}^+)$ ($10^{{-7}}$~cm$^2$/s) & {C['D_li']*1e7:.2f} & {R['D_li']*1e7:.2f} \\\\\n")
        fh.write(f"$D_{{xy}}(\\mathrm{{TFSI}}^-)$ ($10^{{-7}}$~cm$^2$/s) & {C['D_tf']*1e7:.2f} & {R['D_tf']*1e7:.2f} \\\\\n")
        fh.write(f"$D_{{\\mathrm{{TFSI}}}}/\\DLi$ & {C['ratio']:.2f} & {R['ratio']:.2f} \\\\\n")
        fh.write(f"Ideal Nernst--Einstein $t_+$ & {C['tplus']:.2f} & {R['tplus']:.2f} \\\\\n")
        fh.write(f"Li--O first-shell peak (\\AA) & {C['peak']:.2f} & {R['peak']:.2f} \\\\\n")
        fh.write(f"$n_{{\\mathrm{{coord}}}}$ ($r<3.5$~\\AA) & {C['ncoord']:.2f} & {R['ncoord']:.2f} \\\\\n")
        fh.write(f"PEO $P_2$, bulk & {C['P2_bulk']:+.3f} & {R['P2_bulk']:+.3f} \\\\\n")
        fh.write(f"PEO $P_2$, $3$--$8$~\\AA{{}} from sheet & {C['P2_near']:+.3f} & {R['P2_near']:+.3f} \\\\\n")
        fh.write("\\hline\n\\end{tabular}\n\\end{table}\n")

    print("\nWrote control_comparison.txt and table2_control.tex")


if __name__ == '__main__':
    main()
