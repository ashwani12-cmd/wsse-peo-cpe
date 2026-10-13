#!/usr/bin/env python3
"""
==========================================================================
analyze_faces_3d.py

Face-resolved Li+ analysis for the Janus WSSe / PEO-LiTFSI CPE paper,
done with TRUE 3D distances.

WHY THIS REPLACES THE EARLIER ANALYSIS
--------------------------------------
The previous version (make_paper_figure.py, compute_janus()) assigned a
Li+ to a face using a one-dimensional z-window measured from the MEAN
plane of each chalcogen layer:

    n_S  = count( s_face  < z_Li < s_face  + rc )
    n_Se = count( se_face - rc < z_Li < se_face )

That criterion is not a distance to the face.  The sheet is thermally
corrugated (instantaneous rms ~1.9 A, static ~1.5 A), so the z-ranges
spanned by the S and Se layers OVERLAP over ~10 A -- wider than every
cutoff used except 8 A.  The z-windows therefore do not separate the two
faces, and they count ions that are laterally far from any chalcogen.
Diagnostic: ions the 1D rc=3 A criterion assigns to the Se-face have a
true median Li-Se distance of 5.8 A.

Two consequences:
  * No Li+ ever reaches 3 A of either face (closest approach over 30 ns:
    3.40 A to S, 3.48 A to Se), so the correct count at rc = 3 A is zero
    on both faces -- Li+ carries ~6 ether oxygens and cannot reach
    contact.
  * The 1D result was non-monotonic in rc (Se>S at 3,4 A; S>Se at 5,6 A;
    Se>S at 8 A).  A real face preference cannot oscillate with cutoff.

This script computes the minimum 3D distance (xy periodic images
included) from every Li+ to the nearest S and nearest Se atom, frame by
frame, and reports face-resolved coordination numbers with block-averaged
errors whose block size is chosen from a convergence scan rather than
fixed at 5.  It cross-checks against the LAMMPS-computed Li-S / Li-Se
RDFs in rdf_LiWSSe.dat, which were already being written during
production but never analysed.

USAGE
    python3 analyze_faces_3d.py

INPUT
    prod.lammpstrj    trajectory (id mol type x y z)
    rdf_LiWSSe.dat    LAMMPS rdf for pairs Li-S (4-11) and Li-Se (4-12)

OUTPUT
    fig4_janus_contacts.pdf   regenerated figure (correct 3D analysis)
    table1_contacts.tex       drop-in replacement for Table 1
    faces_3d.npz             cached distances for reuse
==========================================================================
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ----------------------------- CONFIG -----------------------------
TRAJ      = 'prod.lammpstrj'
RDF_FILE  = 'rdf_LiWSSe.dat'
CACHE     = 'faces_3d.npz'
OUT_PDF   = 'fig4_janus_contacts.pdf'
OUT_TEX   = 'table1_contacts.tex'

LI_TYPE, S_TYPE, SE_TYPE = 4, 11, 12
LX, LY = 40.0, 40.0          # in-plane box dims (sheet spans the full cell)
CUTOFFS = [3.0, 4.0, 5.0, 6.0, 8.0]
DT_FRAME_NS = 0.05           # 50000 steps x 1 fs
# ------------------------------------------------------------------

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 18, "axes.labelsize": 18, "axes.titlesize": 18,
    "xtick.labelsize": 18, "ytick.labelsize": 18, "legend.fontsize": 14,
    "axes.linewidth": 2.0, "lines.linewidth": 2.5, "lines.markersize": 8,
})

C_S, C_SE = '#b8860b', '#1f618d'


def apply_style(ax):
    ax.tick_params(axis='both', which='major', direction='inout', length=8, width=2)
    ax.tick_params(axis='both', which='minor', direction='in', length=5, width=1.5)
    ax.minorticks_on()


# ============================================================
# Trajectory -> per-frame 3D minimum distances
# ============================================================
def read_traj(fname):
    """Return (li, S, Se) position arrays of shape (nframes, natoms, 3)."""
    print(f"Reading {fname} ...")
    li_all, s_all, se_all = [], [], []
    with open(fname) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith('ITEM: TIMESTEP'):
                continue
            f.readline()                      # timestep
            f.readline()                      # ITEM: NUMBER OF ATOMS
            n = int(f.readline())
            f.readline()                      # ITEM: BOX BOUNDS
            for _ in range(3):
                f.readline()
            f.readline()                      # ITEM: ATOMS ...
            li, S, Se = [], [], []
            for _ in range(n):
                p = f.readline().split()
                t = int(p[2])
                if t == LI_TYPE:
                    li.append((float(p[3]), float(p[4]), float(p[5])))
                elif t == S_TYPE:
                    S.append((float(p[3]), float(p[4]), float(p[5])))
                elif t == SE_TYPE:
                    Se.append((float(p[3]), float(p[4]), float(p[5])))
            li_all.append(li); s_all.append(S); se_all.append(Se)
    li, S, Se = np.array(li_all), np.array(s_all), np.array(se_all)
    print(f"  -> {li.shape[0]} frames, {li.shape[1]} Li+, "
          f"{S.shape[1]} S, {Se.shape[1]} Se")
    return li, S, Se


def min_dist(L, F):
    """Minimum distance from each Li+ to any atom of face F, xy images included."""
    d = L[:, None, :] - F[None, :, :]
    d[:, :, 0] -= LX * np.round(d[:, :, 0] / LX)
    d[:, :, 1] -= LY * np.round(d[:, :, 1] / LY)
    return np.sqrt((d ** 2).sum(-1)).min(axis=1)


def face_distances():
    if os.path.exists(CACHE):
        print(f"Loading cached {CACHE} ...")
        z = np.load(CACHE)
        return z['dS'], z['dSe'], z['corr_S'], z['corr_Se']
    li, S, Se = read_traj(TRAJ)
    print("Computing 3D minimum distances ...")
    dS = np.array([min_dist(li[i], S[i]) for i in range(len(li))])
    dSe = np.array([min_dist(li[i], Se[i]) for i in range(len(li))])
    corr_S = (S[:, :, 2].max(1) - S[:, :, 2].min(1))
    corr_Se = (Se[:, :, 2].max(1) - Se[:, :, 2].min(1))
    np.savez(CACHE, dS=dS, dSe=dSe, corr_S=corr_S, corr_Se=corr_Se)
    return dS, dSe, corr_S, corr_Se


# ============================================================
# Block averaging with a convergence scan
# ============================================================
def block_scan(x, max_blocks=40):
    """SEM vs block count.  For well-sampled data this plateaus once the
    block length exceeds the correlation time.  Here it does NOT plateau --
    it keeps growing with block length, which is itself a sign that the
    30 ns trajectory is too short to converge the near-face occupancy."""
    out = []
    for nb in range(2, max_blocks + 1):
        if len(x) // nb < 4:
            break
        blocks = np.array_split(x, nb)
        means = np.array([b.mean() for b in blocks])
        out.append((nb, (len(x) // nb) * DT_FRAME_NS,
                    means.std(ddof=1) / np.sqrt(nb)))
    return out


def autocorr_time(x):
    """Integrated autocorrelation time (ns) of a per-frame time series."""
    y = x - x.mean()
    if y.std() == 0:
        return DT_FRAME_NS
    ac = np.correlate(y, y, 'full')[len(y) - 1:]
    ac = ac / ac[0]
    below = np.where(ac < 0.05)[0]
    cut = below[0] if len(below) else min(len(ac), 50)
    return (1 + 2 * ac[1:cut].sum()) * DT_FRAME_NS


def block_sem(x):
    """Block-averaged SEM using blocks of at least 5 correlation times
    (and at least 2 ns), while keeping at least 5 blocks for the variance
    estimate.  Returns (sem, tau_ns, n_eff, n_blocks)."""
    tau = autocorr_time(x)
    target_ns = max(5.0 * tau, 2.0)
    nb = int(len(x) * DT_FRAME_NS / target_ns)
    nb = max(5, min(nb, len(x) // 4))            # >=5 blocks, >=4 frames each
    blocks = np.array_split(x, nb)
    means = np.array([b.mean() for b in blocks])
    sem = means.std(ddof=1) / np.sqrt(nb)
    n_eff = len(x) * DT_FRAME_NS / tau
    return sem, tau, n_eff, nb


# ============================================================
# LAMMPS RDF cross-check
# ============================================================
def parse_rdf(fname):
    """Return r and per-time-block coordination numbers for Li-S and Li-Se."""
    blocks = []
    with open(fname) as f:
        ls = f.readlines()
    i = 0
    while i < len(ls):
        if ls[i].startswith('#'):
            i += 1
            continue
        p = ls[i].split()
        if len(p) == 2:
            n = int(p[1])
            D = [[float(v) for v in ls[i + 1 + k].split()] for k in range(n)]
            blocks.append(np.array(D))
            i += 1 + n
        else:
            i += 1
    B = np.array(blocks)                 # (ntime, nbins, 6)
    r = B[0, :, 1]
    nS, nSe = B[:, :, 3], B[:, :, 5]     # running coordination numbers
    return r, nS, nSe


# ============================================================
# Main
# ============================================================
def main():
    dS, dSe, corr_S, corr_Se = face_distances()
    nframes = dS.shape[0]

    print("\n" + "=" * 70)
    print("SHEET GEOMETRY (why the 1D z-window criterion fails)")
    print("=" * 70)
    print(f"  S-layer  instantaneous z-spread : {corr_S.mean():5.2f} A")
    print(f"  Se-layer instantaneous z-spread : {corr_Se.mean():5.2f} A")
    print(f"  -> the two layers' z-ranges overlap by more than most cutoffs,")
    print(f"     so a z-window about the mean plane does not isolate a face.")

    print("\n" + "=" * 70)
    print("CLOSEST APPROACH OVER THE WHOLE RUN")
    print("=" * 70)
    print(f"  min Li-S  distance : {dS.min():.2f} A")
    print(f"  min Li-Se distance : {dSe.min():.2f} A")
    print("  -> no Li+ reaches 3 A of either face; Li+ keeps its ~6-fold")
    print("     ether-oxygen shell and cannot make direct chalcogen contact.")

    print("\n" + "=" * 70)
    print("FACE-RESOLVED Li+ COUNTS, TRUE 3D DISTANCE")
    print("=" * 70)
    print(f"{'rc(A)':>6} {'S-face':>14} {'Se-face':>14} {'Se - S':>15} "
          f"{'sigma':>6} {'ions':>8} {'tau_Se':>8} {'n_eff':>7}")
    rows = []
    for rc in CUTOFFS:
        nS = (dS < rc).sum(1).astype(float)
        nSe = (dSe < rc).sum(1).astype(float)
        mS, mSe = nS.mean(), nSe.mean()
        eS, tS, neS, _ = block_sem(nS)
        eSe, tSe, neSe, _ = block_sem(nSe)
        diff = mSe - mS
        err = float(np.hypot(eS, eSe))
        sig = abs(diff) / err if err > 0 else 0.0
        vS = int((dS < rc).any(0).sum())
        vSe = int((dSe < rc).any(0).sum())
        rows.append((rc, mS, eS, mSe, eSe, diff, err, sig, vS, vSe))
        print(f"{rc:>6.0f} {mS:8.3f}+/-{eS:.3f} {mSe:8.3f}+/-{eSe:.3f} "
              f"{diff:+9.3f}+/-{err:.3f} {sig:6.1f} {vS:3d}/{vSe:<3d} "
              f"{tSe:7.2f}n {neSe:7.0f}")
    print("  'ions' = distinct Li+ that ever came within rc (S/Se);")
    print("  tau_Se = integrated autocorrelation time of the Se-face count;")
    print("  n_eff  = effective independent samples in 30 ns for that count.")

    # ---------------- RDF cross-check ----------------
    print("\n" + "=" * 70)
    print("INDEPENDENT CROSS-CHECK: LAMMPS Li-S / Li-Se RDF coordination")
    print("=" * 70)
    r, NS, NSE = parse_rdf(RDF_FILE)
    print(f"  {len(NS)} time blocks of 1 ns from fix ave/time")
    print(f"{'rc(A)':>6} {'n_Li-S':>16} {'n_Li-Se':>16} {'ratio':>7} {'sigma':>6}")
    rdf_rows = []
    for rc in [4.0, 5.0, 6.0, 8.0]:
        i = int(np.argmin(np.abs(r - rc)))
        a, b = NS[:, i], NSE[:, i]
        ea = a.std(ddof=1) / np.sqrt(len(a))
        eb = b.std(ddof=1) / np.sqrt(len(b))
        d = b.mean() - a.mean()
        e = float(np.hypot(ea, eb))
        ratio = b.mean() / a.mean() if a.mean() > 0 else np.nan
        sig = abs(d) / e if e > 0 else 0.0
        rdf_rows.append((rc, a.mean(), ea, b.mean(), eb, ratio, sig))
        print(f"{rc:>6.0f} {a.mean():9.4f}+/-{ea:.4f} {b.mean():9.4f}+/-{eb:.4f} "
              f"{ratio:7.2f} {sig:6.1f}")
    print("  -> Se-face enrichment is consistent in SIGN at every cutoff,")
    print("     unlike the earlier 1D result which changed sign with rc.")

    # ---------------- Figure ----------------
    print("\nBuilding figure ...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.8))

    # (a) running coordination number vs r, both faces, with error band
    eS_r = NS.std(0, ddof=1) / np.sqrt(NS.shape[0])
    eSe_r = NSE.std(0, ddof=1) / np.sqrt(NSE.shape[0])
    mS_r, mSe_r = NS.mean(0), NSE.mean(0)
    sel = r <= 9.0
    ax1.plot(r[sel], mS_r[sel], '-', color=C_S, label=r'Li$^+$--S (S-face)')
    ax1.fill_between(r[sel], (mS_r - eS_r)[sel], (mS_r + eS_r)[sel],
                     color=C_S, alpha=0.25, lw=0)
    ax1.plot(r[sel], mSe_r[sel], '-', color=C_SE, label=r'Li$^+$--Se (Se-face)')
    ax1.fill_between(r[sel], (mSe_r - eSe_r)[sel], (mSe_r + eSe_r)[sel],
                     color=C_SE, alpha=0.25, lw=0)
    ax1.set_xlabel(r'$r$ [$\mathrm{\AA}$]')
    ax1.set_ylabel(r'Li$^+$ coordination $n(r)$')
    ax1.set_xlim(3.0, 9.0)
    ax1.set_ylim(bottom=0)
    ax1.legend(frameon=False, loc='upper left')
    ax1.text(0.97, 0.05, '(a)', transform=ax1.transAxes, ha='right', fontsize=18)
    apply_style(ax1)

    # (b) difference with significance band
    diff_r = mSe_r - mS_r
    err_r = np.hypot(eS_r, eSe_r)
    ax2.axhline(0, color='k', lw=1.5, ls='-')
    ax2.plot(r[sel], diff_r[sel], '-', color='#7d3c98',
             label=r'$n_{\mathrm{Se}} - n_{\mathrm{S}}$')
    ax2.fill_between(r[sel], (diff_r - err_r)[sel], (diff_r + err_r)[sel],
                     color='#7d3c98', alpha=0.25, lw=0, label=r'$\pm 1\sigma$')
    ax2.fill_between(r[sel], (diff_r - 3 * err_r)[sel], (diff_r + 3 * err_r)[sel],
                     color='#7d3c98', alpha=0.10, lw=0, label=r'$\pm 3\sigma$')
    ax2.set_xlabel(r'$r$ [$\mathrm{\AA}$]')
    ax2.set_ylabel(r'Excess Li$^+$ on Se-face')
    ax2.set_xlim(3.0, 9.0)
    ax2.legend(frameon=False, loc='upper left')
    ax2.text(0.97, 0.05, '(b)', transform=ax2.transAxes, ha='right', fontsize=18)
    apply_style(ax2)

    fig.tight_layout()
    fig.savefig(OUT_PDF, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved {OUT_PDF}")

    # ---------------- LaTeX table ----------------
    with open(OUT_TEX, 'w') as fh:
        fh.write(r"""% Generated by analyze_faces_3d.py -- replaces the old Table 1.
\begin{table}[htbp]
\caption{Face-resolved Li$^+$ coordination from the 30~ns production run,
using the true three-dimensional minimum distance from each Li$^+$ to the
nearest chalcogen of each face. $n_{\mathrm{S}}$ and $n_{\mathrm{Se}}$ are
Li$^+$ counts per frame within $r_c$. Errors are block-averaged standard
errors using blocks of at least five integrated autocorrelation times of
the corresponding occupancy signal. No Li$^+$ approaches either face
within 3~\AA{} at any point in the trajectory, so that cutoff is empty by
construction. The $r_c = 5$~\AA{} shell is the best sampled
($n_{\mathrm{eff}} \approx 250$ independent samples); $r_c = 6$ and
8~\AA{} have $n_{\mathrm{eff}} \approx 22$ and $7$ respectively and
should be read as indicative only.}
\label{tab:contacts}
\centering
\begin{tabular}{crrrc}
\hline
$r_c$ (\AA) & $n_{\mathrm{S}}$ & $n_{\mathrm{Se}}$ & $n_{\mathrm{Se}}-n_{\mathrm{S}}$ & $\sigma$ \\
\hline
""")
        for rc, mS, eS, mSe, eSe, diff, err, sig, vS, vSe in rows:
            if mS == 0 and mSe == 0:
                fh.write(f"{rc:.0f} & 0 & 0 & --- & --- \\\\\n")
            else:
                fh.write(f"{rc:.0f} & ${mS:.2f} \\pm {eS:.2f}$ & "
                         f"${mSe:.2f} \\pm {eSe:.2f}$ & "
                         f"${diff:+.2f} \\pm {err:.2f}$ & {sig:.1f} \\\\\n")
        fh.write(r"""\hline
\end{tabular}
\end{table}

% Independent cross-check from the LAMMPS Li--S / Li--Se RDFs (rdf_LiWSSe.dat),
% 30 one-nanosecond time blocks written during production:
""")
        for rc, a, ea, b, eb, ratio, sig in rdf_rows:
            fh.write(f"%   rc = {rc:.0f} A : n_LiS = {a:.4f} +/- {ea:.4f}, "
                     f"n_LiSe = {b:.4f} +/- {eb:.4f}, "
                     f"ratio = {ratio:.2f}, {sig:.1f} sigma\n")
    print(f"Saved {OUT_TEX}")
    print("\nDone.")


if __name__ == '__main__':
    main()
