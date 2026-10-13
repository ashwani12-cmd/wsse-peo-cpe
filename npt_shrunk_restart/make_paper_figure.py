#!/usr/bin/env python3
"""
==========================================================================
make_paper_figures_pdf.py

Produces FOUR SEPARATE journal-quality PDF figures for the
WSSe/PEO/LiTFSI CPE paper:

    fig1_msd_anisotropy.pdf   - Directional Li+ MSD + diffusivities
    fig2_LiO_coordination.pdf - Li-O RDF + coordination number
    fig3_Li_zdensity.pdf      - Li+ z-density profile around Janus sheet
    fig4_janus_contacts.pdf   - Face-resolved Li+ contact asymmetry

USAGE:
    python3 make_paper_figures_pdf.py

REQUIRED INPUT FILES (same folder):
    prod.lammpstrj   - trajectory dump (id mol type x y z)
    rdf_LiO.dat      - Li-O RDF (LAMMPS ave/time vector format)
==========================================================================
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

# ----------------------------- CONFIG -----------------------------
TRAJ_FILE = 'prod.lammpstrj'
RDF_LIO   = 'rdf_LiO.dat'

ATOM_TYPES = {'Li': 4, 'S': 11, 'Se': 12, 'O': 3}
LX, LY, LZ = 40.0, 40.0, 82.5976
FIT_LO, FIT_HI = 3.0, 10.0
CUTOFFS = [3.0, 4.0, 5.0, 6.0, 8.0]
NBLOCKS = 5
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
    ax.xaxis.set_major_locator(MaxNLocator(6))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    ax.set_aspect("auto")


# ============================================================
# DATA READERS / ANALYSIS
# ============================================================
def read_trajectory(fname):
    print(f"Reading {fname} ...")
    with open(fname) as f:
        lines = f.readlines()
    frames = []
    i = 0
    while i < len(lines):
        if 'ITEM: TIMESTEP' in lines[i]:
            ts = int(lines[i + 1])
            natoms = int(lines[i + 3])
            li_xyz, s_z, se_z = [], [], []
            for k in range(natoms):
                p = lines[i + 9 + k].split()
                t = int(p[2])
                x, y, z = float(p[3]), float(p[4]), float(p[5])
                if t == ATOM_TYPES['Li']:
                    li_xyz.append([x, y, z])
                elif t == ATOM_TYPES['S']:
                    s_z.append(z)
                elif t == ATOM_TYPES['Se']:
                    se_z.append(z)
            frames.append({'ts': ts, 'li': np.array(li_xyz),
                           's_face': np.mean(s_z), 'se_face': np.mean(se_z)})
            i += 9 + natoms
        else:
            i += 1
    print(f"  -> {len(frames)} frames, steps {frames[0]['ts']} to {frames[-1]['ts']}")
    return frames


def compute_anisotropy(frames):
    print("Computing directional MSD ...")
    box = np.array([LX, LY, LZ])
    nf, ni = len(frames), len(frames[0]['li'])
    pos = np.zeros((nf, ni, 3))
    for fi, fr in enumerate(frames):
        pos[fi] = fr['li']
    for j in range(ni):
        for fi in range(1, nf):
            d = pos[fi, j] - pos[fi - 1, j]
            d -= box * np.round(d / box)
            pos[fi, j] = pos[fi - 1, j] + d
    pos -= pos.mean(axis=1, keepdims=True)

    max_lag = nf // 2
    msd_xy = np.zeros(max_lag); msd_z = np.zeros(max_lag)
    exy = np.zeros(max_lag); ez = np.zeros(max_lag)
    for lag in range(1, max_lag):
        dxy, dz = [], []
        for t0 in range(nf - lag):
            d = pos[t0 + lag] - pos[t0]
            dxy.append(np.mean(d[:, 0]**2 + d[:, 1]**2))
            dz.append(np.mean(d[:, 2]**2))
        msd_xy[lag], msd_z[lag] = np.mean(dxy), np.mean(dz)
        exy[lag] = np.std(dxy) / np.sqrt(len(dxy))
        ez[lag] = np.std(dz) / np.sqrt(len(dz))

    lag_ns = (frames[1]['ts'] - frames[0]['ts']) / 1e6
    lag_t = np.arange(max_lag) * lag_ns
    lag_s = lag_t * 1e-9
    mask = (lag_t > FIT_LO) & (lag_t < FIT_HI)
    sxy, _, rxy, _, _ = stats.linregress(lag_s[mask], msd_xy[mask] * 1e-20)
    sz, _, rz, _, _ = stats.linregress(lag_s[mask], msd_z[mask] * 1e-20)
    Dxy, Dz = sxy / 4, sz / 2
    print(f"  D_xy={Dxy*1e4:.3e}  D_z={Dz*1e4:.3e}  ratio={Dxy/Dz:.2f}")
    return dict(lag_t=lag_t, msd_xy=msd_xy, msd_z=msd_z, exy=exy, ez=ez, Dxy=Dxy, Dz=Dz)


def compute_janus(frames):
    print("Computing Janus contacts ...")
    nf = len(frames)

    def block_err(x, nb=NBLOCKS):
        bs = len(x) // nb
        m = [x[j*bs:(j+1)*bs].mean() for j in range(nb)]
        return np.mean(m), np.std(m) / np.sqrt(nb)

    results = {}
    for rc in CUTOFFS:
        ns = np.array([np.sum((fr['li'][:, 2] > fr['s_face']) &
                              (fr['li'][:, 2] < fr['s_face'] + rc)) for fr in frames], float)
        nse = np.array([np.sum((fr['li'][:, 2] < fr['se_face']) &
                               (fr['li'][:, 2] > fr['se_face'] - rc)) for fr in frames], float)
        ms, es = block_err(ns); mse, ese = block_err(nse)
        diff = ms - mse; err = np.hypot(es, ese)
        sigma = abs(diff) / err if err > 0 else 0
        results[rc] = (ms, es, mse, ese, sigma)
        print(f"  {rc:>4}A  S={ms:.2f}+/-{es:.2f}  Se={mse:.2f}+/-{ese:.2f}  {sigma:.1f}sigma")

    rel = np.concatenate([fr['li'][:, 2] - (fr['s_face'] + fr['se_face'])/2 for fr in frames])
    hist, edges = np.histogram(rel, bins=100, range=(-40, 40))
    zc = (edges[:-1] + edges[1:]) / 2
    zh = hist / nf
    s_off = np.mean([fr['s_face'] - (fr['s_face'] + fr['se_face'])/2 for fr in frames])
    se_off = np.mean([fr['se_face'] - (fr['s_face'] + fr['se_face'])/2 for fr in frames])
    return dict(results=results, zc=zc, zh=zh, s_off=s_off, se_off=se_off)


def parse_rdf(fname):
    print(f"Reading {fname} ...")
    blocks = []
    with open(fname) as f:
        ls = f.readlines()
    i = 0
    while i < len(ls):
        if ls[i].startswith('#'):
            i += 1; continue
        p = ls[i].split()
        if len(p) == 2 and p[1].isdigit() and int(p[1]) > 10:
            n = int(p[1])
            D = [[float(x) for x in ls[i + 1 + k].split()] for k in range(n)]
            blocks.append(np.array(D)); i += 1 + n
        else:
            i += 1
    avg = np.array(blocks).mean(axis=0)
    r, g, coord = avg[:, 1], avg[:, 2], avg[:, 3]
    idx = np.argmin(np.abs(r - 3.5)); peak = np.argmax(g)
    print(f"  peak {r[peak]:.2f} A, coord {coord[idx]:.2f}")
    return dict(r=r, g=g, coord=coord, peak_pos=r[peak], cn=coord[idx])


# ============================================================
# FIGURE 1 - Directional MSD
# ============================================================
def figure1_msd(a):
    fig, ax = plt.subplots(figsize=(7, 6))
    lag, mxy, mz = a['lag_t'], a['msd_xy'], a['msd_z']
    m = lag > 0
    ax.plot(lag[m], mxy[m], '-', color='#c0392b', label='In-plane ($xy$)')
    ax.plot(lag[m], mz[m], '-', color='#2980b9', label='Out-of-plane ($z$)')
    ax.fill_between(lag[m], mxy[m]-a['exy'][m], mxy[m]+a['exy'][m], color='#c0392b', alpha=0.2)
    ax.fill_between(lag[m], mz[m]-a['ez'][m], mz[m]+a['ez'][m], color='#2980b9', alpha=0.2)
    ax.set_xlabel('Time lag [ns]')
    ax.set_ylabel(r'MSD [$\mathrm{\AA}^2$]')
    ax.set_xlim(0, lag[m].max())
    ax.set_ylim(0, mxy[m].max() * 1.12)
    Dxy_u = a['Dxy'] * 1e11; Dz_u = a['Dz'] * 1e11
    ax.text(0.05, 0.72,
            f'$D_{{xy}}$ = {Dxy_u:.2f}$\\times10^{{-7}}$ cm²/s\n'
            f'$D_{{z}}$ = {Dz_u:.2f}$\\times10^{{-7}}$ cm²/s\n'
            f'$D_{{xy}}/D_{{z}}$ = {Dxy_u/Dz_u:.2f}',
            transform=ax.transAxes, fontsize=14, va='top',
            bbox=dict(boxstyle='round', fc='#f9e79f', ec='k', lw=1.5))
    ax.legend(frameon=False, loc='upper left')
    apply_style(ax)
    fig.tight_layout()
    fig.savefig('fig1_msd_anisotropy.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Saved fig1_msd_anisotropy.pdf")


# ============================================================
# FIGURE 2 - Li-O RDF + coordination
# ============================================================
def figure2_lio(o):
    fig, ax = plt.subplots(figsize=(7, 6))
    axb = ax.twinx()
    r, g, cn = o['r'], o['g'], o['coord']
    l1, = ax.plot(r, g, color='#16a085', label='$g(r)$')
    l2, = axb.plot(r, cn, color='#d35400', ls='--', label='$n(r)$')
    ax.set_xlim(1.5, 6)
    ax.set_ylim(0, max(g) * 1.15)
    axb.set_ylim(0, 8)
    ax.set_xlabel(r'$r$ [$\mathrm{\AA}$]',fontsize=24)
    ax.set_ylabel(r'$g_{\mathrm{Li-O}}(r)$', color='#16a085',fontsize=24)
    axb.set_ylabel('Coordination number $n(r)$', color='#d35400',fontsize=24)
    ax.axvline(o['peak_pos'], color='gray', ls=':', lw=1.5)
    ax.text(o['peak_pos'] + 0.55, max(g) * 0.52,
            f"peak: {o['peak_pos']:.2f} Å\n$n_{{coord}}$ = {o['cn']:.2f}\n($r < 3.5$ Å)",
            fontsize=14, bbox=dict(boxstyle='round', fc='white', ec='k', lw=1.5))
    ax.legend([l1, l2], ['$g(r)$', '$n(r)$'], frameon=False, loc='center right')
    apply_style(ax)
    axb.tick_params(axis='both', which='major', direction='in', length=8, width=2)
    fig.tight_layout()
    fig.savefig('fig2_LiO_coordination.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Saved fig2_LiO_coordination.pdf")


# ============================================================
# FIGURE 3 - Li z-density profile
# ============================================================
def figure3_zdensity(j):
    fig, ax = plt.subplots(figsize=(7, 6))
    zc, zh = j['zc'], j['zh']
    s_off, se_off = j['s_off'], j['se_off']
    ax.fill_between(zc, zh, color='#8e44ad', alpha=0.35)
    ax.plot(zc, zh, color='#8e44ad')
    ax.axvspan(se_off, s_off, color='#f1c40f', alpha=0.5, label='WSSe sheet')
    ax.axvline(s_off, color='#b7950b', ls='--', lw=1.5)
    ax.axvline(se_off, color='#1f618d', ls='--', lw=1.5)
    ax.text(s_off + 1.5, max(zh) * 0.92, 'S-face', fontsize=15, color='#7d6608')
    ax.text(se_off - 1.8, max(zh) * 1, 'Se-face', fontsize=15, color='#1f618d', ha='right')
    ax.set_xlabel(r'Li$^+$ $z$ rel. to sheet midplane [$\mathrm{\AA}$]')
    ax.set_ylabel('Li$^+$ count / frame / bin')
    ax.set_xlim(-35, 35)
    ax.set_ylim(0, max(zh) * 1.1)
    ax.legend(frameon=False, loc='upper right')
    apply_style(ax)
    fig.tight_layout()
    fig.savefig('fig3_Li_zdensity.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Saved fig3_Li_zdensity.pdf")


# ============================================================
# FIGURE 4 - Janus contacts
# ============================================================
def figure4_janus(j):
    fig, ax = plt.subplots(figsize=(7, 6))
    res = j['results']
    S = [res[c][0] for c in CUTOFFS]
    Se_ = [res[c][1] for c in CUTOFFS]
    SE = [res[c][2] for c in CUTOFFS]
    SEe = [res[c][3] for c in CUTOFFS]
    sig = [res[c][4] for c in CUTOFFS]
    x = np.arange(len(CUTOFFS)); w = 0.38
    ax.bar(x - w/2, S, w, yerr=Se_, capsize=5, color='#b7950b', alpha=0.85,
           edgecolor='k', linewidth=1.5, label='S-face')
    ax.bar(x + w/2, SE, w, yerr=SEe, capsize=5, color='#1f618d', alpha=0.85,
           edgecolor='k', linewidth=1.5, label='Se-face')
    ax.set_xticks(x)
    ax.set_xticklabels([f'< {int(c)} Å' for c in CUTOFFS])
    ax.set_xlabel('Distance from face')
    ax.set_ylabel(r'Li$^+$ contacts / frame')
    ax.legend(frameon=False, loc='upper left')
    ax.set_ylim(0, max(max(S) + max(Se_), max(SE) + max(SEe)) * 1.35)
    for i, (s, c) in enumerate(zip(sig, CUTOFFS)):
        mk = "***" if s > 3 else "**" if s > 2 else "*" if s > 1 else "ns"
        ax.text(i, max(S[i], SE[i]) + max(Se_[i], SEe[i]) + 0.2,
                f'{s:.1f}$\\sigma$\n{mk}', ha='center', fontsize=12, weight='bold')
    apply_style(ax)
    ax.xaxis.set_major_locator(plt.FixedLocator(x))  # keep categorical ticks
    fig.tight_layout()
    fig.savefig('fig4_janus_contacts.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Saved fig4_janus_contacts.pdf")


if __name__ == '__main__':
    frames = read_trajectory(TRAJ_FILE)
    a = compute_anisotropy(frames)
    j = compute_janus(frames)
    o = parse_rdf(RDF_LIO)
    figure1_msd(a)
    figure2_lio(o)
    figure3_zdensity(j)
    figure4_janus(j)
    print("All 4 PDF figures done.")
