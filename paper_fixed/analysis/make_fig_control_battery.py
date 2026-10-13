#!/usr/bin/env python3
"""
make_fig_control_battery.py -- Figure 8: does the Se-face enrichment survive
its own controls?

Plots the Se:S RDF-coordination ratio (Table 3 of main.tex) for the
published composite and four controls, at each of the four cutoffs used
throughout the paper. A ratio of 1 (dashed line) is the null expectation.
Numbers are taken directly from analyze_faces_3d.py's RDF cross-check,
run identically on all five systems (composite, null_se, swap, rep1, rep2).

USAGE
    python3 make_fig_control_battery.py

OUTPUT
    fig8_control_battery.pdf
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import traj_io as T

plt.rcParams.update(T.PLOT_RC)

CUTOFFS = [4, 5, 6, 8]

# system -> (ratio, sigma) at rc = 4, 5, 6, 8 A
DATA = {
    'Composite\n(published)':   [(1.19, 0.6), (1.52, 3.4), (1.35, 4.2), (1.18, 5.0)],
    'Null\n(faces identical)':  [(0.51, 2.7), (1.06, 0.6), (1.08, 1.1), (1.16, 4.8)],
    'Causal\n(chemistry swapped)': [(2.21, 2.4), (1.55, 4.2), (1.36, 4.9), (1.27, 7.0)],
    'Replica 1\n(repacked)':    [(0.49, 1.9), (0.92, 0.7), (1.20, 3.0), (1.35, 9.9)],
    'Replica 2\n(repacked)':    [(0.56, 1.6), (0.59, 4.0), (0.72, 4.5), (0.94, 1.7)],
}
SYSTEMS = list(DATA.keys())

# fixed categorical order, Okabe-Ito (colorblind-safe)
COLORS = {
    'Composite\n(published)':      '#000000',  # black
    'Null\n(faces identical)':     '#0072B2',  # blue
    'Causal\n(chemistry swapped)': '#D55E00',  # vermillion
    'Replica 1\n(repacked)':       '#009E73',  # bluish green
    'Replica 2\n(repacked)':       '#56B4E9',  # sky blue
}

fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), sharey=True)
y = np.arange(len(SYSTEMS))[::-1]

for ax, rc, idx in zip(axes, CUTOFFS, range(4)):
    ax.axvline(1.0, color='#888888', lw=1.6, ls='--', zorder=1)
    for yi, sysname in zip(y, SYSTEMS):
        ratio, sigma = DATA[sysname][idx]
        c = COLORS[sysname]
        ax.scatter([ratio], [yi], s=90, color=c, zorder=3,
                   edgecolor='white', linewidth=1.0)
        ax.text(ratio, yi + 0.30, f'{sigma:.1f}$\\sigma$',
                ha='center', va='bottom', fontsize=11, color=c)
    ax.set_xlim(0.0, 2.6)
    ax.set_xlabel('Se:S ratio')
    ax.set_title(f'$r_c={rc}$' r'$\ \mathrm{\AA}$', fontsize=15)
    T.apply_style(ax)
    ax.set_yticks(y)
    if idx == 0:
        ax.set_yticklabels(SYSTEMS, fontsize=12)
    else:
        ax.tick_params(axis='y', length=0)
    ax.set_ylim(y.min() - 0.6, y.max() + 0.6)

fig.suptitle('Se-face enrichment ratio: published result vs. null, causal, '
             'and replicate controls', fontsize=15, y=1.04)
fig.text(0.5, -0.04,
          'Dashed line = null expectation (ratio = 1, no face preference). '
          'The null control reproduces the published ratio at $r_c=8$'
          r'$\ \mathrm{\AA}$; the causal control (chemistry swapped) does '
          'not move off the composite; the two replicas disagree on sign '
          'at $r_c\\leq6$' r'$\ \mathrm{\AA}$.',
          ha='center', fontsize=11, wrap=True)

fig.tight_layout(rect=[0, 0.02, 1, 0.96])
fig.savefig('fig8_control_battery.pdf', bbox_inches='tight')
print('Wrote fig8_control_battery.pdf')
