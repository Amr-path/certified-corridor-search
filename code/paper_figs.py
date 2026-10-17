#!/usr/bin/env python3
"""Figures for Section 7, from the logged per-instance records only.
Usage: python3 paper_figs.py <outdir>
  fig_ecdf.pdf   : ECDF of certified U/L and true U/C* at k=16 (band-contains-path population, as Table 2)
  fig_rounds.pdf : per-instance expansion ratio CCS / whole-map same inner algorithm, by number of rounds (w=1.05, k0=ceil(l/4))
  fig_pareto.pdf : p99 true factor vs expansions relative to whole-map A*, across w, for WA*, FLS, CCS-A* (both k0)
"""
import os, csv, glob, sys, json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

REPO = os.environ.get('CCS_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ROOT = f'{REPO}/results'
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
INF = float('inf')
NAMES = {'dao': 'DAO', 'bg512': 'BG512', 'sc1': 'SC1', 'random': 'Random', 'room': 'Room', 'street': 'Street', 'maze': 'Maze',
         'wc3maps512': 'WC3', 'weighted': 'Weighted'}
UNI = ['bg512', 'dao', 'sc1', 'street', 'random', 'room', 'maze']
WEI = ['dao', 'sc1', 'wc3maps512', 'weighted']

plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'legend.fontsize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'pdf.fonttype': 42, 'font.family': 'serif'})
C1, C2, C3, C4 = '#1f4e79', '#c0504d', '#6b8e23', '#7f7f7f'

def load(ph, s, cfg, alg_filter=None):
    R = {}
    for f in glob.glob(f'{ROOT}/{ph}/{s}/{cfg}__*.csv'):
        scen = os.path.basename(f).split('__', 1)[1]
        for r in csv.DictReader(open(f)):
            if alg_filter and r['alg'] != alg_filter: continue
            d = {k: float(r[k]) for k in ('U', 'Phi', 'L', 'opt', 'micros')}
            for k in ('expanded', 'generated', 'scanned', 'rounds'): d[k] = int(r[k])
            R[(scen, int(r['id']))] = d
    return R

# ---------------- Fig B: ECDF certified vs true, k=16 ----------------
fig, axes = plt.subplots(2, 7, figsize=(7.0, 3.3), sharey=True)
panels = [('uniform', s) for s in UNI] + [('weighted', s) for s in WEI]
for i, ax in enumerate(axes.flat):
    if i >= len(panels): ax.axis('off'); continue
    ph, s = panels[i]
    base = load(ph, s, 'astar'); opt = {k: r['U'] for k, r in base.items() if 0 < r['U'] < INF}
    R = load(ph, s, 'castar_k16')
    keys = [q for q in opt if q in R and R[q]['U'] < INF]
    cert = np.sort([R[q]['U'] / R[q]['L'] for q in keys]); true = np.sort([R[q]['U'] / opt[q] for q in keys])
    y = np.arange(1, len(keys) + 1) / len(keys)
    ax.step(cert, y, color=C1, lw=1.0, where='post', label='certified $U/L$')
    ax.step(true, y, color=C2, lw=1.0, ls='--', where='post', label='true $U/C^*$')
    ax.set_xscale('log'); ax.set_xlim(1, 12); ax.set_ylim(0, 1.02)
    ax.set_xticks([1, 2, 5, 10]); ax.set_xticklabels(['1', '2', '5', '10']); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_title(f"{NAMES[s]} ({'U' if ph=='uniform' else 'W'})\n$n$={len(keys)}", pad=2)
    ax.grid(alpha=.25, lw=.4)
    if i % 7 == 0: ax.set_ylabel('fraction of instances')
axes[0, 0].legend(loc='lower right', frameon=False)
fig.supxlabel('suboptimality factor (log scale)', y=0.02)
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(f'{out}/fig_ecdf.pdf'); plt.close(fig)

# ---------------- Fig C: expansion ratio by rounds ----------------
def rounds_data(ph, sets, cfg, basecfg):
    groups = {r: [] for r in range(1, 6)}; per_set = {}
    for s in sets:
        base = load(ph, s, basecfg); R = load(ph, s, cfg)
        ok = [q for q in base if 0 < base[q]['U'] < INF and q in R and R[q]['U'] < INF]
        per_set[s] = len(ok)
        for q in ok:
            r = min(5, R[q]['rounds']); groups[r].append(R[q]['expanded'] / max(1, base[q]['expanded']))
    return groups, per_set
fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.3))
specs = [('uniform', UNI, 'ccsAr_w1.05', 'astar', 'CCS\u2013A*, uniform (vs. whole-map A*)'),
         ('uniform', UNI, 'ccsJr_w1.05', 'jps', 'CCS\u2013JPS, uniform (vs. whole-map JPS)'),
         ('weighted', WEI, 'ccsAr_w1.05', 'astar', 'CCS\u2013A*, weighted (vs. whole-map A*)')]
stats = {}
for ax, (ph, sets, cfg, basecfg, title) in zip(axes, specs):
    groups, _ = rounds_data(ph, sets, cfg, basecfg)
    data = [groups[r] for r in range(1, 6)]
    bp = ax.boxplot(data, whis=(5, 95), showfliers=False, widths=.55, patch_artist=True,
                    medianprops=dict(color=C2, lw=1.2), boxprops=dict(facecolor='#dde6f0', edgecolor=C1, lw=.8),
                    whiskerprops=dict(color=C1, lw=.8), capprops=dict(color=C1, lw=.8))
    ax.axhline(1, color=C4, lw=.7, ls=':')
    ax.set_xticks(range(1, 6)); ax.set_xticklabels(['1', '2', '3', '4', '5+'])
    ax.set_xlabel('corridor rounds'); ax.set_title(title.replace('\\ ', ' '), pad=3); ax.grid(axis='y', alpha=.25, lw=.4)
    top = max(np.percentile(d, 95) for d in data if d)
    ax.set_ylim(0, top * 1.18)
    for i, d in enumerate(data):
        ax.text(i + 1, top * 1.08, f'{len(d)}', ha='center', va='bottom', fontsize=6, color=C4)
    stats[title] = {r: (len(d), float(np.median(d)) if d else None) for r, d in zip(range(1, 6), data)}
axes[0].set_ylabel('expansions / whole-map')
fig.tight_layout(); fig.savefig(f'{out}/fig_rounds.pdf'); plt.close(fig)
json.dump(stats, open(f'{out}/fig_rounds_stats.json', 'w'), indent=1)

# ---------------- Fig D: Pareto, 4 panels ----------------
U = json.load(open(f'{REPO}/analysis_uniform/summary.json'))
W = json.load(open(f'{REPO}/analysis_weighted/summary.json'))
fig, axes = plt.subplots(1, 4, figsize=(7.0, 2.3))
for ax, (D, s, tag) in zip(axes, [(U, 'dao', 'U'), (U, 'random', 'U'), (W, 'dao', 'W'), (W, 'weighted', 'W')]):
    rows = {r['w']: r for r in D['T3'] if r['set'] == s}
    ws = ['1.05', '1.2', '1.5', '2.0']
    for alg, lab, col, mk, ls in [('wastar', 'WA*', C2, 'o', '-'), ('fls', 'FLS', C3, 's', '-'),
                                  ('ccsAr', 'CCS--A*, $k_0{=}\\lceil\\ell/4\\rceil$', C1, '^', '-'),
                                  ('ccsA', 'CCS--A*, $k_0{=}2$', C1, 'v', ':')]:
        xs = [rows[w][alg]['exp_rel_astar'] for w in ws]; ys = [rows[w][alg]['p99_factor'] for w in ws]
        ax.plot(xs, ys, marker=mk, ms=3.5, color=col, ls=ls, lw=.9, label=lab.replace('--', '–'))
        for w, x, y in zip(ws, xs, ys):
            if alg == 'wastar': ax.annotate(w, (x, y), textcoords='offset points', xytext=(3, -7), fontsize=5.5, color=C4)
    ax.axvline(1, color=C4, lw=.7, ls=':'); ax.set_xscale('log')
    ax.set_title(f'{NAMES[s]} ({tag})', pad=3); ax.set_xlabel('expansions / whole-map A*'); ax.grid(alpha=.25, lw=.4)
    ax.set_xticks([0.1, 0.3, 1, 3]); ax.set_xticklabels(['0.1', '0.3', '1', '3']); ax.xaxis.set_minor_formatter(NullFormatter())
axes[0].set_ylabel('p99 of true factor $U/C^*$')
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig(f'{out}/fig_pareto.pdf', bbox_inches='tight'); plt.close(fig)
print('ok')
