#!/usr/bin/env python3
"""Emit the four tables of Section 7 (compact layout) from the analysis JSONs. Numbers only from logs."""
import json, math, sys
import os
REPO = os.environ.get('CCS_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
U = json.load(open(f'{REPO}/analysis_uniform/summary.json'))
WF = json.load(open(f'{REPO}/analysis_uniform/wfls.json'))
W = json.load(open(f'{REPO}/analysis_weighted/summary.json'))
out = sys.argv[1]
NAMES = {'dao': 'DAO', 'bg512': 'BG512', 'sc1': 'SC1', 'random': 'Random', 'room': 'Room', 'street': 'Street', 'maze': 'Maze',
         'wc3maps512': 'WC3', 'weighted': 'Weighted'}
def f(x, d=3): return '--' if (x is None or (isinstance(x, float) and math.isnan(x))) else f'{x:.{d}f}'

# ---- Table 1: certification rate, corridor A*, k in {4,16,64}, w in {1.05, 2}
rows = []
for ph, D in (('uniform', U), ('weighted', W)):
    for s in D['summary']:
        r = {k: next(t for t in D['T1'] if t[0] == s and t[1] == 'castar' and t[2] == k) for k in ('4', '16', '64')}
        rows.append((ph, s, D['summary'][s]['instances'], r))
with open(f'{out}/tab_cert.tex', 'w') as fh:
    fh.write('\\begin{tabular}{llr rr rr rr}\\toprule\n & & & \\multicolumn{2}{c}{$k=4$} & \\multicolumn{2}{c}{$k=16$} & \\multicolumn{2}{c}{$k=64$}\\\\\n')
    fh.write(' & Set & $n$ & $w{=}1.05$ & $w{=}2$ & $w{=}1.05$ & $w{=}2$ & $w{=}1.05$ & $w{=}2$\\\\\\midrule\n')
    for ph, s, n, r in rows:
        fh.write(f"{'U' if ph=='uniform' else 'W'} & {NAMES[s]} & {n} & " + ' & '.join(f"{r[k][5]['1.05']:.2f} & {r[k][5]['2.0']:.2f}" for k in ('4','16','64')) + '\\\\\n')
    fh.write('\\bottomrule\\end{tabular}\n')

# ---- Table 2: certified vs true factor at k=16 (uniform: corridor A*; weighted: corridor A*), plus k=64 for weighted
with open(f'{out}/tab_factor.tex', 'w') as fh:
    fh.write('\\begin{tabular}{llr rr rr r}\\toprule\n & & & \\multicolumn{2}{c}{certified $U/L$} & \\multicolumn{2}{c}{true $U/C^*$} & \\\\\n & Set & $k$ & geo.\\ mean & p99 & geo.\\ mean & p99 & exp.\\ / A*\\\\\\midrule\n')
    for ph, D in (('uniform', U), ('weighted', W)):
        for s in D['summary']:
            for k in (('16',) if ph == 'uniform' else ('16', '64')):
                t = next(t for t in D['T2'] if t[0] == s and t[1] == 'castar' and t[2] == k)
                fh.write(f"{'U' if ph=='uniform' else 'W'} & {NAMES[s]} & {k} & {f(t[3])} & {f(t[5])} & {f(t[6])} & {f(t[8])} & {t[9]:.2f}\\\\\n")
    fh.write('\\bottomrule\\end{tabular}\n')

# ---- Table 3: bounded-suboptimal comparison at w = 1.05 and 1.5, k0 = l/4 variants of CCS
with open(f'{out}/tab_bounded.tex', 'w') as fh:
    fh.write('\\begin{tabular}{llr rrr rrr rrrr rrrr rrr}\\toprule\n & & & \\multicolumn{3}{c}{WA*} & \\multicolumn{3}{c}{FLS} & \\multicolumn{4}{c}{CCS--A*} & \\multicolumn{4}{c}{CCS--JPS} & \\multicolumn{3}{c}{CJPS--FLS}\\\\\n')
    fh.write(' & Set & $w$ & gm & p99 & exp & gm & p99 & exp & gm & p99 & exp & rounds & gm & p99 & exp & exp/JPS & gm & p99 & exp/JPS2\\\\\\midrule\n')
    for ph, D in (('uniform', U), ('weighted', W)):
        for row in D['T3']:
            if row['w'] not in ('1.05', '1.5'): continue
            c = [('U' if ph=='uniform' else 'W'), NAMES[row['set']], row['w']]
            for a in ('wastar', 'fls'):
                r = row[a]; c += [f(r['gm_factor']), f(r['p99_factor']), f(r['exp_rel_astar'], 2)]
            r = row['ccsAr']; c += [f(r['gm_factor']), f(r['p99_factor']), f(r['exp_rel_astar'], 2), f(r['mean_rounds'], 1)]
            r = row.get('ccsJr'); c += [f(r['gm_factor']), f(r['p99_factor']), f(r['exp_rel_astar'], 2), f(r['exp_rel_jps'], 2)] if r else ['--'] * 4
            wf = WF.get(row['set'], {}).get(f"cjps-fl_w{row['w']}") if ph == 'uniform' else None
            c += [f(wf['gm']), f(wf['p99']), f(wf['exp_rel_jps2'], 2)] if wf else ['--'] * 3
            fh.write(' & '.join(c) + '\\\\\n')
    fh.write('\\bottomrule\\end{tabular}\n')

# ---- Table 4: JPS_C vs corridor A* on M_C, paired back-to-back, k=16
with open(f'{out}/tab_paired.tex', 'w') as fh:
    fh.write('\\begin{tabular}{lr rr rr r rr rr r}\\toprule\n & & \\multicolumn{2}{c}{expanded} & \\multicolumn{2}{c}{generated} & & \\multicolumn{2}{c}{time (ms)} & \\multicolumn{2}{c}{time ratio A*/JPS} & \\\\\nSet & $n$ & A*$_{M_C}$ & JPS$_C$ & A*$_{M_C}$ & JPS$_C$ & JPS$_C$ scanned & A*$_{M_C}$ & JPS$_C$ & of means & per inst. & $L$ equal\\\\\\midrule\n')
    for t in U['T4mc']:
        if t[1] != '16': continue
        fh.write(f"{NAMES[t[0]]} & {t[2]} & {t[3]:.0f} & {t[4]:.0f} & {t[9]:.0f} & {t[10]:.0f} & {t[5]:.0f} & {t[6]/1000:.2f} & {t[7]/1000:.2f} & {t[11]:.2f} & {t[12]:.2f} & {t[8]:.3f}\\\\\n")
    fh.write('\\bottomrule\\end{tabular}\n')
print('ok')
