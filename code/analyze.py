#!/usr/bin/env python3
"""Aggregate per-instance CSVs into the tables of Section 7. Writes LaTeX tables + JSON summary.
Usage: python3 analyze.py <results_root/uniform|weighted> <outdir>
All numbers come from the logged records; nothing is estimated."""
import csv, glob, os, sys, json, math, statistics as stt
from collections import defaultdict

root, outdir = sys.argv[1], sys.argv[2]
os.makedirs(outdir, exist_ok=True)
INF = float('inf')
sets = sorted(d for d in os.listdir(root) if os.path.isdir(f'{root}/{d}'))

def load(setname, cfg):
    rows = {}
    for f in glob.glob(f'{root}/{setname}/{cfg}__*.csv'):
        scen = os.path.basename(f).split('__', 1)[1]
        for r in csv.DictReader(open(f)):
            key = (scen, int(r['id']))
            rows[key] = {k: (float(v) if k in ('micros','U','Phi','L','opt') else v) for k, v in r.items()}
            for k in ('expanded','generated','scanned','rounds'): rows[key][k] = int(r[k])
    return rows

def gmean(xs):
    xs = [x for x in xs if x > 0 and x < INF]
    return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float('nan')
def pct(xs, p):
    xs = sorted(x for x in xs if x < INF)
    if not xs: return float('nan')
    return xs[min(len(xs)-1, int(round(p * (len(xs)-1))))]

WS = ['1.05', '1.2', '1.5', '2.0']; KS = ['4', '16', '64']
summary = {}
T1 = []; T2 = []; T3 = []; T4 = []
for s in sets:
    base = load(s, 'astar')
    if not base: continue
    jps = load(s, 'jps')
    # reference optimal cost = astar U (also check against scen opt where available)
    opt = {k: r['U'] for k, r in base.items() if 0 < r['U'] < INF}
    n = len(opt)
    summary[s] = {'instances': n, 'astar_exp': stt.mean(base[k]['expanded'] for k in opt), 'astar_us': stt.mean(base[k]['micros'] for k in opt)}
    if jps: summary[s]['jps_exp'] = stt.mean(jps[k]['expanded'] for k in opt if k in jps); summary[s]['jps_us'] = stt.mean(jps[k]['micros'] for k in opt if k in jps)
    # ---- Table 1/2: single-band certification and factors
    for alg in (['castar', 'cjps'] if jps else ['castar']):
        for k in KS:
            R = load(s, f'{alg}_k{k}')
            if not R: continue
            keys = [q for q in opt if q in R]
            solved = [q for q in keys if R[q]['U'] < INF]
            viol = sum(1 for q in keys if opt[q] < R[q]['L'] - 1e-6)
            certw = {w: sum(1 for q in solved if R[q]['U'] <= float(w) * R[q]['Phi'] + 1e-9) / max(1, len(keys)) for w in WS}
            cf = [R[q]['U'] / R[q]['L'] for q in solved]; tf = [R[q]['U'] / opt[q] for q in solved]
            T1.append((s, alg, k, len(keys), len(solved)/len(keys), certw, viol))
            T2.append((s, alg, k, gmean(cf), pct(cf, .5), pct(cf, .99), gmean(tf), pct(tf, .5), pct(tf, .99),
                       stt.mean(R[q]['expanded'] for q in solved) / stt.mean(base[q]['expanded'] for q in solved), stt.mean(R[q]['micros'] for q in solved)))
    # ---- Table 3: bounded-suboptimal comparison at matched w
    for w in WS:
        row = {'set': s, 'w': w}
        for alg in (['wastar', 'fls', 'ccsA', 'ccsJ', 'ccsAr', 'ccsJr'] if jps else ['wastar', 'fls', 'ccsA', 'ccsAr']):
            R = load(s, f'{alg}_w{w}')
            if not R: continue
            keys = [q for q in opt if q in R and R[q]['U'] < INF]
            bv = sum(1 for q in keys if R[q]['U'] > float(w) * opt[q] + 1e-6)
            tf = [R[q]['U'] / opt[q] for q in keys]
            row[alg] = {'n': len(keys), 'boundviol': bv, 'gm_factor': gmean(tf), 'p99_factor': pct(tf, .99), 'max_factor': max(tf) if tf else float('nan'), 'frac_subopt': sum(1 for x in tf if x > 1 + 1e-9) / max(1, len(tf)),
                        'exp_rel_astar': stt.mean(R[q]['expanded'] for q in keys) / summary[s]['astar_exp'],
                        'us_rel_astar': stt.mean(R[q]['micros'] for q in keys) / summary[s]['astar_us'],
                        'mean_rounds': stt.mean(R[q]['rounds'] for q in keys),
                        'cert_gm': gmean([R[q]['U'] / R[q]['L'] for q in keys]) if alg.startswith('ccs') else None}
            if jps and alg in ('ccsJ','ccsJr'):
                row[alg]['exp_rel_jps'] = stt.mean(R[q]['expanded'] for q in keys) / summary[s]['jps_exp']
                row[alg]['us_rel_jps'] = stt.mean(R[q]['micros'] for q in keys) / summary[s]['jps_us']
        T3.append(row)
    # ---- Table 4: cjps vs castar paired (same band), tightness check
    if jps:
        for k in KS:
            A = load(s, f'castar_k{k}'); J = load(s, f'cjps_k{k}')
            keys = [q for q in opt if q in A and q in J]
            if not keys: continue
            # note: castar is on G_C^+, cjps on M_C -> L may legitimately differ; report both and the fraction equal
            eqL = sum(1 for q in keys if abs(A[q]['L'] - J[q]['L']) < 1e-6) / len(keys)
            T4.append((s, k, len(keys), stt.mean(A[q]['expanded'] for q in keys), stt.mean(J[q]['expanded'] for q in keys),
                       stt.mean(J[q]['scanned'] for q in keys), stt.mean(A[q]['micros'] for q in keys), stt.mean(J[q]['micros'] for q in keys), eqL))

T4mc = []
def load_pair(setname, k):
    A = {}; J = {}
    for f in glob.glob(f'{root}/{setname}/pair_k{k}__*.csv'):
        scen = os.path.basename(f).split('__', 1)[1]
        for r in csv.DictReader(open(f)):
            d = {kk: (float(v) if kk in ('micros','U','Phi','L','opt') else v) for kk, v in r.items()}
            for kk in ('expanded','generated','scanned'): d[kk] = int(r[kk])
            (A if r['alg'] == 'pmc' else J)[(scen, int(r['id']))] = d
    return A, J
for s in sets:
    base = load(s, 'astar')
    opt = {k: r['U'] for k, r in base.items() if 0 < r['U'] < INF}
    for k in KS:
        A, J = load_pair(s, k)
        keys = [q for q in opt if q in A and q in J]
        if not keys: continue
        eqL = sum(1 for q in keys if abs(A[q]['L'] - J[q]['L']) < 1e-6 or (A[q]['L']==INF and J[q]['L']==INF)) / len(keys)
        ratio_mean = stt.mean(A[q]['micros'] for q in keys) / stt.mean(J[q]['micros'] for q in keys)
        per_inst = stt.mean(A[q]['micros'] / max(1e-9, J[q]['micros']) for q in keys)
        T4mc.append((s, k, len(keys), stt.mean(A[q]['expanded'] for q in keys), stt.mean(J[q]['expanded'] for q in keys),
                     stt.mean(J[q]['scanned'] for q in keys), stt.mean(A[q]['micros'] for q in keys), stt.mean(J[q]['micros'] for q in keys), eqL,
                     stt.mean(A[q]['generated'] for q in keys), stt.mean(J[q]['generated'] for q in keys), ratio_mean, per_inst))
json.dump({'summary': summary, 'T1': T1, 'T2': T2, 'T3': T3, 'T4': T4, 'T4mc': T4mc}, open(f'{outdir}/summary.json', 'w'), indent=1, default=str)

def fmt(x, d=3):
    return '--' if (isinstance(x, float) and (math.isnan(x) or x == INF)) else (f'{x:.{d}f}' if isinstance(x, float) else str(x))

with open(f'{outdir}/table1_certification.tex', 'w') as f:
    f.write('\\begin{tabular}{llr r rrrr r}\\toprule\nSet & Inner & $k$ & solved & $w{=}1.05$ & $1.2$ & $1.5$ & $2$ & viol.\\\\\\midrule\n')
    for (s, alg, k, n, solved, certw, viol) in T1:
        f.write(f"{s} & {alg} & {k} & {fmt(solved,2)} & " + ' & '.join(fmt(certw[w],2) for w in WS) + f" & {viol}\\\\\n")
    f.write('\\bottomrule\\end{tabular}\n')
with open(f'{outdir}/table2_factors.tex', 'w') as f:
    f.write('\\begin{tabular}{llr rrr rrr rr}\\toprule\n & & & \\multicolumn{3}{c}{certified $U/L$} & \\multicolumn{3}{c}{true $U/C^*$} & & \\\\\nSet & Inner & $k$ & geomean & median & p99 & geomean & median & p99 & exp. & $\\mu$s\\\\\\midrule\n')
    for row in T2:
        s, alg, k = row[:3]; f.write(f"{s} & {alg} & {k} & " + ' & '.join(fmt(x) for x in row[3:9]) + f" & {fmt(row[9],0)} & {fmt(row[10],0)}\\\\\n")
    f.write('\\bottomrule\\end{tabular}\n')
with open(f'{outdir}/table3_bounded.tex', 'w') as f:
    f.write('\\begin{tabular}{ll rrr rrr rrr rrrr}\\toprule\n & & \\multicolumn{3}{c}{WA*} & \\multicolumn{3}{c}{FLS} & \\multicolumn{3}{c}{CCS-A*} & \\multicolumn{4}{c}{CCS-JPS}\\\\\nSet & $w$ & gm & p99 & exp & gm & p99 & exp & gm & p99 & exp & gm & p99 & exp & exp/JPS\\\\\\midrule\n')
    for row in T3:
        cells = [row['set'], row['w']]
        for alg in ['wastar', 'fls', 'ccsA']:
            r = row.get(alg); cells += [fmt(r['gm_factor']), fmt(r['p99_factor']), fmt(r['exp_rel_astar'], 2)] if r else ['--']*3
        r = row.get('ccsJ'); cells += [fmt(r['gm_factor']), fmt(r['p99_factor']), fmt(r['exp_rel_astar'], 2), fmt(r['exp_rel_jps'], 2)] if r else ['--']*4
        f.write(' & '.join(cells) + '\\\\\n')
    f.write('\\bottomrule\\end{tabular}\n')
with open(f'{outdir}/table4_paired.tex', 'w') as f:
    f.write('\\begin{tabular}{lr r rr r rr r}\\toprule\nSet & $k$ & $n$ & A* exp. & JPS exp. & JPS scanned & A* $\\mu$s & JPS $\\mu$s & $L$ equal\\\\\\midrule\n')
    for row in T4:
        f.write(f"{row[0]} & {row[1]} & {row[2]} & {fmt(row[3],0)} & {fmt(row[4],0)} & {fmt(row[5],0)} & {fmt(row[6],0)} & {fmt(row[7],0)} & {fmt(row[8],3)}\\\\\n")
    f.write('\\bottomrule\\end{tabular}\n')
print(json.dumps(summary, indent=1, default=str))
