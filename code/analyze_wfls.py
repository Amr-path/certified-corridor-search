#!/usr/bin/env python3
"""Aggregate the warthog-fls baseline (focal-list JPS variants of Carlson et al. 2024).
Per set and w: gm and p99 of pcost/opt, max, bound violations, mean expansions relative to warthog jps2 on the same instances,
plus a cross-check of jps2 against the scenario optimum and against our own JPS expansions.
Usage: python3 analyze_wfls.py <outjson>"""
import os, csv, glob, sys, json, math, statistics as stt
ROOT = os.environ.get('CCS_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); R = f'{ROOT}/results_wfls'; SUB = f'{ROOT}/bench_sub25'; OURS = f'{ROOT}/results/uniform'
SETS = ['bg512', 'dao', 'sc1', 'street', 'random', 'room', 'maze']; WS = ['1.05', '1.2', '1.5', '2.0']
def gmean(xs): xs = [x for x in xs if x > 0]; return math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float('nan')
def pct(xs, p): xs = sorted(xs); return xs[min(len(xs) - 1, int(round(p * (len(xs) - 1))))] if xs else float('nan')
def load(s, tag):
    out = {}
    for f in glob.glob(f'{R}/{s}/{tag}__*.tsv'):
        scen = os.path.basename(f)[len(tag) + 2:-4]
        for r in csv.DictReader(open(f), delimiter='\t'):
            out[(scen, int(r['id']))] = (float(r['pcost']), int(r['expanded']), int(r['generated']))
    return out
def opt_of(s):
    out = {}
    for f in glob.glob(f'{SUB}/{s}/*.scen'):
        scen = os.path.basename(f); rows = [l for l in open(f).read().split('\n')[1:] if l.strip()]
        for i, l in enumerate(rows): out[(scen, i)] = float(l.split('\t')[8])
    return out
def ours_jps(s):
    out = {}
    for f in glob.glob(f'{OURS}/{s}/jps__*.csv'):
        scen = os.path.basename(f).split('__', 1)[1][:-4]
        for r in csv.DictReader(open(f)): out[(scen, int(r['id']) // 25)] = int(r['expanded'])   # our id is the scen row index; sub file row = id/25
    return out
res = {}
for s in SETS:
    opt = opt_of(s); base = load(s, 'jps2'); mine = ours_jps(s)
    keys = [k for k in opt if opt[k] > 0 and k in base]
    chk = sum(1 for k in keys if abs(base[k][0] - opt[k]) > 1e-4)
    common = [k for k in keys if k in mine]
    res[s] = {'n': len(keys), 'jps2_not_opt': chk, 'jps2_exp': stt.mean(base[k][1] for k in keys),
              'ours_jps_exp_on_common': stt.mean(mine[k] for k in common), 'jps2_exp_on_common': stt.mean(base[k][1] for k in common), 'n_common': len(common)}
    for alg in ('cjps-fl', 'bjps2-fl'):
        for w in WS:
            D = load(s, f'{alg}_w{w}'); kk = [k for k in keys if k in D]
            tf = [D[k][0] / opt[k] for k in kk]
            res[s][f'{alg}_w{w}'] = {'n': len(kk), 'boundviol': sum(1 for x in tf if x > float(w) + 1e-6), 'gm': gmean(tf), 'p99': pct(tf, .99), 'max': max(tf),
                                     'exp_rel_jps2': stt.mean(D[k][1] for k in kk) / res[s]['jps2_exp'], 'unsolved': sum(1 for k in keys if k not in D or D[k][0] <= 0)}
json.dump(res, open(sys.argv[1], 'w'), indent=1)
for s in SETS:
    r = res[s]; print(s, 'n', r['n'], 'jps2!=opt', r['jps2_not_opt'], 'jps2 exp %.1f ours %.1f (n=%d)' % (r['jps2_exp_on_common'], r['ours_jps_exp_on_common'], r['n_common']))
    for alg in ('cjps-fl', 'bjps2-fl'):
        print('  ', alg, ' | '.join(f"w{w}: gm {r[alg+'_w'+w]['gm']:.3f} p99 {r[alg+'_w'+w]['p99']:.3f} max {r[alg+'_w'+w]['max']:.2f} exp/jps2 {r[alg+'_w'+w]['exp_rel_jps2']:.2f} bv {r[alg+'_w'+w]['boundviol']} uns {r[alg+'_w'+w]['unsolved']}" for w in WS))
