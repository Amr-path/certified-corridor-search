#!/usr/bin/env python3
"""Baseline from the authors' code: focal-list JPS variants of Carlson, Harabor & Stuckey (2024), run with the
warthog-fls binary on the same every-25th instance sample. Writes one TSV per (set, alg, w, scen).
Usage: python3 run_wfls.py"""
import os, sys, subprocess, glob, time, concurrent.futures as cf
ROOT = os.environ.get('CCS_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root; override with CCS_ROOT
BIN = os.environ.get('WARTHOG_BIN', f'{ROOT}/warthog-fls/warthog/build/fast/bin/warthog')
BENCH = f'{ROOT}/bench'; OUT = f'{ROOT}/results_wfls'; SUB = f'{ROOT}/bench_sub25'
SETS = ['dao', 'bg512', 'sc1', 'random', 'room', 'street', 'maze']
WS = ['1.05', '1.2', '1.5', '2.0']; EVERY = 25

def subsample(setname, scen):
    d = f'{SUB}/{setname}'; os.makedirs(d, exist_ok=True)
    out = f'{d}/{os.path.basename(scen)}'
    if not os.path.exists(out):
        lines = open(scen, newline='').read().replace('\r\n', '\n').split('\n')
        head, rows = lines[0], [l for l in lines[1:] if l.strip()]
        keep = [rows[i] for i in range(0, len(rows), EVERY)]
        open(out, 'w').write(head + '\n' + '\n'.join(keep) + '\n')
    return out

def job(args):
    setname, scen, alg, w = args
    outdir = f'{OUT}/{setname}'; os.makedirs(outdir, exist_ok=True)
    tag = f'{alg}_w{w}' if w else alg
    outf = f'{outdir}/{tag}__{os.path.basename(scen)}.tsv'
    if os.path.exists(outf) and os.path.getsize(outf) > 0: return outf, 'skip'
    sub = subsample(setname, scen)
    cmd = [BIN, '--alg', alg, '--scen', sub] + (['--w', w] if w else [])
    r = subprocess.run(cmd, cwd=f'{BENCH}/{setname}', stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0: return outf, 'FAIL ' + r.stderr[:200]
    open(outf + '.tmp', 'w').write(r.stdout); os.rename(outf + '.tmp', outf)
    return outf, 'ok'

def main():
    jobs = []
    for s in SETS:
        for scen in sorted(glob.glob(f'{BENCH}/{s}/*.scen')):
            jobs.append((s, scen, 'jps2', None))
            for w in WS:
                jobs.append((s, scen, 'cjps-fl', w)); jobs.append((s, scen, 'bjps2-fl', w))
    print(f'{len(jobs)} jobs', flush=True); t0 = time.time(); done = fails = 0
    with cf.ThreadPoolExecutor(2) as ex:
        for outf, st in ex.map(job, jobs):
            done += 1
            if st.startswith('FAIL'): fails += 1; print('FAIL', outf, st, flush=True)
            if done % 200 == 0: print(f'{done}/{len(jobs)} {fails} fails {time.time()-t0:.0f}s', flush=True)
    print(f'finished {done} jobs, {fails} fails, {time.time()-t0:.0f}s', flush=True)
if __name__ == '__main__': main()
