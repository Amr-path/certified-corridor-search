#!/usr/bin/env python3
"""Experiment driver. Runs every (set, scen, config) and writes one CSV per (set, config).
Sampling: every EVERY-th instance of each scenario file (instances are sorted by bucket, so this spreads over difficulty).
Usage: python3 run_all.py <phase>   phase in {uniform, weighted, all}
"""
import os, sys, subprocess, glob, time, concurrent.futures as cf

ROOT = os.environ.get('CCS_ROOT', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root; override with CCS_ROOT
BIN = f'{ROOT}/ccs/ccs'
BENCH = f'{ROOT}/bench'
OUT = f'{ROOT}/results'
EVERY = 25
WORKERS = 2

UNI_SETS = ['dao', 'bg512', 'sc1', 'random', 'room', 'street', 'maze']
WEI_SETS = ['dao', 'sc1', 'wc3maps512', 'weighted']
WS = [1.05, 1.2, 1.5, 2.0]
KS = [4, 16, 64]

def configs(weighted):
    c = []
    c.append(('astar', ['astar']))
    if not weighted: c.append(('jps', ['jps']))
    for w in WS:
        c.append((f'wastar_w{w}', ['wastar', '--w', str(w)]))
        c.append((f'fls_w{w}', ['fls', '--w', str(w)]))
    for k in KS:
        c.append((f'castar_k{k}', ['castar', '--k', str(k)]))
        if not weighted: c.append((f'pair_k{k}', ['pair', '--k', str(k)]))
        if not weighted: c.append((f'cjps_k{k}', ['cjps', '--k', str(k)]))
    for w in WS:
        c.append((f'ccsA_w{w}', ['ccs', '--w', str(w), '--inner', 'castar']))
        if not weighted: c.append((f'ccsJ_w{w}', ['ccs', '--w', str(w), '--inner', 'cjps']))
    for w in WS:
        c.append((f'ccsAr_w{w}', ['ccs', '--w', str(w), '--inner', 'castar', '--k0rel', '0.25']))
        if not weighted: c.append((f'ccsJr_w{w}', ['ccs', '--w', str(w), '--inner', 'cjps', '--k0rel', '0.25']))
    return c

def job(args):
    setname, scen, cname, cargs, weighted = args
    outdir = f'{OUT}/{"weighted" if weighted else "uniform"}/{setname}'
    os.makedirs(outdir, exist_ok=True)
    outf = f'{outdir}/{cname}__{os.path.basename(scen)}.csv'
    if os.path.exists(outf) and os.path.getsize(outf) > 0: return outf, 'skip'
    cmd = [BIN] + cargs[:1] + [scen, f'{BENCH}/{setname}'] + cargs[1:] + ['--every', str(EVERY)] + ((['--letters'] if setname == 'weighted' else ['--weighted']) if weighted else [])
    with open(outf + '.tmp', 'w') as f:
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        os.remove(outf + '.tmp'); return outf, 'FAIL ' + r.stderr[:200]
    os.rename(outf + '.tmp', outf)
    return outf, 'ok'

def main():
    phase = sys.argv[1] if len(sys.argv) > 1 else 'all'
    jobs = []
    if phase in ('uniform', 'all'):
        for s in UNI_SETS:
            for scen in sorted(glob.glob(f'{BENCH}/{s}/*.scen')):
                for cname, cargs in configs(False): jobs.append((s, scen, cname, cargs, False))
    if phase in ('weighted', 'all'):
        for s in WEI_SETS:
            for scen in sorted(glob.glob(f'{BENCH}/{s}/*.scen')):
                for cname, cargs in configs(True): jobs.append((s, scen, cname, cargs, True))
    print(f'{len(jobs)} jobs', flush=True)
    t0 = time.time(); done = 0; fails = 0
    with cf.ThreadPoolExecutor(WORKERS) as ex:
        for outf, status in ex.map(job, jobs):
            done += 1
            if status.startswith('FAIL'): fails += 1; print('FAIL', outf, status, flush=True)
            if done % 200 == 0: print(f'{done}/{len(jobs)} done, {fails} fails, {time.time()-t0:.0f}s', flush=True)
    print(f'finished {done} jobs, {fails} fails, {time.time()-t0:.0f}s', flush=True)

if __name__ == '__main__':
    main()
