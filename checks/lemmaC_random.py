import os
import random, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lemmaC_exhaustive import check
from jps_touch import Grid, dijkstra
seed, trials, H, W = map(int, sys.argv[1:5])
rng = random.Random(seed); fails = {}; n=0; sp_goal=0; sp_inst=0
for tr in range(trials):
    pb = rng.choice([0.1,0.2,0.3,0.4,0.5])
    cells = [[rng.random() > pb for _ in range(W)] for _ in range(H)]
    s=(rng.randrange(W),rng.randrange(H)); t=(rng.randrange(W),rng.randrange(H))
    cells[s[1]][s[0]]=True; cells[t[1]][t[0]]=True
    if s==t: continue
    G=Grid(cells); dist=dijkstra(G,s)
    for mode in ('goal','exh'):
        if mode=='goal' and t not in dist: continue
        r = check(G,s,t,dist,mode=='exh'); n+=1
        if r=="NONOPT" or r[0]!="OK":
            k = r if r=="NONOPT" else r[0]; fails[k]=fails.get(k,0)+1
            if sum(fails.values())<=3:
                print("FAIL",mode,r,"s=",s,"t=",t); [print(''.join('.' if c else '#' for c in row)) for row in cells]
        elif mode=='goal' and r[1]>0: sp_goal+=r[1]; sp_inst+=1
print(f"seed={seed} {H}x{W} runs={n} fails={fails} goal-mode suboptimal-g pops={sp_goal} in {sp_inst} instances")
