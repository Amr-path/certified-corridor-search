import os
import random, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lemmaC_exhaustive import jps_instr
from jps_touch import Grid, dijkstra, octile
rng = random.Random(7)
for tr in range(20000):
    H=W=rng.choice([5,6])
    pb = rng.choice([0.2,0.3,0.4])
    cells = [[rng.random() > pb for _ in range(W)] for _ in range(H)]
    s=(rng.randrange(W),rng.randrange(H)); t=(rng.randrange(W),rng.randrange(H))
    cells[s[1]][s[0]]=True; cells[t[1]][t[0]]=True
    if s==t: continue
    G=Grid(cells); dist=dijkstra(G,s)
    if t not in dist: continue
    res, visit, closed, pops = jps_instr(G,s,t,False)
    bad=[(n,g) for (n,g) in pops if g>dist[n]+1e-9]
    if bad:
        n,g=bad[0]
        print("s=",s,"t=",t,"C*=%.3f"%dist[t],"popped",n,"with g=%.3f"%g,"g*=%.3f"%dist[n],"f*=%.3f"%(dist[n]+octile(n,t)))
        print("optimal visit of",n,": min scan g=%.3f, enqueued by that scan? %s"%(visit[n][0], visit[n][2]))
        for y in range(H):
            print(''.join('S' if (x,y)==s else 'T' if (x,y)==t else 'X' if (x,y)==n else ('.' if cells[y][x] else '#') for x in range(W)))
        print("pop order:", [(p, round(gg,3)) for p,gg in pops])
        break
