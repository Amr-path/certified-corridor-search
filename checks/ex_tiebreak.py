import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phi_random_corr3 import *
rows=["Tc.cc","wcccc",".cccS","ccccc","ccwcc"]
H,W=5,5; t=[[1]*W for _ in range(H)]
C=set();Cp=set()
for y,r in enumerate(rows):
    for x,ch in enumerate(r):
        if ch in 'cST': C.add((x,y)); Cp.add((x,y))
        if ch=='w': Cp.add((x,y))
s=(4,2); g=(0,0)
# static id order: (3,1) last, everything else row-major
order=[(x,y) for y in range(H) for x in range(W) if (x,y)!=(3,1)]+[(3,1)]
rng=random.Random(0); rng.perm={n:i for i,n in enumerate(order)}
for name,CC in [('C',C),("C'",Cp)]:
    for tb in ['rid','lowg']:
        U,Phi,EXIT,cl,gv,pa=astar_gc(t,s,g,CC,tb,rng)
        print(f"{name:3s} tiebreak={tb:5s}: U={U.val():.4f} Phi={Phi.val():.4f} (3,1) closed? {(3,1) in cl}; EXIT={{{', '.join(f'{b}:{v.val():.3f}' for b,v in EXIT.items())}}}")
