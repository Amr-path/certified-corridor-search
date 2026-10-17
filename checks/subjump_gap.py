import os
"""Count cells u with f*(u)<C* whose FIRST optimal-cost visit is the terminal cell of an orthogonal sub-jump
(inside a diagonal scan) that returned u because u has a forced successor -- u is NOT enqueued by that scan
(the diagonal cell is). The draft's Lemma C case split (a)/(b) does not cover this situation."""
import heapq, math, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jps_touch import Grid, dijkstra, octile, SQ2

def run(G, s, goal):
    h=lambda n: octile(n,goal); gval={s:0.0}; parent={s:None}; closed=set(); OPEN=[(h(s),0.0,s)]
    visit={}  # cell -> (g, kind) kind in {'top','sub-pass','sub-stop'}
    def note(c,g,kind):
        if g < visit.get(c,(math.inf,))[0]-1e-12: visit[c]=(g,kind)
    def jump(x,y,dx,dy,g,sub):
        while True:
            if not G.can_move(x,y,dx,dy): return None
            x,y=x+dx,y+dy; g+= SQ2 if dx and dy else 1.0
            note((x,y),g,'sub-pass' if sub else 'top')
            if (x,y)==goal: return (x,y),g
            if dx and dy:
                if jump(x,y,dx,0,g,True) is not None or jump(x,y,0,dy,g,True) is not None: return (x,y),g
            else:
                stop = ((not G.free(x-dx,y+1) and G.free(x,y+1)) or (not G.free(x-dx,y-1) and G.free(x,y-1))) if dx else \
                       ((not G.free(x+1,y-dy) and G.free(x+1,y)) or (not G.free(x-1,y-dy) and G.free(x-1,y)))
                if stop:
                    if sub: visit[(x,y)] = (visit[(x,y)][0], 'sub-stop') if abs(visit[(x,y)][0]-g)<1e-12 else visit[(x,y)]
                    return (x,y),g
    def pruned_dirs(n,p):
        x,y=n
        if p is None: return [(dx,dy) for dx in (-1,0,1) for dy in (-1,0,1) if dx or dy]
        dx=(x-p[0]);dy=(y-p[1]);dx=(dx>0)-(dx<0);dy=(dy>0)-(dy<0);dirs=[]
        if dx and dy: dirs+=[(dx,0),(0,dy),(dx,dy)]
        elif dx:
            dirs.append((dx,0))
            if not G.free(x-dx,y+1) and G.free(x,y+1): dirs+=[(0,1),(dx,1)]
            if not G.free(x-dx,y-1) and G.free(x,y-1): dirs+=[(0,-1),(dx,-1)]
        else:
            dirs.append((0,dy))
            if not G.free(x+1,y-dy) and G.free(x+1,y): dirs+=[(1,0),(1,dy)]
            if not G.free(x-1,y-dy) and G.free(x-1,y): dirs+=[(-1,0),(-1,dy)]
        return dirs
    while OPEN:
        f,g,n=heapq.heappop(OPEN)
        if n in closed or g>gval[n]+1e-12: continue
        closed.add(n)
        if n==goal: return g,visit
        for d in pruned_dirs(n,parent[n]):
            r=jump(n[0],n[1],d[0],d[1],g,False)
            if r is None: continue
            m,gm=r
            if gm<gval.get(m,math.inf)-1e-12: gval[m]=gm;parent[m]=n;heapq.heappush(OPEN,(gm+h(m),gm,m))
    return math.inf,visit

rng=random.Random(11); cnt=0; inst=0; shown=0
for tr in range(3000):
    H=W=rng.choice([6,8,10]); pb=rng.choice([0.2,0.3,0.4])
    cells=[[rng.random()>pb for _ in range(W)] for _ in range(H)]
    s=(rng.randrange(W),rng.randrange(H)); t=(rng.randrange(W),rng.randrange(H))
    cells[s[1]][s[0]]=True; cells[t[1]][t[0]]=True
    if s==t: continue
    G=Grid(cells); dist=dijkstra(G,s)
    if t not in dist: continue
    C,visit=run(G,s,t); hit=[]
    for u,gu in dist.items():
        if u!=s and gu+octile(u,t)<C-1e-9 and u in visit and abs(visit[u][0]-gu)<1e-9 and visit[u][1]=='sub-stop':
            hit.append(u)
    if hit:
        cnt+=len(hit); inst+=1
        if shown<1 and H==6:
            shown+=1; print("example: s=",s,"t=",t,"C*=%.3f"%C,"cells u (optimal first visit = terminal cell of a sub-jump, not enqueued):",hit)
            for y in range(H): print(''.join('S' if (x,y)==s else 'T' if (x,y)==t else 'U' if (x,y) in hit else ('.' if cells[y][x] else '#') for x in range(W)))
print("cells falling outside draft's (a)/(b) dichotomy:",cnt,"in",inst,"instances")
