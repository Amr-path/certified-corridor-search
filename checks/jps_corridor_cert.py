"""Test: JPS inside a corridor (outside treated as blocked) + frontier certificate from scan touches.
Compare with A*-in-corridor certificate (same semantics). Uniform grids.
Claims tested:
  T1: OPT >= min(U_J, Phi_J)
  T2: min(U_J, Phi_J) == min(U_A, Phi_A)   (JPS certificate is as tight as the A* one)
  T3: U_J == U_A (both = optimal within corridor-with-outside-blocked)
"""
import math, random, sys, heapq
sys.path.insert(0, '.')
from jps_touch import Grid, dijkstra, jps, octile, SQ2

DIRS = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]

def band(H, W, s, g, k):
    x0,y0 = s; x1,y1 = g; line=[]
    dx, dy = abs(x1-x0), abs(y1-y0); sx = 1 if x0<x1 else -1; sy = 1 if y0<y1 else -1
    err = dx-dy; x,y = x0,y0
    while True:
        line.append((x,y))
        if (x,y)==(x1,y1): break
        e2 = 2*err
        if e2 > -dy: err -= dy; x += sx
        if e2 < dx: err += dx; y += sy
    C=set()
    for (x,y) in line:
        for ax in range(-k,k+1):
            for ay in range(-k,k+1):
                if 0<=x+ax<W and 0<=y+ay<H: C.add((x+ax,y+ay))
    return C

def exit_candidates(G, C, goal):
    """for each corridor cell c, the minimal extra term over moves (c,v) legal in G but not inside the
    corridor-with-outside-blocked graph: d(c,v) + h(v). Returns dict c -> term."""
    h = lambda n: octile(n, goal)
    out = {}
    for (x,y) in C:
        best = math.inf
        for dx,dy in DIRS:
            if not G.can_move(x,y,dx,dy): continue
            v = (x+dx,y+dy)
            inside = v in C and (not (dx and dy) or ((x+dx,y) in C and (x,y+dy) in C))
            if not inside:
                best = min(best, (SQ2 if dx and dy else 1.0) + h(v))
        if best < math.inf: out[(x,y)] = best
    return out

def astar_corridor(Gc, s, goal, C, cand):
    """A* on corridor grid Gc (outside blocked). Returns U, Phi computed from closed/generated g of boundary cells
    at goal close (uses g of cells in OPEN too, as in Lemma 1: shallowest unexpanded node)."""
    h = lambda n: octile(n, goal)
    gval={s:0.0}; closed=set(); OPEN=[(h(s),0.0,s)]
    while OPEN:
        f,g,n = heapq.heappop(OPEN)
        if n in closed or g>gval[n]+1e-12: continue
        closed.add(n)
        if n==goal:
            U=g; break
        x,y=n
        for dx,dy in DIRS:
            if Gc.can_move(x,y,dx,dy):
                m=(x+dx,y+dy); ng=g+(SQ2 if dx and dy else 1.0)
                if ng<gval.get(m,math.inf)-1e-12:
                    gval[m]=ng; heapq.heappush(OPEN,(ng+h(m),ng,m))
    else:
        U=math.inf
    Phi = min((gval[c]+cand[c] for c in cand if c in gval), default=math.inf)
    return U, Phi

def main(seed, trials, H, W):
    rng=random.Random(seed); t1=t2=t3=0; tested=0; certs=0
    for tr in range(trials):
        pb=rng.choice([0.1,0.2,0.3,0.4])
        cells=[[rng.random()>pb for _ in range(W)] for _ in range(H)]
        s=(rng.randrange(W),rng.randrange(H)); g=(rng.randrange(W),rng.randrange(H))
        cells[s[1]][s[0]]=True; cells[g[1]][g[0]]=True
        if s==g: continue
        G=Grid(cells); dist=dijkstra(G,s)
        if g not in dist: continue
        OPT=dist[g]; tested+=1
        for k in [0,1,2,4,8]:
            C=band(H,W,s,g,k)
            cc=[[cells[y][x] and (x,y) in C for x in range(W)] for y in range(H)]
            Gc=Grid(cc)
            cand=exit_candidates(G,C,g)
            touch={s:0.0}
            UJ=jps(Gc,s,g,touch)
            PhiJ=min((touch[c]+cand[c] for c in cand if c in touch), default=math.inf)
            UA,PhiA=astar_corridor(Gc,s,g,C,cand)
            if OPT < min(UJ,PhiJ)-1e-9: t1+=1; print("T1 FAIL",tr,k,OPT,UJ,PhiJ)
            if abs(min(UJ,PhiJ)-min(UA,PhiA))>1e-9: t2+=1; print("T2 FAIL",tr,k,"J:",UJ,PhiJ,"A:",UA,PhiA)
            if abs(UJ-UA)>1e-9 and not (UJ==math.inf and UA==math.inf): t3+=1; print("T3 FAIL",tr,k,UJ,UA)
            if UJ<math.inf and UJ<=1.2*PhiJ+1e-12: certs+=1
    print(f"seed={seed} tested={tested} T1fail={t1} T2fail={t2} T3fail={t3} certified(w=1.2) over (inst,k)={certs}")

if __name__=="__main__":
    main(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]))
