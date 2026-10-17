"""Checks for C1 v2 open items.
 (A) Goal-close variant of corridor A* on weighted grids: is Phi (min f over EXIT at goal-close) monotone
     nondecreasing when the corridor widens?  Also L = min(U,Phi) monotone?
 (B) Prop 2: uniform grid, s=(0,0), t=(l,0), axis band |y|<=k: min over cells b with |y|=k+1 of oct(s,b)+oct(b,t)
     equals oct((l,2(k+1))) = max(l,2(k+1)) + (sqrt2-1) min(l,2(k+1))?  (pure geometry, brute force over b)
"""
import math, random, sys, heapq
sys.path.insert(0, '.')
from check_certificate import make_grid, neighbours, tmin, octile, dijkstra, bresenham_band, SQ2

def corridor_astar_goalclose(t, s, g, C):
    k = tmin(t); h = lambda n: k*octile(n, g)
    gval = {s: 0.0}; closed = set(); OPEN = [(h(s), 0.0, s)]; EXIT = {}
    while OPEN:
        f, gn, n = heapq.heappop(OPEN)
        if n in closed or gn > gval[n] + 1e-12: continue
        closed.add(n)
        if n == g:
            Phi = min(EXIT.values()) if EXIT else math.inf
            return gn, Phi
        for m, c in neighbours(t, n):
            ng = gn + c
            if ng < gval.get(m, math.inf) - 1e-12:
                gval[m] = ng; fm = ng + h(m)
                if m in C: heapq.heappush(OPEN, (fm, ng, m))
                else: EXIT[m] = fm
    Phi = min(EXIT.values()) if EXIT else math.inf
    return math.inf, Phi

def partA(seed, trials, H, W):
    rng = random.Random(seed); phi_dec = 0; L_dec = 0; tested = 0; viol = 0
    for tr in range(trials):
        pb = rng.choice([0.1,0.2,0.3,0.4])
        t = make_grid(H, W, pb, rng, nterr=rng.choice([1,2,4,8]))
        s=(0,0); g=(W-1,H-1)
        if t[0][0]==math.inf: t[0][0]=1.0
        if t[H-1][W-1]==math.inf: t[H-1][W-1]=1.0
        dist = dijkstra(t, s)
        if g not in dist: continue
        OPT = dist[g]; tested += 1
        prevPhi = None; prevL = None
        for k in [0,1,2,3,5,8,13,21,34]:
            C = bresenham_band(H, W, s, g, k)
            U, Phi = corridor_astar_goalclose(t, s, g, C)
            L = min(U, Phi)
            if OPT < L - 1e-9: viol += 1; print("THM VIOLATION", tr, k, OPT, U, Phi)
            if prevPhi is not None and Phi < prevPhi - 1e-9:
                phi_dec += 1
                if phi_dec <= 3: print(f"Phi decreased: trial {tr} k={k} prevPhi={prevPhi:.4f} Phi={Phi:.4f} U={U:.4f} OPT={OPT:.4f}")
            if prevL is not None and L < prevL - 1e-9:
                L_dec += 1
                if L_dec <= 3: print(f"L decreased: trial {tr} k={k} prevL={prevL:.4f} L={L:.4f}")
            prevPhi, prevL = Phi, L
    print(f"[A] seed={seed} tested={tested} theorem violations={viol} Phi decreases={phi_dec} L decreases={L_dec}")

def partB():
    bad = 0
    for l in range(1, 60):
        for k in range(0, 30):
            best = math.inf
            for x in range(-40, l+41):
                for y in (k+1, -(k+1)):
                    best = min(best, octile((0,0),(x,y)) + octile((x,y),(l,0)))
            pred = max(l, 2*(k+1)) + (SQ2-1)*min(l, 2*(k+1))  # = oct((l, 2(k+1)))
            if abs(best - pred) > 1e-9:
                bad += 1
                if bad <= 5: print(f"[B] l={l} k={k} brute={best:.5f} formula={pred:.5f}")
    print(f"[B] mismatches={bad}")

if __name__ == "__main__":
    partA(int(sys.argv[1]), int(sys.argv[2]), 40, 40)
    partB()
