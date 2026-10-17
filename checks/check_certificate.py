"""Independent computational check of Draft C1 (corridor certificate).

Model: weighted 8-connected grid, AAAI-2023 cost model (orthogonal = avg of 2 cells,
diagonal = sqrt2 * avg of 4 touched cells), corner cutting forbidden.
Heuristic h = t_min * octile.

Checks:
 (1) h is consistent on every generated instance (all edges).
 (2) Exhaustive corridor A* (exit cells generated but never expanded):
     OPT >= min(U, Phi).  (Theorem 1a)
 (3) Early-stopping variant: stop once t closed and min f over OPEN+EXIT >= U/w:
     OPT >= min(U, Phi').  (Section 6)
 (4) Monotonicity of U over nested corridors; record whether Phi is monotone.
"""
import heapq, math, random, sys

SQ2 = math.sqrt(2)
DIRS = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]

def octile(a, b):
    dx, dy = abs(a[0]-b[0]), abs(a[1]-b[1])
    return max(dx,dy) + (SQ2-1)*min(dx,dy)

def make_grid(H, W, p_block, rng, nterr=4, maxcost=5.0):
    t = [[rng.uniform(1.0, maxcost) if rng.random() > p_block else math.inf for _ in range(W)] for _ in range(H)]
    if nterr == 1:  # uniform cost
        return [[1.0 if v != math.inf else math.inf for v in row] for row in t]
    if nterr:  # quantise to a few terrain types
        for y in range(H):
            for x in range(W):
                if t[y][x] != math.inf:
                    t[y][x] = 1.0 + round((t[y][x]-1.0)/(maxcost-1.0)*(nterr-1))*(maxcost-1.0)/(nterr-1)
    return t

def edge_cost(t, n, m):
    H, W = len(t), len(t[0])
    (x0,y0),(x1,y1) = n, m
    if not (0 <= x1 < W and 0 <= y1 < H): return math.inf
    if t[y0][x0] == math.inf or t[y1][x1] == math.inf: return math.inf
    dx, dy = x1-x0, y1-y0
    if dx and dy:
        a, b = t[y0][x1], t[y1][x0]
        if a == math.inf or b == math.inf: return math.inf  # corner cutting forbidden
        return SQ2 * (t[y0][x0] + t[y1][x1] + a + b) / 4.0
    return (t[y0][x0] + t[y1][x1]) / 2.0

def neighbours(t, n):
    for dx, dy in DIRS:
        m = (n[0]+dx, n[1]+dy)
        c = edge_cost(t, n, m)
        if c < math.inf:
            yield m, c

def tmin(t):
    return min(v for row in t for v in row if v != math.inf)

def dijkstra(t, s):
    dist = {s: 0.0}
    pq = [(0.0, s)]
    while pq:
        d, n = heapq.heappop(pq)
        if d > dist[n]: continue
        for m, c in neighbours(t, n):
            nd = d + c
            if nd < dist.get(m, math.inf):
                dist[m] = nd
                heapq.heappush(pq, (nd, m))
    return dist

def check_consistency(t, goal):
    k = tmin(t)
    H, W = len(t), len(t[0])
    worst = 0.0
    for y in range(H):
        for x in range(W):
            n = (x,y)
            if t[y][x] == math.inf: continue
            hn = k*octile(n, goal)
            for m, c in neighbours(t, n):
                hm = k*octile(m, goal)
                worst = max(worst, hn - (c + hm))
    return worst  # must be <= 1e-9

def bresenham_band(H, W, s, g, k):
    """Chebyshev band of half-width k around the Bresenham line s->g."""
    x0,y0 = s; x1,y1 = g
    line = []
    dx, dy = abs(x1-x0), abs(y1-y0)
    sx = 1 if x0 < x1 else -1; sy = 1 if y0 < y1 else -1
    err = dx - dy
    x, y = x0, y0
    while True:
        line.append((x,y))
        if (x,y) == (x1,y1): break
        e2 = 2*err
        if e2 > -dy: err -= dy; x += sx
        if e2 < dx: err += dx; y += sy
    C = set()
    for (x,y) in line:
        for ax in range(-k, k+1):
            for ay in range(-k, k+1):
                if 0 <= x+ax < W and 0 <= y+ay < H:
                    C.add((x+ax, y+ay))
    return C

def corridor_astar(t, s, g, C, w=None):
    """A* on C with exit bookkeeping. If w is None: exhaustive. Else early stop.
    Returns (U, Phi, expansions, stopped_early)."""
    k = tmin(t)
    h = lambda n: k*octile(n, g)
    gval = {s: 0.0}
    closed = set()
    OPEN = [(h(s), 0.0, s)]
    EXIT = {}   # node -> f
    U = math.inf
    exps = 0
    while OPEN:
        f, gn, n = heapq.heappop(OPEN)
        if n in closed or gn > gval[n] + 1e-12: continue
        if n not in C:
            continue  # never expand exits (they are tracked in EXIT dict)
        closed.add(n); exps += 1
        if n == g:
            U = gn
        if w is not None and U < math.inf:
            # early stop check: min f over OPEN (valid entries) and EXIT
            fmin = math.inf
            for ff, gg, nn in OPEN:
                if nn in closed: continue
                if gg > gval[nn] + 1e-12: continue
                fmin = min(fmin, ff)
            if EXIT: fmin = min(fmin, min(EXIT.values()))
            if fmin >= U / w:
                return U, fmin, exps, True
        for m, c in neighbours(t, n):
            ng = gn + c
            if ng < gval.get(m, math.inf) - 1e-12:
                gval[m] = ng
                fm = ng + h(m)
                if m in C:
                    heapq.heappush(OPEN, (fm, ng, m))
                else:
                    EXIT[m] = fm
    Phi = min(EXIT.values()) if EXIT else math.inf
    return U, Phi, exps, False

def main(seed=0, trials=400, H=40, W=40):
    rng = random.Random(seed)
    viol = 0; cons_viol = 0; early_viol = 0; mono_viol = 0; phi_nonmono = 0; n_tested = 0
    cert_rate = {1.0:0, 1.05:0, 1.2:0, 1.5:0}
    for tr in range(trials):
        pb = rng.choice([0.1, 0.2, 0.3, 0.4])
        t = make_grid(H, W, pb, rng, nterr=rng.choice([1,2,4,8]))
        s = (0,0); g = (W-1,H-1)
        t[0][0] = 1.0 if t[0][0]==math.inf else t[0][0]
        t[H-1][W-1] = 1.0 if t[H-1][W-1]==math.inf else t[H-1][W-1]
        if check_consistency(t, g) > 1e-9: cons_viol += 1
        dist = dijkstra(t, s)
        if g not in dist: continue
        OPT = dist[g]; n_tested += 1
        prevU = math.inf; prevPhi = None
        for k in [0,1,2,3,5,8,13,21,34]:
            C = bresenham_band(H, W, s, g, k)
            U, Phi, _, _ = corridor_astar(t, s, g, C)
            if OPT < min(U, Phi) - 1e-9: viol += 1; print("VIOLATION exhaustive", tr, k, OPT, U, Phi)
            if U > prevU + 1e-9: mono_viol += 1; print("U non-monotone", tr, k, prevU, U)
            if prevPhi is not None and Phi < prevPhi - 1e-9: phi_nonmono += 1
            prevU, prevPhi = U, Phi
            for w in cert_rate:
                if U < math.inf and U <= w*Phi + 1e-12:
                    cert_rate[w] += 1
                    if U > w*OPT + 1e-9: viol += 1; print("CERT VIOLATION", tr, k, w, OPT, U, Phi)
            for w in [1.05, 1.2, 1.5, 2.0]:
                U2, Phi2, _, st = corridor_astar(t, s, g, C, w=w)
                if OPT < min(U2, Phi2) - 1e-9: early_viol += 1; print("VIOLATION early", tr, k, w, OPT, U2, Phi2)
                if st and U2 > w*OPT + 1e-9: early_viol += 1; print("EARLY CERT VIOLATION", tr, k, w, OPT, U2, Phi2)
    print(f"instances tested={n_tested}  consistency violations={cons_viol}")
    print(f"Theorem1 violations={viol}  early-stop violations={early_viol}  U-monotonicity violations={mono_viol}")
    print(f"Phi decreased when widening (non-monotone) count={phi_nonmono}")
    print("certified (U<=w*Phi) counts over all (instance,k):", cert_rate)

if __name__ == "__main__":
    main(seed=int(sys.argv[1]) if len(sys.argv) > 1 else 0)
