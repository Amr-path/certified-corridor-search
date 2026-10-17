"""Computational probe of the 'touch property' of canonical JPS (uniform 8-connected grid,
no corner cutting, Harabor & Grastien 2014 pruning rules).

Property P: for every cell v with f*(v) = g*(v) + h(v) < C* (C* = optimal s-t cost), some
jump scan visits v with accumulated g equal to g*(v) BEFORE the goal is popped.

We also log the weaker property P': v is touched with SOME g (any value) before goal pop.
If P fails we print a counterexample.
"""
import heapq, math, random, sys
SQ2 = math.sqrt(2)

def octile(a, b):
    dx, dy = abs(a[0]-b[0]), abs(a[1]-b[1]); return max(dx,dy)+(SQ2-1)*min(dx,dy)

class Grid:
    def __init__(self, cells):
        self.c = cells; self.H = len(cells); self.W = len(cells[0])
    def free(self, x, y):
        return 0 <= x < self.W and 0 <= y < self.H and self.c[y][x]
    def can_move(self, x, y, dx, dy):
        """legal move from (x,y) by (dx,dy): destination free; diagonal needs both side cells free."""
        if not self.free(x+dx, y+dy): return False
        if dx and dy:
            return self.free(x+dx, y) and self.free(x, y+dy)
        return True

def dijkstra(G, s):
    dist = {s: 0.0}; pq = [(0.0, s)]
    while pq:
        d, (x,y) = heapq.heappop(pq)
        if d > dist[(x,y)]: continue
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                if dx==0 and dy==0: continue
                if G.can_move(x,y,dx,dy):
                    nd = d + (SQ2 if dx and dy else 1.0)
                    if nd < dist.get((x+dx,y+dy), math.inf):
                        dist[(x+dx,y+dy)] = nd; heapq.heappush(pq,(nd,(x+dx,y+dy)))
    return dist

def jps(G, s, goal, touch):
    """Canonical JPS. touch: dict cell -> min g with which a scan visited the cell (before goal pop)."""
    h = lambda n: octile(n, goal)
    gval = {s: 0.0}; parent = {s: None}; closed = set()
    OPEN = [(h(s), 0.0, s)]

    def note(cell, g):
        if g < touch.get(cell, math.inf): touch[cell] = g

    def jump(x, y, dx, dy, g):
        """jump from (x,y) in direction (dx,dy) with g value g at (x,y). returns (cell, g) or None."""
        while True:
            if not G.can_move(x, y, dx, dy): return None
            x, y = x+dx, y+dy; g = g + (SQ2 if dx and dy else 1.0)
            note((x,y), g)
            if (x,y) == goal: return (x,y), g
            if dx and dy:
                # forced neighbours for diagonal under no-corner-cutting: none (2014 rules); check orthogonal jumps
                if jump(x, y, dx, 0, g) is not None or jump(x, y, 0, dy, g) is not None:
                    return (x,y), g
            else:
                if dx:  # horizontal: forced perpendicular neighbour when the cell behind-perpendicular is blocked
                    if (not G.free(x-dx, y+1) and G.free(x, y+1)) or (not G.free(x-dx, y-1) and G.free(x, y-1)):
                        return (x,y), g
                else:   # vertical
                    if (not G.free(x+1, y-dy) and G.free(x+1, y)) or (not G.free(x-1, y-dy) and G.free(x-1, y)):
                        return (x,y), g

    def pruned_dirs(n, p):
        x, y = n
        if p is None:
            return [(dx,dy) for dx in (-1,0,1) for dy in (-1,0,1) if (dx or dy)]
        dx = (x - p[0]); dy = (y - p[1])
        dx = (dx>0)-(dx<0); dy = (dy>0)-(dy<0)
        dirs = []
        if dx and dy:
            dirs += [(dx,0), (0,dy), (dx,dy)]
        elif dx:
            dirs.append((dx,0))
            if not G.free(x-dx, y+1) and G.free(x, y+1): dirs.append((0,1)); dirs.append((dx,1))
            if not G.free(x-dx, y-1) and G.free(x, y-1): dirs.append((0,-1)); dirs.append((dx,-1))
        else:
            dirs.append((0,dy))
            if not G.free(x+1, y-dy) and G.free(x+1, y): dirs.append((1,0)); dirs.append((1,dy))
            if not G.free(x-1, y-dy) and G.free(x-1, y): dirs.append((-1,0)); dirs.append((-1,dy))
        return dirs

    while OPEN:
        f, g, n = heapq.heappop(OPEN)
        if n in closed or g > gval[n] + 1e-12: continue
        closed.add(n)
        if n == goal: return g
        for (dx,dy) in pruned_dirs(n, parent[n]):
            r = jump(n[0], n[1], dx, dy, g)
            if r is None: continue
            m, gm = r
            if gm < gval.get(m, math.inf) - 1e-12:
                gval[m] = gm; parent[m] = n
                heapq.heappush(OPEN, (gm + h(m), gm, m))
    return math.inf

def main(seed, trials, H, W):
    rng = random.Random(seed)
    P_fail = P2_fail = 0; tested = 0; opt_mismatch = 0
    for tr in range(trials):
        pb = rng.choice([0.1,0.2,0.3,0.35,0.4])
        cells = [[rng.random() > pb for _ in range(W)] for _ in range(H)]
        s = (rng.randrange(W), rng.randrange(H)); g = (rng.randrange(W), rng.randrange(H))
        cells[s[1]][s[0]] = True; cells[g[1]][g[0]] = True
        if s == g: continue
        G = Grid(cells)
        dist = dijkstra(G, s)
        if g not in dist: continue
        Cstar = dist[g]; tested += 1
        touch = {s: 0.0}
        res = jps(G, s, g, touch)
        if abs(res - Cstar) > 1e-9: opt_mismatch += 1; print("JPS NOT OPTIMAL?!", tr, res, Cstar)
        bad = None
        for v, gv in dist.items():
            if gv + octile(v, g) < Cstar - 1e-9:
                if v not in touch:
                    P2_fail += 1; bad = (v, gv, None); break
                if touch[v] > gv + 1e-9:
                    P_fail += 1; bad = (v, gv, touch[v]); break
        if bad:
            print(f"trial {tr}: property fails at cell {bad[0]} g*={bad[1]:.3f} touched_g={bad[2]} C*={Cstar:.3f} s={s} g={g} pb={pb}")
            if P_fail + P2_fail <= 2:
                for y in range(H):
                    row = ''
                    for x in range(W):
                        ch = '.' if cells[y][x] else '#'
                        if (x,y) == s: ch = 'S'
                        if (x,y) == g: ch = 'G'
                        if (x,y) == bad[0]: ch = 'X'
                        row += ch
                    print(row)
    print(f"seed={seed} tested={tested} optimality mismatches={opt_mismatch} P-fail(touched with worse g)={P_fail} P'-fail(never touched)={P2_fail}")

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
