import os
"""Exhaustive check of Lemma C (all clauses) on all small grids.
Instrumented canonical JPS identical in rules to check/jps_touch.py.
For every cell v != s with f*(v) < C* (or every reachable v when goal unreachable / exhaustive mode):
  C1: some scan visits v with cost == g*(v)
  C2: the FIRST such visit occurs during expansion of j with f(j) <= f*(v)   (draft: 'some' such scan; we test the first)
  C3: if that scan enqueues v (returns v as jump point) then v is popped with g*(v)
  C4: v is never popped (closed) with g > g*(v) BEFORE the optimal visit (stronger than draft; informational)
Also counts: pops with suboptimal g (known JPS pathology), to document that Lemma C must tolerate it.
"""
import heapq, math, itertools, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jps_touch import Grid, dijkstra, octile, SQ2

def jps_instr(G, s, goal, exhaustive=False):
    h = lambda n: octile(n, goal)
    gval = {s: 0.0}; parent = {s: None}; closed = {}
    OPEN = [(h(s), 0.0, s)]
    visit = {}   # cell -> (min g, f of expanding node at that visit, enqueued_by_that_scan?)
    pops = []    # (cell, g)
    cur = [None, None]  # current expanding node, its f
    def note(cell, g):
        if g < visit.get(cell, (math.inf,))[0] - 1e-12:
            visit[cell] = [g, cur[1], False]
    def jump(x, y, dx, dy, g):
        while True:
            if not G.can_move(x, y, dx, dy): return None
            x, y = x+dx, y+dy; g = g + (SQ2 if dx and dy else 1.0)
            note((x,y), g)
            if (x,y) == goal: return (x,y), g
            if dx and dy:
                if jump(x, y, dx, 0, g) is not None or jump(x, y, 0, dy, g) is not None:
                    return (x,y), g
            else:
                if dx:
                    if (not G.free(x-dx, y+1) and G.free(x, y+1)) or (not G.free(x-dx, y-1) and G.free(x, y-1)):
                        return (x,y), g
                else:
                    if (not G.free(x+1, y-dy) and G.free(x+1, y)) or (not G.free(x-1, y-dy) and G.free(x-1, y)):
                        return (x,y), g
    def pruned_dirs(n, p):
        x, y = n
        if p is None:
            return [(dx,dy) for dx in (-1,0,1) for dy in (-1,0,1) if (dx or dy)]
        dx = (x - p[0]); dy = (y - p[1]); dx = (dx>0)-(dx<0); dy = (dy>0)-(dy<0)
        dirs = []
        if dx and dy: dirs += [(dx,0), (0,dy), (dx,dy)]
        elif dx:
            dirs.append((dx,0))
            if not G.free(x-dx, y+1) and G.free(x, y+1): dirs.append((0,1)); dirs.append((dx,1))
            if not G.free(x-dx, y-1) and G.free(x, y-1): dirs.append((0,-1)); dirs.append((dx,-1))
        else:
            dirs.append((0,dy))
            if not G.free(x+1, y-dy) and G.free(x+1, y): dirs.append((1,0)); dirs.append((1,dy))
            if not G.free(x-1, y-dy) and G.free(x-1, y): dirs.append((-1,0)); dirs.append((-1,dy))
        return dirs
    res = math.inf
    while OPEN:
        f, g, n = heapq.heappop(OPEN)
        if n in closed or g > gval[n] + 1e-12: continue
        closed[n] = g; pops.append((n, g))
        if n == goal and not exhaustive: res = g; break
        cur[0], cur[1] = n, f
        for (dx,dy) in pruned_dirs(n, parent[n]):
            r = jump(n[0], n[1], dx, dy, g)
            if r is None: continue
            m, gm = r
            if abs(gm - visit[m][0]) < 1e-12: visit[m][2] = True
            if gm < gval.get(m, math.inf) - 1e-12:
                gval[m] = gm; parent[m] = n
                heapq.heappush(OPEN, (gm + h(m), gm, m))
    return res, visit, closed, pops

def check(G, s, goal, dist, exhaustive):
    res, visit, closed, pops = jps_instr(G, s, goal, exhaustive)
    Cstar = dist.get(goal, math.inf)
    if not exhaustive and abs(res - Cstar) > 1e-9: return "NONOPT"
    bound = math.inf if exhaustive else Cstar
    subopt_pops = sum(1 for (n,g) in pops if g > dist[n] + 1e-9)
    for v, gv in dist.items():
        if v == s: continue
        fv = gv + octile(v, goal)
        if fv < bound - 1e-9:
            if v not in visit or visit[v][0] > gv + 1e-9: return ("C1", v)
            if visit[v][1] > fv + 1e-9: return ("C2", v, visit[v][1], fv)
            if visit[v][2]:
                if v not in closed or abs(closed[v] - gv) > 1e-9: return ("C3", v)
    return ("OK", subopt_pops)

def main(H, W, mode):
    fails = {}; n = 0; sp = 0; grids_with_subopt = 0
    cells_n = H*W
    for bits in itertools.product([True, False], repeat=cells_n):
        cells = [list(bits[r*W:(r+1)*W]) for r in range(H)]
        G = Grid(cells)
        free = [(x,y) for y in range(H) for x in range(W) if cells[y][x]]
        for s in free:
            dist = dijkstra(G, s)
            for t in free:
                if t == s: continue
                exhaustive = (mode == 'exh')
                r = check(G, s, t, dist, exhaustive)
                n += 1
                if r == "NONOPT" or r[0] != "OK":
                    fails[r[0] if r != "NONOPT" else r] = fails.get(r[0] if r!="NONOPT" else r, 0) + 1
                    if sum(fails.values()) <= 3:
                        print("FAIL", r, "s=",s,"t=",t); [print(''.join('.' if c else '#' for c in row)) for row in cells]
                else:
                    if r[1] > 0: sp += r[1]; grids_with_subopt += 1
    print(f"{H}x{W} mode={mode}: instances={n} fails={fails} suboptimal-g pops={sp} in {grids_with_subopt} instances")

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3])
