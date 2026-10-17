import os
"""Enumerate every 3x3 neighbourhood (centre free) x parent direction.
Compute Carlson F2 neighbourhood successors (uniform cost, no corner cutting, strict tiebroken comparison)
by brute force over all simple paths from p to n inside the 3x3 block, and compare with
(a) pruned_dirs of check/jps_touch.py (restricted to legal moves), and
(b) the jump-stop rule (straight: stop iff forced successor exists; diagonal: never stops for forced reasons).
"""
import itertools, math, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jps_touch import Grid, jps  # jps not used; we re-implement pruned_dirs textually below
SQ2 = math.sqrt(2)
DIRS = [(dx,dy) for dx in (-1,0,1) for dy in (-1,0,1) if dx or dy]

def pruned_dirs(G, n, p):
    # verbatim copy of jps_touch.pruned_dirs
    x, y = n
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

def straight_stop(G, x, y, dx, dy):
    if dx:
        return (not G.free(x-dx, y+1) and G.free(x, y+1)) or (not G.free(x-dx, y-1) and G.free(x, y-1))
    return (not G.free(x+1, y-dy) and G.free(x+1, y)) or (not G.free(x-1, y-dy) and G.free(x-1, y))

def best_tb(G, src, dst, cells):
    """min tiebroken cost (c, |last|) over all simple paths src->dst using only cells in `cells` (3x3 block)."""
    best = (math.inf, math.inf)
    def rec(cur, visited, cost, last):
        nonlocal best
        if cur == dst:
            if (cost, last) < best: best = (cost, last)
            return
        for d in DIRS:
            nxt = (cur[0]+d[0], cur[1]+d[1])
            if nxt not in cells or nxt in visited: continue
            if not G.can_move(cur[0], cur[1], d[0], d[1]): continue
            w = SQ2 if d[0] and d[1] else 1.0
            rec(nxt, visited|{nxt}, cost+w, w)
    rec(src, {src}, 0.0, 0.0)
    return best

mism = 0; stop_mism = 0; total = 0
for bits in itertools.product([True, False], repeat=8):
    # 5x5 grid, centre (2,2) free; 3x3 neighbourhood around it
    cells = [[False]*5 for _ in range(5)]
    cells[2][2] = True
    k = 0
    for (ox, oy) in DIRS:
        cells[2+oy][2+ox] = bits[k]; k += 1
    G = Grid(cells)
    x = (2,2)
    nb = {(2+ox, 2+oy) for (ox,oy) in DIRS if G.free(2+ox,2+oy)} | {x}
    for d in DIRS:
        p = (2-d[0], 2-d[1])
        if not G.can_move(p[0], p[1], d[0], d[1]):  # parent move must be legal
            continue
        total += 1
        f2 = set()
        for m in DIRS:
            n = (2+m[0], 2+m[1])
            if not G.can_move(2, 2, m[0], m[1]): continue
            w = SQ2 if m[0] and m[1] else 1.0
            via_x = ((SQ2 if d[0] and d[1] else 1.0) + w, w)
            alt = best_tb(G, p, n, nb)
            if not (alt < via_x):  # pruned iff strictly lower alternative exists
                f2.add(m)
        impl = {m for m in pruned_dirs(G, x, p) if G.can_move(2,2,m[0],m[1])}
        if f2 != impl:
            mism += 1
            if mism <= 5:
                print("MISMATCH d=",d,"bits=",bits,"F2=",sorted(f2),"impl=",sorted(impl))
        # stop rule equivalence
        nat = {d} if not (d[0] and d[1]) else {d,(d[0],0),(0,d[1])}
        forced = f2 - nat
        if d[0] and d[1]:
            if forced: stop_mism += 1; print("DIAG has forced under F2:", bits, forced)
        else:
            if bool(forced) != straight_stop(G, 2, 2, d[0], d[1]):
                stop_mism += 1; print("STRAIGHT stop mismatch", d, bits, forced)
print(f"cases={total} pruning mismatches={mism} stop-rule mismatches={stop_mism}")
