// ccs.cpp — Certified Corridor Search experimental framework.
// Single-file C++17. Implements, on MovingAI grid maps:
//   astar       : A* (octile, 8-connected, no corner cutting), whole map
//   jps         : canonical JPS, cell-by-cell scans, exactly the rules of Sec. 5.1
//   castar      : corridor A* on G_C^+ with EXIT bookkeeping (Sec. 4)
//   cjps        : corridor JPS on M_C with tau bookkeeping (Sec. 5.5)
//   wastar      : Weighted A* without re-expansion (whole map)
//   fls         : Focal List Search (Carlson, Harabor, Stuckey 2024), whole map
//   ccs         : Certified Corridor Search, geometric widening, inner = castar or cjps
// Per-instance CSV output. Uniform-cost and weighted (terrain) modes.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>
#include <limits>
#include <cstdint>

static const double SQ2 = std::sqrt(2.0);
static const double INF = std::numeric_limits<double>::infinity();

struct Map {
    int W = 0, H = 0;
    std::vector<uint8_t> pass;    // 1 traversable
    std::vector<double> cost;     // terrain cost (weighted mode); 1 in uniform mode
    double tmin = 1.0;
    bool weighted = false;
    inline int id(int x, int y) const { return y * W + x; }
    inline bool in(int x, int y) const { return x >= 0 && y >= 0 && x < W && y < H; }
    inline bool free_(int x, int y) const { return in(x, y) && pass[id(x, y)]; }
};

// Terrain conventions.
// uniform: '.' and 'G' traversable (warthog convention), everything else blocked.
// weighted (Carlson et al. 2024 assignment): '.'/'G' 1, 'T' 1.5, 'S' 2, 'W' 4, '@'/'O' 10;
//   letters 'A'..'P' (MovingAI "weighted" set) -> 1..16. Nothing is blocked in weighted mode
//   except characters outside this list.
static bool g_letters = false;
static bool load_map(const std::string& fn, Map& m, bool weighted) {
    std::ifstream f(fn);
    if (!f) return false;
    std::string tok;
    f >> tok;               // type
    f >> tok;               // octile
    f >> tok >> m.H;        // height
    f >> tok >> m.W;        // width
    f >> tok;               // map
    m.pass.assign(m.W * m.H, 0);
    m.cost.assign(m.W * m.H, 1.0);
    m.weighted = weighted;
    m.tmin = INF;
    for (int y = 0; y < m.H; y++) {
        std::string line;
        do { std::getline(f, line); if (!line.empty() && line.back() == '\r') line.pop_back(); } while (line.empty() && f);
        for (int x = 0; x < m.W && x < (int)line.size(); x++) {
            char c = line[x];
            int i = m.id(x, y);
            if (!weighted) {
                m.pass[i] = (c == '.' || c == 'G');
            } else {
                double t = -1;
                if (g_letters) { if (c >= 'A' && c <= 'P') t = 1 + (c - 'A'); }
                else if (c == '.' || c == 'G') t = 1;
                else if (c == 'T') t = 1.5;
                else if (c == 'S') t = 2;
                else if (c == 'W') t = 4;
                else if (c == '@' || c == 'O') t = 10;
                if (t > 0) { m.pass[i] = 1; m.cost[i] = t; m.tmin = std::min(m.tmin, t); }
            }
        }
    }
    if (!weighted) m.tmin = 1.0;
    if (m.tmin == INF) m.tmin = 1.0;
    return true;
}

struct Inst { int sx, sy, gx, gy; double opt; };
static bool load_scen(const std::string& fn, std::vector<Inst>& v, std::string& mapname) {
    std::ifstream f(fn);
    if (!f) return false;
    std::string line;
    std::getline(f, line); // version
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        int b, w, h; Inst in; std::string mn;
        ss >> b >> mn >> w >> h >> in.sx >> in.sy >> in.gx >> in.gy >> in.opt;
        mapname = mn; v.push_back(in);
    }
    return true;
}

static inline double octile(int dx, int dy) {
    dx = std::abs(dx); dy = std::abs(dy);
    return (dx > dy) ? dx + (SQ2 - 1) * dy : dy + (SQ2 - 1) * dx;
}

// Legal move test (no corner cutting) on an arbitrary passability bitmap.
struct Grid {
    const Map* m;
    const std::vector<uint8_t>* pass; // may be a corridor-masked bitmap
    inline bool fr(int x, int y) const { return m->in(x, y) && (*pass)[m->id(x, y)]; }
    inline bool legal(int x, int y, int dx, int dy) const {
        if (!fr(x + dx, y + dy)) return false;
        if (dx && dy) return fr(x + dx, y) && fr(x, y + dy);
        return true;
    }
    inline double d(int x, int y, int dx, int dy) const {
        if (!m->weighted) return (dx && dy) ? SQ2 : 1.0;
        int a = m->id(x, y), b = m->id(x + dx, y + dy);
        if (dx && dy) {
            int c = m->id(x + dx, y), e = m->id(x, y + dy);
            return SQ2 * (m->cost[a] + m->cost[b] + m->cost[c] + m->cost[e]) / 4.0;
        }
        return (m->cost[a] + m->cost[b]) / 2.0;
    }
};


// ---------------- persistent workspace (epoch-stamped) ----------------
struct WS {
    std::vector<uint32_t> stamp; std::vector<double> gv, tau; std::vector<uint8_t> closed; std::vector<int> par; std::vector<double> gexp;
    uint32_t epoch = 0; int n = 0;
    void init(int N) { if (N != n) { n = N; stamp.assign(N, 0); gv.assign(N, INF); tau.assign(N, INF); closed.assign(N, 0); par.assign(N, -1); gexp.assign(N, INF); epoch = 0; } }
    void reset() { epoch++; if (epoch == 0) { std::fill(stamp.begin(), stamp.end(), 0); epoch = 1; } }
    inline void touch(int i) { if (stamp[i] != epoch) { stamp[i] = epoch; gv[i] = INF; tau[i] = INF; closed[i] = 0; par[i] = -1; gexp[i] = INF; } }
    inline double G(int i) { touch(i); return gv[i]; }
};
static WS ws;

static const int DX[8] = {1, -1, 0, 0, 1, 1, -1, -1};
static const int DY[8] = {0, 0, 1, -1, 1, -1, 1, -1};

struct Stats {
    long expanded = 0, generated = 0, scanned = 0; // scanned: cells visited by JPS scans
    double U = INF, Phi = INF, L = INF;
    int rounds = 0;
};

// ---------------- Corridor: Bresenham band of half-width k, as per-row intervals ----------------
struct Band {
    int H; std::vector<int> lo, hi; // per row; lo>hi means empty
    bool full = false;
    inline bool in(int x, int y) const { return full || (y >= 0 && y < H && x >= lo[y] && x <= hi[y]); }
};
static Band make_band(const Map& m, int sx, int sy, int gx, int gy, int k) {
    Band b; b.H = m.H; b.lo.assign(m.H, 1 << 30); b.hi.assign(m.H, -1);
    if (k >= std::max(m.W, m.H)) { b.full = true; return b; }
    int x0 = sx, y0 = sy, x1 = gx, y1 = gy;
    int dx = std::abs(x1 - x0), dy = std::abs(y1 - y0);
    int stx = x0 < x1 ? 1 : -1, sty = y0 < y1 ? 1 : -1;
    int err = dx - dy, x = x0, y = y0;
    auto mark = [&](int cx, int cy) {
        for (int yy = std::max(0, cy - k); yy <= std::min(m.H - 1, cy + k); yy++) {
            b.lo[yy] = std::min(b.lo[yy], std::max(0, cx - k));
            b.hi[yy] = std::max(b.hi[yy], std::min(m.W - 1, cx + k));
        }
    };
    while (true) {
        mark(x, y);
        if (x == x1 && y == y1) break;
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; x += stx; }
        if (e2 < dx) { err += dx; y += sty; }
    }
    return b;
}
// corridor-masked passability bitmap (cells outside band treated as blocked)
// persistent mask buffer: rows of the previous band are cleared before the next band is written
static std::vector<uint8_t> g_mask; static Band g_mask_band; static bool g_mask_used = false;
static const std::vector<uint8_t>& mask_band(const Map& m, const Band& b) {
    if (b.full) return m.pass;
    if ((int)g_mask.size() != m.W * m.H) { g_mask.assign(m.W * m.H, 0); g_mask_used = false; }
    if (g_mask_used) for (int y = 0; y < m.H; y++) for (int x = std::max(0, g_mask_band.lo[y]); x <= g_mask_band.hi[y]; x++) g_mask[m.id(x, y)] = 0;
    for (int y = 0; y < m.H; y++)
        for (int x = std::max(0, b.lo[y]); x <= b.hi[y]; x++)
            g_mask[m.id(x, y)] = m.pass[m.id(x, y)];
    g_mask_band = b; g_mask_used = true;
    return g_mask;
}

// ---------------- generic bucketed-free binary heap over node ids with (f, g) keys ----------------
struct PQ {
    struct E { double f, g; int id; };
    std::vector<E> h;
    // tie-break: smaller f; then larger g (standard) — configurable
    static bool less(const E& a, const E& b) { if (a.f != b.f) return a.f > b.f; return a.g < b.g; } // for max-heap semantics via std::push_heap we invert
    void push(double f, double g, int id) { h.push_back({f, g, id}); std::push_heap(h.begin(), h.end(), less); }
    E pop() { std::pop_heap(h.begin(), h.end(), less); E e = h.back(); h.pop_back(); return e; }
    bool empty() const { return h.empty(); }
    const E& top() const { return h.front(); }
};

// ---------------- A* family ----------------
// corridor: if band != nullptr, run on G_C^+ : nodes inside band expanded, nodes outside generated into EXIT.
// mode 0: plain A* (band may be null). mode 1: WA* without re-expansion (w>1). mode 2: FLS.
struct Search {
    const Map& m; Grid g;
    struct GV { double operator[](int i) const { ws.touch(i); return ws.gv[i]; } } gvr;
    Search(const Map& mm, const std::vector<uint8_t>* pass) : m(mm) { g.m = &mm; g.pass = pass; ws.init(m.W * m.H); ws.reset(); }
    inline double getg(int i) { ws.touch(i); return ws.gv[i]; }
    inline void setg(int i, double v) { ws.touch(i); ws.gv[i] = v; }
    inline bool isclosed(int i) { ws.touch(i); return ws.closed[i]; }
    inline void close(int i) { ws.touch(i); ws.closed[i] = 1; }
    inline double getgexp(int i) { ws.touch(i); return ws.gexp[i]; }
    inline void setgexp(int i, double v) { ws.touch(i); ws.gexp[i] = v; }
};

static Stats run_astar(const Map& m, int sx, int sy, int gx, int gy, const Band* band, double w, int mode) {
    // mode 0: A* (optionally corridor on G_C^+). mode 1: WA* without re-expansion. mode 2: FLS (Carlson et al. 2024, Alg. 1).
    Stats st;
    Search S(m, &m.pass);
        auto h = [&](int x, int y) { return m.tmin * octile(x - gx, y - gy); };
    int s = m.id(sx, sy), t = m.id(gx, gy);
    PQ open, focal;
    S.setg(s, 0); open.push(h(sx, sy), 0, s); st.generated++;
    if (mode == 2) focal.push(w * h(sx, sy), 0, s);
    double Phi = INF;
    while (true) {
        // drop stale entries
        while (!open.empty() && (open.top().g > S.getg(open.top().id) + 1e-12 || open.top().g >= S.getgexp(open.top().id) - 1e-12)) open.pop();
        if (mode == 2) while (!focal.empty() && (focal.top().g > S.getg(focal.top().id) + 1e-12 || S.isclosed(focal.top().id))) focal.pop();
        if (open.empty()) break;
        int n; double gn;
        if (mode == 2 && !focal.empty() && focal.top().f <= w * open.top().f + 1e-12) { auto e = focal.pop(); n = e.id; gn = e.g; }
        else { auto e = open.pop(); n = e.id; gn = e.g; }
        if (mode != 2 && S.isclosed(n)) continue; // A*/WA*: never re-expand (A*: consistent h; WA*: by design)
        S.close(n); S.setgexp(n, gn); st.expanded++;
        if (n == t) { st.U = gn; break; }
        int x = n % m.W, y = n / m.W;
        for (int k = 0; k < 8; k++) {
            int dx = DX[k], dy = DY[k];
            if (!S.g.legal(x, y, dx, dy)) continue;
            int nx = x + dx, ny = y + dy, nn = m.id(nx, ny);
            double ng = gn + S.g.d(x, y, dx, dy);
            if (ng < S.getg(nn) - 1e-12) {
                if (mode == 1 && S.isclosed(nn)) continue;   // WA* without re-expansion ignores improvements to closed nodes
                S.setg(nn, ng); st.generated++;
                double hf = h(nx, ny);
                if (band && !band->in(nx, ny)) { Phi = std::min(Phi, ng + hf); continue; } // exit cell: generated, never expanded
                if (mode == 1) open.push(ng + w * hf, ng, nn);
                else { open.push(ng + hf, ng, nn); if (mode == 2 && !S.isclosed(nn)) focal.push(ng + w * hf, ng, nn); }
            }
        }
    }
    st.Phi = Phi; st.L = std::min(st.U, Phi);
    return st;
}

// ---------------- canonical JPS (Sec. 5.1 rules), with optional corridor map + tau ----------------
struct JPS {
    const Map& m; Grid g; int gx, gy;
    bool usetau; long scanned = 0;
    JPS(const Map& mm, const std::vector<uint8_t>* pass, int gx_, int gy_, bool usetau_)
        : m(mm), gx(gx_), gy(gy_), usetau(usetau_) { g.m = &mm; g.pass = pass; ws.init(m.W*m.H); ws.reset(); }
    inline void visit(int x, int y, double gg) { scanned++; if (usetau) { int i = m.id(x, y); ws.touch(i); if (gg < ws.tau[i]) ws.tau[i] = gg; } }
    inline bool blocked(int x, int y) const { return !g.fr(x, y); }
    // forced-direction test for straight arrival d at cell (x,y): some perpendicular e with x-d+e blocked and x+e free
    inline bool has_forced(int x, int y, int dx, int dy) const {
        if (dx) return (blocked(x - dx, y + 1) && g.fr(x, y + 1)) || (blocked(x - dx, y - 1) && g.fr(x, y - 1));
        return (blocked(x + 1, y - dy) && g.fr(x + 1, y)) || (blocked(x - 1, y - dy) && g.fr(x - 1, y));
    }
    // Jump: returns true and sets (jx,jy,jg) if a jump point is returned.
    bool jump(int x, int y, int dx, int dy, double gg, int& jx, int& jy, double& jg) {
        while (true) {
            if (!g.legal(x, y, dx, dy)) return false;
            x += dx; y += dy; gg += (dx && dy) ? SQ2 : 1.0;
            visit(x, y, gg);
            if (x == gx && y == gy) { jx = x; jy = y; jg = gg; return true; }
            if (dx && dy) {
                int a, b; double c;
                if (jump(x, y, dx, 0, gg, a, b, c) || jump(x, y, 0, dy, gg, a, b, c)) { jx = x; jy = y; jg = gg; return true; }
            } else {
                if (has_forced(x, y, dx, dy)) { jx = x; jy = y; jg = gg; return true; }
            }
        }
    }
    // successor directions for arrival direction (dx,dy); (0,0) = start
    int dirs(int x, int y, int dx, int dy, int* ox, int* oy) {
        int n = 0;
        if (!dx && !dy) { for (int k = 0; k < 8; k++) { ox[n] = DX[k]; oy[n] = DY[k]; n++; } return n; }
        if (dx && dy) { ox[n]=dx; oy[n]=dy; n++; ox[n]=dx; oy[n]=0; n++; ox[n]=0; oy[n]=dy; n++; return n; }
        ox[n]=dx; oy[n]=dy; n++;
        if (dx) {
            if (blocked(x - dx, y + 1) && g.fr(x, y + 1)) { ox[n]=0; oy[n]=1; n++; ox[n]=dx; oy[n]=1; n++; }
            if (blocked(x - dx, y - 1) && g.fr(x, y - 1)) { ox[n]=0; oy[n]=-1; n++; ox[n]=dx; oy[n]=-1; n++; }
        } else {
            if (blocked(x + 1, y - dy) && g.fr(x + 1, y)) { ox[n]=1; oy[n]=0; n++; ox[n]=1; oy[n]=dy; n++; }
            if (blocked(x - 1, y - dy) && g.fr(x - 1, y)) { ox[n]=-1; oy[n]=0; n++; ox[n]=-1; oy[n]=dy; n++; }
        }
        return n;
    }
    double search(int sx, int sy, Stats& st) {
        int s = m.id(sx, sy), t = m.id(gx, gy);
        PQ open; ws.touch(s); ws.gv[s] = 0; open.push(octile(sx - gx, sy - gy), 0, s); st.generated++;
        if (usetau) ws.tau[s] = 0;
        while (!open.empty()) {
            auto e = open.pop(); int n = e.id; double gn = e.g;
            ws.touch(n); if (gn > ws.gv[n] || ws.closed[n]) continue;
            ws.closed[n] = 1; st.expanded++;
            if (n == t) return gn;
            int x = n % m.W, y = n / m.W;
            int dx = 0, dy = 0;
            if (ws.par[n] >= 0) { int px = ws.par[n] % m.W, py = ws.par[n] / m.W; dx = (x > px) - (x < px); dy = (y > py) - (y < py); }
            int ox[8], oy[8]; int nd = dirs(x, y, dx, dy, ox, oy);
            for (int i = 0; i < nd; i++) {
                int jx, jy; double jg;
                if (!jump(x, y, ox[i], oy[i], gn, jx, jy, jg)) continue;
                int jn = m.id(jx, jy);
                ws.touch(jn); if (jg < ws.gv[jn] - 1e-12) { ws.gv[jn] = jg; ws.par[jn] = n; st.generated++; open.push(jg + octile(jx - gx, jy - gy), jg, jn); }
            }
        }
        return INF;
    }
};

// corridor JPS: on M_C (band mask), tau over boundary cells, beta via original-map moves not legal in M_C
static Stats run_cjps(const Map& m, int sx, int sy, int gx, int gy, const Band* band, double* PhiA_out) {
    Stats st;
    const std::vector<uint8_t>& mask = band ? mask_band(m, *band) : m.pass;
    JPS J(m, &mask, gx, gy, true);
    st.U = J.search(sx, sy, st);
    st.scanned = J.scanned;
    // Phi_J = min over cells c with tau<inf of tau(c)+beta(c)
    double Phi = INF;
    if (band && !band->full) {
        Grid gfull{&m, &m.pass}, gmask{&m, &mask};
        for (int y = 0; y < m.H; y++) { if (band->lo[y] > band->hi[y]) continue; for (int x = band->lo[y]; x <= band->hi[y]; x++) {
            int i = m.id(x, y); if (ws.stamp[i] != ws.epoch || ws.tau[i] == INF) continue;
            for (int k = 0; k < 8; k++) {
                int dx = DX[k], dy = DY[k];
                if (gfull.legal(x, y, dx, dy) && !gmask.legal(x, y, dx, dy))
                    Phi = std::min(Phi, ws.tau[i] + gfull.d(x, y, dx, dy) + octile(x + dx - gx, y + dy - gy));
            }
        }}
    }
    st.Phi = Phi; st.L = std::min(st.U, Phi);
    return st;
}
// corridor A* with the SAME corridor-map semantics (for the tightness check): A* on mask with beta terms
static Stats run_castar_mc(const Map& m, int sx, int sy, int gx, int gy, const Band* band) {
    Stats st;
    const std::vector<uint8_t>& mask = band ? mask_band(m, *band) : m.pass;
    Search S(m, &mask);
    auto h = [&](int x, int y) { return m.tmin * octile(x - gx, y - gy); };
    int s = m.id(sx, sy), t = m.id(gx, gy);
    PQ open; S.setg(s, 0); open.push(h(sx, sy), 0, s);
    while (!open.empty()) {
        auto e = open.pop(); int n = e.id; double gn = e.g;
        if (gn > S.getg(n) || S.isclosed(n)) continue;
        S.close(n); st.expanded++;
        if (n == t) { st.U = gn; break; }
        int x = n % m.W, y = n / m.W;
        for (int k = 0; k < 8; k++) {
            int dx = DX[k], dy = DY[k];
            if (!S.g.legal(x, y, dx, dy)) continue;
            int nn = m.id(x + dx, y + dy); double ng = gn + S.g.d(x, y, dx, dy);
            if (ng < S.getg(nn) - 1e-12) { S.setg(nn, ng); st.generated++; open.push(ng + h(x + dx, y + dy), ng, nn); }
        }
    }
    double Phi = INF;
    if (band && !band->full) {
        Grid gfull{&m, &m.pass}, gmask{&m, &mask};
        for (int y = 0; y < m.H; y++) { if (band->lo[y] > band->hi[y]) continue; for (int x = band->lo[y]; x <= band->hi[y]; x++) {
            int i = m.id(x, y); if (ws.stamp[i] != ws.epoch || ws.gv[i] == INF) continue;
            for (int k = 0; k < 8; k++) {
                int dx = DX[k], dy = DY[k];
                if (gfull.legal(x, y, dx, dy) && !gmask.legal(x, y, dx, dy))
                    Phi = std::min(Phi, ws.gv[i] + gfull.d(x, y, dx, dy) + m.tmin * octile(x + dx - gx, y + dy - gy));
            }
        }}
    }
    st.Phi = Phi; st.L = std::min(st.U, Phi);
    return st;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: ccs <alg> <scen> <mapdir> [--weighted] [--k K] [--w W] [--every N] [--inner castar|cjps] [--k0 K0]\n"
                        "alg: astar jps castar cjps castar_mc wastar fls ccs pair\n");
        return 1;
    }
    std::string alg = argv[1], scen = argv[2], mapdir = argv[3];
    bool weighted = false; int k = 8, every = 1, k0 = 2; double w = 1.2, k0rel = 0; std::string inner = "castar";
    for (int i = 4; i < argc; i++) {
        std::string a = argv[i];
        if (a == "--weighted") weighted = true;
        else if (a == "--k") k = atoi(argv[++i]);
        else if (a == "--w") w = atof(argv[++i]);
        else if (a == "--every") every = atoi(argv[++i]);
        else if (a == "--inner") inner = argv[++i];
        else if (a == "--k0") k0 = atoi(argv[++i]);
        else if (a == "--k0rel") k0rel = atof(argv[++i]);
        else if (a == "--letters") { g_letters = true; weighted = true; }
    }
    std::vector<Inst> insts; std::string mapname;
    if (!load_scen(scen, insts, mapname)) { fprintf(stderr, "bad scen\n"); return 1; }
    if (weighted && (alg == "jps" || alg == "cjps" || inner == "cjps")) { fprintf(stderr, "JPS variants are uniform-cost only\n"); return 1; }
    Map m; if (!load_map(mapdir + "/" + mapname, m, weighted)) { fprintf(stderr, "bad map %s\n", mapname.c_str()); return 1; }
    printf("id,alg,map,k,w,expanded,generated,scanned,micros,U,Phi,L,opt,rounds,extra\n");
    for (size_t i = 0; i < insts.size(); i += every) {
        const Inst& in = insts[i];
        if (!m.pass[m.id(in.sx, in.sy)] || !m.pass[m.id(in.gx, in.gy)]) continue;
        Stats st; double extra = 0;
        auto t0 = std::chrono::steady_clock::now();
        if (alg == "astar") st = run_astar(m, in.sx, in.sy, in.gx, in.gy, nullptr, 1.0, 0);
        else if (alg == "wastar") st = run_astar(m, in.sx, in.sy, in.gx, in.gy, nullptr, w, 1);
        else if (alg == "fls") st = run_astar(m, in.sx, in.sy, in.gx, in.gy, nullptr, w, 2);
        else if (alg == "castar") { Band b = make_band(m, in.sx, in.sy, in.gx, in.gy, k); st = run_astar(m, in.sx, in.sy, in.gx, in.gy, &b, 1.0, 0); }
        else if (alg == "castar_mc") { Band b = make_band(m, in.sx, in.sy, in.gx, in.gy, k); st = run_castar_mc(m, in.sx, in.sy, in.gx, in.gy, &b); }
        else if (alg == "jps") { Stats s2; JPS J(m, &m.pass, in.gx, in.gy, false); s2.U = J.search(in.sx, in.sy, s2); s2.scanned = J.scanned; s2.L = s2.U; st = s2; }
        else if (alg == "cjps") { Band b = make_band(m, in.sx, in.sy, in.gx, in.gy, k); st = run_cjps(m, in.sx, in.sy, in.gx, in.gy, &b, nullptr); }
        else if (alg == "ccs") {
            int kk = k0; if (k0rel > 0) kk = std::max(k0, (int)std::ceil(k0rel * octile(in.gx - in.sx, in.gy - in.sy))); int rounds = 0; Stats acc; double bestU = INF, bestL = 0;
            while (true) {
                rounds++;
                Band b = make_band(m, in.sx, in.sy, in.gx, in.gy, kk);
                Stats r = (inner == "cjps") ? run_cjps(m, in.sx, in.sy, in.gx, in.gy, &b, nullptr)
                                            : run_astar(m, in.sx, in.sy, in.gx, in.gy, &b, 1.0, 0);
                acc.expanded += r.expanded; acc.generated += r.generated; acc.scanned += r.scanned;
                if (r.U < INF) { bestU = std::min(bestU, r.U); }
                bestL = std::max(bestL, r.L);
                if (r.U < INF && r.U <= w * r.Phi) { acc.U = r.U; acc.Phi = r.Phi; acc.L = r.L; break; }
                if (b.full) { acc.U = r.U; acc.Phi = r.Phi; acc.L = r.L; break; }
                kk *= 2;
            }
            acc.rounds = rounds; st = acc; extra = kk;
        }
        else if (alg == "pair") { // paired, back-to-back timing of corridor A* on M_C and JPS_C on the same instance and band
            Band b = make_band(m, in.sx, in.sy, in.gx, in.gy, k);
            auto a0 = std::chrono::steady_clock::now(); Stats sa = run_castar_mc(m, in.sx, in.sy, in.gx, in.gy, &b);
            auto a1 = std::chrono::steady_clock::now(); Stats sj = run_cjps(m, in.sx, in.sy, in.gx, in.gy, &b, nullptr);
            auto a2 = std::chrono::steady_clock::now();
            printf("%zu,%s,%s,%d,%g,%ld,%ld,%ld,%.1f,%.6f,%.6f,%.6f,%.6f,%d,%g\n", i, "pmc", mapname.c_str(), k, w, sa.expanded, sa.generated, sa.scanned,
                   std::chrono::duration<double, std::micro>(a1 - a0).count(), sa.U, sa.Phi, sa.L, in.opt, 0, 0.0);
            printf("%zu,%s,%s,%d,%g,%ld,%ld,%ld,%.1f,%.6f,%.6f,%.6f,%.6f,%d,%g\n", i, "pjps", mapname.c_str(), k, w, sj.expanded, sj.generated, sj.scanned,
                   std::chrono::duration<double, std::micro>(a2 - a1).count(), sj.U, sj.Phi, sj.L, in.opt, 0, 0.0);
            continue;
        }
        auto t1 = std::chrono::steady_clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        printf("%zu,%s,%s,%d,%g,%ld,%ld,%ld,%.1f,%.6f,%.6f,%.6f,%.6f,%d,%g\n", i, alg.c_str(), mapname.c_str(), k, w,
               st.expanded, st.generated, st.scanned, us, st.U, st.Phi, st.L, in.opt, st.rounds, extra);
    }
    return 0;
}
