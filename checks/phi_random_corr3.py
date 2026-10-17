import os, sys, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from phi_goalclose import *
import heapq

def astar_gc(t,s,g,C,tb,rng):
    k=tmin(t); h=lambda n: Q(k*octile(n,g).a, k*octile(n,g).b)
    gval={s:Q(0)}; closed=set(); cnt=[0]
    def key(f,gg):
        cnt[0]+=1
        if tb=='lowg': return (f.val(), gg.val(), cnt[0])
        if tb=='highg': return (f.val(), -gg.val(), cnt[0])
        if tb=='lifo': return (f.val(), -cnt[0])
        if tb=='fifo': return (f.val(), cnt[0])
        if tb=='rand': return (f.val(), rng.random())
        if tb=='id': return (f.val(), 0)  # placeholder, replaced below
    PERM=getattr(rng,'perm',None)
    def key2(f,gg,n):
        if tb=='id': return (f.val(), n[1]*100+n[0])
        if tb=='rid': return (f.val(), PERM[n])
        if tb=='idhg': return (f.val(), -gg.val(), n[1]*100+n[0])
        if tb=='ridhg': return (f.val(), -gg.val(), PERM[n])
        return key(f,gg)
    OPEN=[(key2(h(s),Q(0),s),s,Q(0))]; EXIT={}; parent={s:None}
    while OPEN:
        _,n,gn=heapq.heappop(OPEN)
        if n in closed or not gn==gval[n]: continue
        closed.add(n)
        if n==g:
            return gn,(min(EXIT.values()) if EXIT else INF),EXIT,closed,gval,parent
        for m,c in nbrs(t,n):
            ng=gn+c
            if m not in gval or ng<gval[m]:
                gval[m]=ng; parent[m]=n
                if m in C: heapq.heappush(OPEN,(key2(ng+h(m),ng,m),m,ng))
                else: EXIT[m]=ng+h(m)
    return INF,(min(EXIT.values()) if EXIT else INF),EXIT,closed,gval,parent

def main(seed,trials,H,W,tb,weighted):
    rng=random.Random(seed); found=0; tested=0
    for tr in range(trials):
        pb=rng.choice([0.0,0.1,0.2,0.3])
        t=[[ (None if rng.random()<pb else (rng.choice([1,2]) if weighted else 1)) for _ in range(W)] for _ in range(H)]
        s=(rng.randrange(W),rng.randrange(H)); g=(rng.randrange(W),rng.randrange(H))
        if s==g: continue
        t[s[1]][s[0]]=1; t[g[1]][g[0]]=1
        trav=[(x,y) for y in range(H) for x in range(W) if t[y][x] is not None]
        allc=[(x,y) for y in range(H) for x in range(W)]; pr=allc[:]; rng.shuffle(pr); rng.perm={n:i for i,n in enumerate(pr)}
        p1=rng.uniform(0.3,0.8); p2=rng.uniform(0.0,0.7)
        C={n for n in trav if rng.random()<p1}|{s,g}
        Cp=C|{n for n in trav if rng.random()<p2}
        U,Phi,EXIT,cl,gv,pa=astar_gc(t,s,g,C,tb,rng)
        Up,Php,EXITp,clp,gvp,pap=astar_gc(t,s,g,Cp,tb,rng)
        tested+=1
        if Up>U: print("U NONMONO?!",tr)
        if min(U,Phi)>min(Up,Php): print("L NONMONO (Lemma2ii VIOLATED)",tr,U,Phi,Up,Php)
        if Php<Phi:
            found+=1
            if found<=3:
                print(f"--- seed={seed} trial={tr} tb={tb}: Phi'<Phi : U={U} U'={Up} Phi={Phi} Phi'={Php}")
                print(show(t,s,g,C,Cp,EXIT,EXITp))
                bp=min(EXITp,key=lambda b:EXITp[b].val()); ch=[];x=bp
                while x is not None: ch.append(x); x=pap[x]
                print("b'=",bp,"rho=",ch[::-1]); print("EXIT(C)=",{b:round(v.val(),3) for b,v in EXIT.items()})
                print("OPEN-ish(C) unclosed generated in C:",{n:round((gv[n]+Q(tmin(t)*octile(n,g).a,tmin(t)*octile(n,g).b)).val(),3) for n in gv if n in C and n not in cl})
    print(f"seed={seed} tb={tb} weighted={weighted} {H}x{W} tested={tested}: Phi decreases={found}")

if __name__=="__main__":
    main(int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5],sys.argv[6]=='1')
