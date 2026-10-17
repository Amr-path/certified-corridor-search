import math, heapq
W,H=12,7
s=(0,1); t=(11,4); k=1
def bres(a,b):
    (x0,y0),(x1,y1)=a,b; dx=abs(x1-x0); dy=abs(y1-y0); sx=1 if x1>x0 else -1; sy=1 if y1>y0 else -1; err=dx-dy; out=[]
    while True:
        out.append((x0,y0))
        if (x0,y0)==(x1,y1): break
        e2=2*err
        if e2>-dy: err-=dy; x0+=sx
        if e2<dx: err+=dx; y0+=sy
    return out
line=bres(s,t)
C={(x+i,y+j) for (x,y) in line for i in range(-k,k+1) for j in range(-k,k+1) if 0<=x+i<W and 0<=y+j<H}
blocked={(5,2),(5,3)}
C-=blocked
def inb(c): return 0<=c[0]<W and 0<=c[1]<H and c not in blocked
DIRS=[(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]
def legal(a,d,blk):
    b=(a[0]+d[0],a[1]+d[1])
    if not (0<=b[0]<W and 0<=b[1]<H) or b in blk: return False
    if d[0] and d[1]:
        if (a[0]+d[0],a[1]) in blk or (a[0],a[1]+d[1]) in blk: return False
    return True
X={ (a[0]+d[0],a[1]+d[1]) for a in C for d in DIRS if legal(a,d,blocked) and (a[0]+d[0],a[1]+d[1]) not in C}
# diagonal inside C whose side cell is outside C
diag=[(a,d) for a in C for d in DIRS if d[0] and d[1] and legal(a,d,blocked) and (a[0]+d[0],a[1]+d[1]) in C and (((a[0]+d[0],a[1]) not in C) or ((a[0],a[1]+d[1]) not in C))]
print(len(C),len(X),diag[:6])
# in-corridor optimal path (G[C])
def dij(src,nodes,blk,extra_block=set()):
    dist={src:0}; pq=[(0,src)]; par={}
    while pq:
        d,u=heapq.heappop(pq)
        if d>dist[u]: continue
        for dd in DIRS:
            if not legal(u,dd,blk|extra_block): continue
            v=(u[0]+dd[0],u[1]+dd[1])
            if v not in nodes: continue
            nd=d+(math.sqrt(2) if dd[0] and dd[1] else 1)
            if nd<dist.get(v,1e9)-1e-9: dist[v]=nd; par[v]=u; heapq.heappush(pq,(nd,v))
    return dist,par
distC,parC=dij(s,C,blocked)
path=[t]; 
while path[-1]!=s: path.append(parC[path[-1]])
path=path[::-1]
U=distC[t]
oct=lambda a,b: max(abs(a[0]-b[0]),abs(a[1]-b[1]))+(math.sqrt(2)-1)*min(abs(a[0]-b[0]),abs(a[1]-b[1]))
# exits generated: f = g(n)+d + h(b) for n in C with g=distC
EX={}
for a in C:
    for d in DIRS:
        if legal(a,d,blocked):
            b=(a[0]+d[0],a[1]+d[1])
            if b in X:
                g=distC[a]+(math.sqrt(2) if d[0] and d[1] else 1)
                if g<EX.get(b,(1e9,))[0]: EX[b]=(g,a)
Phi=min(g+oct(b,t) for b,(g,a) in EX.items()); bstar=min(EX,key=lambda b:EX[b][0]+oct(b,t))
print('U',U,'Phi',Phi,'b*',bstar,EX[bstar])
# whole map optimum
distG,_=dij(s,{(x,y) for x in range(W) for y in range(H)}-blocked,blocked)
print('C*',distG[t])
# write tikz
def cell(x,y,style): return f"\\fill[{style}] ({x},{y}) rectangle ({x+1},{y+1});\n"
out=[]
def panel(mc):
    o=[]
    for x in range(W):
        for y in range(H):
            c=(x,y)
            if c in blocked: o.append(cell(x,y,'black!70'))
            elif c in C: o.append(cell(x,y,'blue!12'))
            elif c in X and not mc: o.append(cell(x,y,'pattern=north east lines,pattern color=orange!80!black'))
            elif mc: o.append(cell(x,y,'black!25'))
    o.append(f"\\draw[black!40,very thin] (0,0) grid[step=1] ({W},{H});\n")
    # bresenham line
    o.append(f"\\draw[black!60,dashed,thin] ({s[0]+.5},{s[1]+.5}) -- ({t[0]+.5},{t[1]+.5});\n")
    # path
    pts=' -- '.join(f"({x+.5},{y+.5})" for x,y in path)
    o.append(f"\\draw[blue!60!black,very thick] {pts};\n")
    o.append(f"\\node[font=\\footnotesize] at ({s[0]+.5},{s[1]+.5}) {{$s$}};\\node[font=\\footnotesize] at ({t[0]+.5},{t[1]+.5}) {{$t$}};\n")
    return ''.join(o)
a,d=((6,1),(1,1)); b=(a[0]+d[0],a[1]+d[1]); side=[(a[0]+d[0],a[1]),(a[0],a[1]+d[1])]
sideout=[c for c in side if c not in C][0]
pa=panel(False)
g,n=EX[bstar]
pa+=f"\\draw[red!70!black,very thick,->] ({n[0]+.5},{n[1]+.5}) -- ({bstar[0]+.5},{bstar[1]+.5});\n"
pa+=f"\\draw[red!70!black,thick,dotted] ({bstar[0]+.5},{bstar[1]+.5}) -- ({t[0]+.5},{t[1]+.5});\n"
pa+=f"\\node[font=\\footnotesize,red!70!black,anchor=south] at ({bstar[0]+.5},{bstar[1]+1}) {{$b$}};\n"
pa+=f"\\draw[green!50!black,very thick,->] ({a[0]+.5},{a[1]+.5}) -- ({b[0]+.5},{b[1]+.5});\n"
outside={(x,y) for x in range(W) for y in range(H)}-C
distM,parM=dij(s,C,blocked,outside)
pm=[t]
while pm[-1]!=s: pm.append(parM[pm[-1]])
pm=pm[::-1]; print('U_MC',distM[t],pm)
pb=panel(True).replace(' -- '.join(f"({x+.5},{y+.5})" for x,y in path), ' -- '.join(f"({x+.5},{y+.5})" for x,y in pm))
pb+=f"\\draw[green!50!black,very thick] ({a[0]+.5},{a[1]+.5}) -- ({b[0]+.5},{b[1]+.5});\n\\draw[red,thick] ({(a[0]+b[0])/2+.5-.18},{(a[1]+b[1])/2+.5-.18}) -- ++(.36,.36) ({(a[0]+b[0])/2+.5-.18},{(a[1]+b[1])/2+.5+.18}) -- ++(.36,-.36);\n"
# boundary cell c: pick n (the parent of b*) as boundary cell; scan ray through it: horizontal ray from a cell to the left
cc=n
pb+=f"\\draw[orange!80!black,very thick,->] ({cc[0]-2+.5},{cc[1]+.5}) -- ({cc[0]+1.5},{cc[1]+.5});\n"
pb+=f"\\node[font=\\footnotesize,orange!80!black,anchor=north] at ({cc[0]+.5},{cc[1]}) {{$c$}};\n"
pb+=f"\\draw[red!70!black,thick,dotted] ({cc[0]+.5},{cc[1]+.5}) -- ({bstar[0]+.5},{bstar[1]+.5}) -- ({t[0]+.5},{t[1]+.5});\n"
tex=r"""\begin{tikzpicture}[x=0.3cm,y=0.3cm]
\begin{scope}
"""+pa+r"""\node[anchor=north,font=\small] at (6,-0.4) {(a) $G_C^+$: corridor A*};
\end{scope}
\begin{scope}[xshift=4.35cm]
"""+pb+r"""\node[anchor=north,font=\small] at (6,-0.4) {(b) $\MC$: JPS$_C$};
\end{scope}
\end{tikzpicture}"""
open('fig_corridor.tex','w').write(tex)
print('b*',bstar,'parent',n,'diag',a,b,'sideout',sideout,'path',path)
