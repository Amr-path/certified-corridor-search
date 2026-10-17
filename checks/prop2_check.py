import math
SQ2=math.sqrt(2)
def oct(a,b):
    dx,dy=abs(a[0]-b[0]),abs(a[1]-b[1]); return max(dx,dy)+(SQ2-1)*min(dx,dy)
bad=0; bad2=0; n=0
for l in range(1,80):
    for k in range(0,40):
        n+=1
        pred=oct((0,0),(l,2*(k+1)))
        best=min(oct((0,0),(x,y))+oct((x,y),(l,0)) for x in range(-100,l+101) for y in (k+1,-(k+1)))
        if abs(best-pred)>1e-9: bad+=1
        # also over all |y|>=k+1 up to 3(k+1)+5
        best2=min(oct((0,0),(x,y))+oct((x,y),(l,0)) for x in range(-100,l+101) for yy in range(k+1,3*(k+1)+6) for y in (yy,-yy))
        if abs(best2-pred)>1e-9: bad2+=1
print("pairs",n,"mismatch(|y|=k+1)",bad,"mismatch(|y|>=k+1)",bad2)
# norm check: oct(v) = (sqrt2-1)|v|_1 + (2-sqrt2)|v|_inf
import random
for _ in range(10000):
    x,y=random.randint(-50,50),random.randint(-50,50)
    assert abs(oct((0,0),(x,y)) - ((SQ2-1)*(abs(x)+abs(y)) + (2-SQ2)*max(abs(x),abs(y))))<1e-9
print("norm identity ok")
# design rule check: U=l (straight path) certified at eps=0 for every k?  l <= Lambda_k always
print("l<=Lambda_k for all pairs:", all(l <= oct((0,0),(l,2*(k+1)))+1e-12 for l in range(1,80) for k in range(40)))
# counterexample to 'any inner algorithm': uniform free grid l=10,k=5, DFS-like in-band path cost 15
l,k=10,5; Lam=oct((0,0),(l,2*(k+1))); print("l=10,k=5 Lambda=",round(Lam,3),"U_dfs=15 sp=10 -> min(U,Lambda)=",min(15,Lam),"> sp")
