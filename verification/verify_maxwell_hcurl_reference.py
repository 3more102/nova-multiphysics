"""Independent NumPy/SciPy reference, not execution of NOVA C++ source."""
import itertools
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from numpy.polynomial.legendre import leggauss

EDGES = list(itertools.combinations(range(4),2))
LOCAL = [(0,1,3,7),(0,3,2,7),(0,2,6,7),(0,6,4,7),(0,4,5,7),(0,5,1,7)]
Q, W = leggauss(5); Q=(Q+1)/2; W=W/2
ABC=np.array(list(itertools.product(Q,repeat=3)))
ww=np.array([a*b*c for a,b,c in itertools.product(W,repeat=3)])
a,b,c=ABC.T
L=np.column_stack(((1-a)*(1-b)*(1-c),a,(1-a)*b,(1-a)*(1-b)*c))
weights=6*ww*(1-a)**2*(1-b)
EPS=2+0.6j; K0=2.3

def mesh(n):
    x=np.array([(i/n,j/n,k/n) for k in range(n+1) for j in range(n+1) for i in range(n+1)])
    idx=lambda i,j,k:i+(n+1)*(j+(n+1)*k)
    tet=[]
    for k,j,i in itertools.product(range(n),repeat=3):
        v=[idx(i+dx,j+dy,k+dz) for dz,dy,dx in itertools.product(range(2),repeat=3)]
        tet.extend([[v[t] for t in cell] for cell in LOCAL])
    return x,np.array(tet)

def run(n):
    x,tets=mesh(n)
    ed=sorted({tuple(sorted((t[i],t[j]))) for t in tets for i,j in EDGES})
    free=[e for e in ed if not any(x[e[0],d]==x[e[1],d] and x[e[0],d] in (0,1) for d in range(3))]
    ids={e:i for i,e in enumerate(free)}
    rows=[];cols=[];vals=[];rhs=np.zeros(len(free),complex); cache=[]
    exact=lambda p: np.sin(np.pi*p[...,0])*np.sin(np.pi*p[...,2])
    for t in tets:
        p=x[t]; J=(p[1:]-p[0]).T; V=abs(np.linalg.det(J))/6
        g=np.vstack((-np.ones(3),np.eye(3)))@np.linalg.inv(J)
        signs=np.array([1 if t[i]<t[j] else -1 for i,j in EDGES])
        curls=np.array([2*s*np.cross(g[i],g[j]) for (i,j),s in zip(EDGES,signs)])
        basis=lambda l: np.stack([s*(l[...,i,None]*g[j]-l[...,j,None]*g[i]) for (i,j),s in zip(EDGES,signs)],axis=-2)
        B=basis(L); M=V*np.einsum('q,qid,qjd->ij',weights,B,B)
        A=V*curls@curls.T-K0**2*EPS*M
        F=(2*np.pi**2-K0**2*EPS)*exact(p.mean(axis=0))
        f=V*basis(np.full(4,0.25))[:,1]*F
        local=[ids.get(tuple(sorted((t[i],t[j]))),-1) for i,j in EDGES]
        for i,ii in enumerate(local):
            if ii<0:continue
            rhs[ii]+=f[i]
            for j,jj in enumerate(local):
                if jj>=0: rows.append(ii);cols.append(jj);vals.append(A[i,j])
        cache.append((p,V,B,local))
    A=coo_matrix((vals,(rows,cols)),shape=(len(free),len(free))).tocsr()
    sol=spsolve(A,rhs)
    err=0; norm=0
    for p,V,B,local in cache:
        coeff=np.array([sol[i] if i>=0 else 0 for i in local]); E=np.einsum('i,qid->qd',coeff,B)
        ref=exact(L@p); E[:,1]-=ref
        err+=V*np.dot(weights,np.sum(abs(E)**2,axis=1)); norm+=V*np.dot(weights,ref**2)
    print(f'n={n}, tetrahedra={len(tets)}, free_edges={len(free)}, relative_L2={np.sqrt(err/norm):.9g}, relative_residual={np.linalg.norm(A@sol-rhs)/np.linalg.norm(rhs):.3g}',flush=True)
    return np.sqrt(err/norm)
if __name__=='__main__':
    errs=[run(n) for n in (2,4,8)]
    print('Observed L2 orders:',np.log2(np.array(errs[:-1])/errs[1:]))
