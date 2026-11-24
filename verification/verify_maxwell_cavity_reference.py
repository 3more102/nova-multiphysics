"""Independent SciPy PEC cube cavity reference, not NOVA C++ execution."""
import numpy as np
from scipy.linalg import eigh
from verify_maxwell_hcurl_reference import mesh, EDGES, L, weights

def cavity(n):
    x,tets=mesh(n)
    edges=sorted({tuple(sorted((t[i],t[j]))) for t in tets for i,j in EDGES})
    free=[e for e in edges if not any(x[e[0],d]==x[e[1],d] and x[e[0],d] in (0,1) for d in range(3))]
    ids={e:i for i,e in enumerate(free)}
    K=np.zeros((len(free),len(free))); M=K.copy()
    for t in tets:
        p=x[t]; J=(p[1:]-p[0]).T; V=abs(np.linalg.det(J))/6
        g=np.vstack((-np.ones(3),np.eye(3)))@np.linalg.inv(J)
        signs=[1 if t[i]<t[j] else -1 for i,j in EDGES]
        curls=np.array([2*s*np.cross(g[i],g[j]) for (i,j),s in zip(EDGES,signs)])
        B=np.stack([s*(L[...,i,None]*g[j]-L[...,j,None]*g[i]) for (i,j),s in zip(EDGES,signs)],axis=-2)
        m=V*np.einsum('q,qid,qjd->ij',weights,B,B); k=V*curls@curls.T
        loc=[ids.get(tuple(sorted((t[i],t[j]))),-1) for i,j in EDGES]
        for a,i in enumerate(loc):
            if i<0: continue
            for b,j in enumerate(loc):
                if j>=0: K[i,j]+=k[a,b]; M[i,j]+=m[a,b]
    val,_=eigh(K,M)
    zero=abs(val)<1e-8; physical=val[~zero][:3]
    assert zero.sum()==(n-1)**3
    print(n,len(free),zero.sum(),physical,'relative frequency errors',
          np.sqrt(physical/(2*np.pi**2))-1,flush=True)
    return physical
if __name__=='__main__':
    for n in (2,3,4): cavity(n)
