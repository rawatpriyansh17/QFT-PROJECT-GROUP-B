"""Profile likelihood for a logarithmic clock with linear nuisance coefficients."""
import numpy as np
from scipy.optimize import minimize_scalar,brentq

def design(r,mu,extra=False):
    cols=[np.ones(len(r)),r*r,r**1.5*np.cos(mu*np.log(r)),r**1.5*np.sin(mu*np.log(r))]
    if extra:cols.append(r**4)
    return np.array(cols).T

def profile(r,y,sigma,extra=False):
    r=np.asarray(r);y=np.asarray(y);sigma=np.broadcast_to(sigma,y.shape)
    grid=np.arange(.1,4.0001,.01)
    def solve(mu):
        A=design(r,mu,extra)/sigma[:,None];b=y/sigma
        c=np.linalg.lstsq(A,b,rcond=None)[0]
        return float(np.sum((A@c-b)**2)),c
    chi=np.array([solve(mu)[0] for mu in grid]);candidates=[]
    for j in range(len(grid)):
        if (j==0 or chi[j]<=chi[j-1]) and (j==len(grid)-1 or chi[j]<=chi[j+1]):
            a=grid[max(0,j-1)];b=grid[min(len(grid)-1,j+1)]
            opt=minimize_scalar(lambda m:solve(m)[0],bounds=(a,b),method='bounded',options={'xatol':1e-10})
            candidates.append((float(opt.fun),float(opt.x)))
    best,mu=min(candidates);coeff=solve(mu)[1]
    fine=np.unique(np.r_[grid,mu,np.linspace(max(.1,mu-.05),min(4,mu+.05),101)])
    values=np.array([solve(m)[0] for m in fine]);allowed=values<=best+3.841459
    intervals=[];start=None
    for j in range(len(fine)):
        if allowed[j] and start is None:
            start=fine[j] if j==0 else brentq(lambda m:solve(m)[0]-best-3.841459,fine[j-1],fine[j])
        if start is not None and (not allowed[j] or j==len(fine)-1):
            end=fine[j] if allowed[j] else brentq(lambda m:solve(m)[0]-best-3.841459,fine[j-1],fine[j])
            intervals.append([float(start),float(end)]);start=None
    return {'mu':mu,'chi2':best,'coefficients':coeff.tolist(),'intervals_95':intervals,'competitive_minima':[[m,c-best] for c,m in sorted(candidates) if c-best<10],'grid':fine,'delta_chi2':values-best}
