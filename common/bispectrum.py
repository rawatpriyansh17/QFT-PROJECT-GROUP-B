"""Leading spectator bispectrum with common final time and regulated quadrature."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from .physics import evaluate

def mesh(zf,zmax,r,refine=1,order=12):
    t,w=leggauss(order)
    lo=np.linspace(np.log(zf),0,int(np.ceil(-np.log(zf)/(.08/refine)))+1)
    hi=np.linspace(1,zmax,int(np.ceil((zmax-1)/(np.pi/(4*(2+r))/refine)))+1)
    z=[];weights=[]
    for edges,logarithmic in [(lo,True),(hi,False)]:
        left=edges[:-1,None];right=edges[1:,None]
        x=(left+right)/2+(right-left)/2*t
        ws=(right-left)/2*w
        if logarithmic:ws=ws*np.exp(x);x=np.exp(x)
        z.append(x.ravel());weights.append(ws.ravel())
    return np.concatenate(z),np.concatenate(weights)

def integral(m,r=.03,zf=.01,zmax=320.,epsilons=(.1,.05,.025),lam=.01,refine=1,order=12):
    z,w=mesh(zf,zmax,r,refine,order)
    qf=evaluate(m,[zf])[0,:,0].astype(np.clongdouble)
    ql=evaluate(m,[r*zf])[0,:,0].astype(np.clongdouble)
    sig=evaluate(m,z)[2].astype(np.clongdouble)
    sigl=evaluate(m,r*z)[2].astype(np.clongdouble)
    cs=np.sum(qf[:,None]*sig.conj(),axis=0)
    cl=np.sum(ql[:,None]*sigl.conj(),axis=0)
    complex_integrand=cl*cs*cs/z**4
    pl=float(np.sum(np.abs(ql)**2));ps=float(np.sum(np.abs(qf)**2))
    values=[]
    for eps in epsilons:
        damp=np.exp(-eps*(2+r)*z)
        value=2*lam*np.sum(w*damp*complex_integrand.imag,dtype=np.longdouble)/(pl*ps)
        values.append(float(value))
    return np.array(values),{'z':z,'weights':w,'integrand':np.asarray(complex_integrand,complex),'power_long_dimensionless':pl,'power_short_dimensionless':ps}

def extrapolate(eps,values):
    eps=np.asarray(eps);values=np.asarray(values)
    # Same last three points for the two extrapolation orders.
    lin=np.polyfit(eps[-3:],values[-3:],1)[-1]
    quad=np.polyfit(eps[-3:],values[-3:],2)[-1]
    return float(quad),float(abs(quad-lin))

def adaptive_check(m,r=.03,zf=.01,zmax=640,epsilon=.025,lam=.01):
    from scipy.integrate import quad
    qf=evaluate(m,[zf])[0,:,0].astype(np.clongdouble)
    ql=evaluate(m,[r*zf])[0,:,0].astype(np.clongdouble)
    norm=2*lam/(float(sum(abs(qf)**2))*float(sum(abs(ql)**2)))
    def f(z):
        s=evaluate(m,[z])[2,:,0].astype(np.clongdouble)
        sl=evaluate(m,[r*z])[2,:,0].astype(np.clongdouble)
        cs=sum(qf*s.conj());cl=sum(ql*sl.conj())
        return float((cl*cs*cs).imag/z**4*np.exp(-epsilon*(2+r)*z))
    val,err=quad(lambda y:np.exp(y)*f(np.exp(y)),np.log(zf),0,epsabs=1e-15,epsrel=1e-8,limit=200)
    edges=np.linspace(1,zmax,int(np.ceil((zmax-1)/(np.pi/(4*(2+r)))))+1)
    for a,b in zip(edges[:-1],edges[1:]):
        v,e=quad(f,a,b,epsabs=1e-16,epsrel=1e-8,limit=100);val+=v;err+=e
    return norm*val,norm*err
