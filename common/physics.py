"""Numerical conventions shared by Project B steps (H=k=1, nats/real sector)."""
from pathlib import Path
import json
import platform
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OMEGA=np.kron(np.eye(2),np.array([[0.,1.],[-1.,0.]]))
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})

def entropy(nu):
    nu=np.asarray(nu,dtype=float)
    if np.min(nu)<.5-1e-7: raise ArithmeticError('Unphysical reduced covariance')
    n=np.maximum(nu-.5,0)
    return (n+1)*np.log1p(n)-np.where(n>0,n*np.log(np.maximum(n,1e-300)),0)

def oscillator(c):
    K=np.array([[1.,c],[c,2.]])
    w,U=np.linalg.eigh(K)
    if np.min(w)<=0: raise ValueError('K must be positive')
    X=(U*(.5/np.sqrt(w)))@U.T
    P=(U*(.5*np.sqrt(w)))@U.T
    nu=np.sqrt(X[0,0]*P[0,0])
    return X,P,float(nu),float(entropy(nu))

def bd(x,M):
    mp.mp.dps=35
    z=mp.mpf(float(x));nu=mp.sqrt(mp.mpf('2.25')-M*M)
    pref=mp.sqrt(mp.pi)/2*mp.exp(mp.j*mp.pi*(nu/2+mp.mpf('.25')))
    h=mp.hankel1(nu,z)
    dh=(mp.hankel1(nu-1,z)-mp.hankel1(nu+1,z))/2
    sig=pref*z**mp.mpf('1.5')*h
    ds=pref*(mp.mpf('1.5')*mp.sqrt(z)*h+z**mp.mpf('1.5')*dh)
    q=(1-1j*x)*np.exp(1j*x)/np.sqrt(2)
    pq=-np.exp(1j*x)/(np.sqrt(2)*x)
    return q,pq,complex(sig),complex(-ds/z**2)

def initial(X,M,g,state='product',order=6):
    if state=='adiabatic':
        I=np.eye(2);J=np.array([[0.,1.],[-1.,0.]])
        D=np.array([[0.,-3*g],[0.,-M*M]])
        coeff={1:I+g*J/2}
        for n in range(2,order+1):
            cross=sum((coeff[j]@coeff[n-j] for j in range(1,n)),start=np.zeros((2,2),complex))
            coeff[n]=(((n+1)*I+g*J)@coeff[n-1]-cross+(D if n==2 else 0))/(2j)
        R=1j*I+sum((coeff[n]/X**n for n in coeff),start=np.zeros((2,2),complex))
        A=-R/X**2+g*np.array([[0.,1.],[0.,0.]])/X**3
        assert np.max(abs(A-A.T))<1e-13
        C=.5*X**2*np.linalg.inv(R.imag)
        w,U=np.linalg.eigh(C)
        if np.min(w)<=0:raise ValueError('Adiabatic covariance not positive')
        F=(U*np.sqrt(w))@U.T;vel=-X*R@F
        return np.array([F[0],vel[0],F[1],vel[1]])
    q,pq,s,ps=bd(X,M)
    u=np.zeros((4,2),complex)
    u[:,0]=[q,X**3*pq,0,0]
    u[:,1]=[0,-g*s,s,X**3*ps]
    return u

def rhs_N(N,y,M,g):
    q,v,s,w=y.reshape(4,2);x=np.exp(-N)
    return np.array([v,-3*v-x*x*q-g*(w+3*s),w,-3*w-(x*x+M*M)*s+g*v]).ravel()

def rhs_x(x,y,M,g):
    q,p,s,t=y.reshape(4,2)
    return np.array([g*s/x-x*x*p,q/x**2,-x*x*t,s/x**2+(M*M+g*g)*s/x**4-g*p/x]).ravel()

def modes(M=2.,g=.1,X=400.,xf=1e-4,rtol=1e-10,atol=1e-12,form='N',state='product',order=6):
    u=initial(X,M,g,state,order)
    if form=='N':
        sol=solve_ivp(lambda n,y:rhs_N(n,y,M,g),(-np.log(X),-np.log(xf)),u.ravel(),method='DOP853',rtol=rtol,atol=atol,dense_output=True)
    else:
        c=u.copy();c[1]=(u[1]+g*u[2])/X**3;c[3]=u[3]/X**3
        sol=solve_ivp(lambda x,y:rhs_x(x,y,M,g),(X,xf),c.ravel(),method='DOP853',rtol=rtol,atol=atol,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return {'sol':sol,'M':M,'g':g,'X':X,'xf':xf,'rtol':rtol,'atol':atol,'form':form,'state':state,'adiabatic_order':order}

def evaluate(m,x,canonical=True):
    x=np.atleast_1d(np.asarray(x,float))
    if np.min(x)<m['xf']*(1-1e-10) or np.max(x)>m['X']*(1+1e-10):raise ValueError('Outside integrated mode interval')
    y=m['sol'].sol(-np.log(x) if m['form']=='N' else x).reshape(4,2,-1)
    if m['form']=='N' and canonical:
        y=y.copy();y[1]=(y[1]+m['g']*y[2])/x**3;y[3]=y[3]/x**3
    return y

def diagnostics(m,x):
    x=np.atleast_1d(x);u=evaluate(m,x)
    out=[]
    for j,z in enumerate(x):
        a=u[:,:,j].astype(np.clongdouble)
        raw=a@a.conj().T
        comm=float(np.max(np.abs(2*raw.imag-OMEGA)))
        a[1]+=a[0]/z
        # Local canonical scalings balance the full covariance for eigensolving.
        v=(a@a.conj().T).real
        scales=np.ones(4,dtype=np.longdouble)
        for k in [0,2]:
            b=(v[k+1,k+1]/v[k,k])**.25;scales[k]=b;scales[k+1]=1/b
        vb=np.asarray(v*scales[:,None]*scales[None,:],float)
        eig=np.sort(np.abs(np.linalg.eigvals(1j*OMEGA@vb)))[::2]
        det=v[0,0]*v[1,1]-v[0,1]**2
        nu=float(np.sqrt(max(det,0)))
        S=float(entropy(nu)) if nu>=.5-1e-7 else np.nan
        out.append([z,float(raw[0,0].real),nu,S,comm,float(np.max(np.abs(eig-.5)))])
    return np.array(out)

def output(step):
    p=ROOT/f'step{step}';(p/'results').mkdir(parents=True,exist_ok=True);(p/'figures').mkdir(exist_ok=True)
    return p

def save_table(path,header,rows):
    np.savetxt(path,np.asarray(rows),delimiter=',',header=header,comments='',fmt='%.14g')

def save_json(path,data):
    path.write_text(json.dumps(data,indent=2,default=lambda x:float(x))+'\n')

def figure(p,name):
    from io import BytesIO
    for extension in ['png','pdf']:
        buf=BytesIO()
        plt.savefig(buf,format=extension,dpi=160)
        data=buf.getvalue()
        if extension=='pdf' and b'%%EOF' not in data[-50:]:raise IOError('Incomplete PDF export')
        target=p/'figures'/f'{name}.{extension}'
        temporary=target.with_suffix(target.suffix+'.part')
        temporary.write_bytes(data);temporary.replace(target)
    plt.close()


def save_modes(path,m,x):
    u=evaluate(m,x)
    np.savez_compressed(path,x=x,modes=u,M=m['M'],g=m['g'],X=m['X'],rtol=m['rtol'],atol=m['atol'],initial_state=m['state'],adiabatic_order=m['adiabatic_order'],rows='Q,Pi_q,Sigma,Pi_sigma')

def software():
    return {'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__,'matplotlib':matplotlib.__version__}
