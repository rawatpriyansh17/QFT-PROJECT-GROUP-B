"""Symbolic Legendre transform and comparison of two independent ODE forms."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *
import sympy as sp

SETTINGS={'M': 2, 'g': 0.1, 'X': 100, 'xf': 0.001, 'rtol': 1e-10, 'atol': 1e-12, 'formulations': ['canonical_x', 'velocity_N']}

def main():
    p=output(3);a,k,m,r=sp.symbols('a k m rho',positive=True)
    q,s,v,w,P,T=sp.symbols('q sigma v w P T')
    L=a*a*(v*v+w*w-k*k*(q*q+s*s))/2-a**4*m*m*s*s/2+a**3*r*v*s
    vel={v:(P-a**3*r*s)/a**2,w:T/a**2}
    H=sp.expand((P*v+T*w-L).subs(vel))
    expected=((P-a**3*r*s)**2+T*T)/(2*a*a)+a*a*k*k*(q*q+s*s)/2+a**4*m*m*s*s/2
    assert sp.simplify(H-expected)==0
    eqs=[sp.diff(expected,P),-sp.diff(expected,q),sp.diff(expected,T),-sp.diff(expected,s)]
    (p/'results/hamilton_equations.txt').write_text('Hamiltonian: '+str(H)+'\nHamilton equations: '+str(eqs)+'\n')
    x=np.geomspace(100,1e-3,100)
    n=modes(X=100,xf=1e-3);c=modes(X=100,xf=1e-3,form='x')
    un=evaluate(n,x);uc=evaluate(c,x)
    rel=np.max(np.abs(un-uc)/(1+np.abs(un)),axis=(0,1))
    dn=diagnostics(n,[1,.1,.01,.001]);dc=diagnostics(c,[1,.1,.01,.001])
    save_table(p/'results/formulation_agreement.csv','x,power_N,power_x,entropy_N,entropy_x',np.c_[dn[:,0],dn[:,1],dc[:,1],dn[:,3],dc[:,3]])
    plt.figure(figsize=(6.4,3.4));plt.loglog(x,rel);plt.gca().invert_xaxis();plt.xlabel('x');plt.ylabel('Max scaled column difference');figure(p,'formulation_agreement')
    save_json(p/'results/summary.json',{'settings':SETTINGS,'symbolic_legendre_check':True,'max_scaled_mode_difference':float(max(rel)),'max_entropy_difference':float(max(abs(dn[:,3]-dc[:,3]))),'software':software()})
if __name__=='__main__':main()
