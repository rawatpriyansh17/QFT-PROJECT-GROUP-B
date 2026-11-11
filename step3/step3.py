"""Symbolic Legendre transform and comparison of two independent ODE forms."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *
import sympy as sp

SETTINGS={'M': 2, 'g': 0.1, 'X': 100, 'xf': 0.001, 'rtol': 1e-10, 'atol': 1e-12, 'formulations': ['canonical_x', 'velocity_N']}

def main():
    p=output(3);a,k,m,r=sp.symbols('a k m rho',positive=True)
    q,s,v,w,P,T=sp.symbols('q sigma v w pi_q pi_sigma')
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
    # Repeat the comparison with tolerances tightened by a factor of 10
    n_tight = modes(
        M=2, g=0.1, X=100, xf=1e-3,
        rtol=1e-11, atol=1e-13, form='N'
    )

    c_tight = modes(
        M=2, g=0.1, X=100, xf=1e-3,
        rtol=1e-11, atol=1e-13, form='x'
    )

    # Evaluate both tighter solutions on the same x grid
    un_tight = evaluate(n_tight, x)
    uc_tight = evaluate(c_tight, x)

    # Maximum scaled difference at tighter tolerances
    rel_tight = np.max(
        np.abs(un_tight - uc_tight) / (1 + np.abs(un_tight)),
        axis=(0,1)
    )

    # Compare entanglement entropy at tighter tolerances
    dn_tight = diagnostics(n_tight, [1,.1,.01,.001])
    dc_tight = diagnostics(c_tight, [1,.1,.01,.001])

    entropy_diff_tight = np.max(
        np.abs(dn_tight[:,3] - dc_tight[:,3])
    )

    print("Original maximum column difference:", np.max(rel))
    print("Tight maximum column difference:", np.max(rel_tight))
    print("Tight entropy difference:", entropy_diff_tight)
    
    save_table(p/'results/formulation_agreement.csv','x,power_N,power_x,entropy_N,entropy_x',np.c_[dn[:,0],dn[:,1],dc[:,1],dn[:,3],dc[:,3]])
    plt.figure(figsize=(6.4,3.4))
    plt.loglog(x, rel, label='Original tolerances')
    plt.loglog(x, rel_tight, label='Tighter tolerances')
    plt.gca().invert_xaxis()
    plt.xlabel('x')
    plt.ylabel('Max scaled column difference')
    plt.legend()
    plt.tight_layout()
    figure(p,'formulation_agreement')

    save_json(p/'results/summary.json',{
        'settings': SETTINGS,
        'symbolic_legendre_check': True,

        'max_scaled_mode_difference': float(np.max(rel)),
        'max_entropy_difference': float(np.max(np.abs(dn[:,3]-dc[:,3]))),

        'tight_settings': {
            'rtol': 1e-11,
            'atol': 1e-13
        },
        'max_scaled_mode_difference_tight': float(np.max(rel_tight)),
        'max_entropy_difference_tight': float(entropy_diff_tight),

        'software': software()
    })
    
if __name__=='__main__':main()
