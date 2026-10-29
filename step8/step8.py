"""Separate initial-state, regulator, quadrature and final-time errors."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *
from common.bispectrum import integral,extrapolate,adaptive_check

def product_control(p):
    rows=[]
    for X in [800,1600]:
        for state in ['product','adiabatic']:
            m=modes(X=X,xf=1e-5,rtol=1e-12,atol=1e-15,state=state)
            vals,_=integral(m,.03,zmax=640,epsilons=[.1,.05,.025,.0125])
            f,e=extrapolate([.1,.05,.025,.0125],vals)
            rows.append([X,int(state=='adiabatic'),f,e])
    save_table(p/'results/initial_state_control.csv','X,adiabatic,F_extrapolated,regulator_error',rows)

def independent_quadrature(p):
    m=modes(X=1600,xf=1e-5,rtol=1e-12,atol=1e-15,state='adiabatic')
    v,_=integral(m,.03,zmax=640,epsilons=[.025])
    a,e=adaptive_check(m)
    save_table(p/'results/adaptive_quadrature.csv','Gauss_Legendre,adaptive_Gauss_Kronrod,difference,adaptive_error_estimate',[[v[0],a,abs(v[0]-a),e]])

SETTINGS={'M': 2, 'g': 0.1, 'lambda3': 0.01, 'X': [800, 1600], 'xf': 1e-06, 'epsilons': [0.1, 0.05, 0.025, 0.0125], 'zf': [0.01, 0.003, 0.001], 'zmax_full': 640, 'initial_state': 'adiabatic_order_6', 'rtol': 1e-12, 'atol': 1e-15, 'quadrature_Gauss_orders': [12, 16]}

def main():
    p=output(8);eps=np.array([.1,.05,.025,.0125]);rs=[.2,.1,.03,.01,.003,.001]
    baseline=modes(X=1600,xf=1e-6,rtol=1e-12,atol=1e-15,state='adiabatic')
    other=modes(X=800,xf=1e-6,rtol=1e-12,atol=1e-15,state='adiabatic')
    loose=modes(X=1600,xf=1e-6,rtol=1e-11,atol=1e-14,state='adiabatic')
    rows=[];summary=[];variation_rows=[]
    for r in rs:
        vals,_=integral(baseline,r,zmax=640,epsilons=eps)
        variation_rows.extend([[r,0,e,v] for e,v in zip(eps,vals)])
        fixed=[]
        for e,L in zip([.1,.05,.025],[80,160,320]):
            short,_=integral(baseline,r,zmax=L,epsilons=[e]);long,_=integral(baseline,r,zmax=2*L,epsilons=[e])
            fixed.append(abs(short[0]-long[0]));rows.append([r,0,e,L,short[0],long[0]])
        central,regerr=extrapolate(eps,vals)
        vx,_=integral(other,r,zmax=640,epsilons=eps);fx,_=extrapolate(eps,vx)
        vt,_=integral(loose,r,zmax=640,epsilons=eps);ft,_=extrapolate(eps,vt)
        vq,_=integral(baseline,r,zmax=640,epsilons=eps,refine=2,order=16);fq,_=extrapolate(eps,vq)
        variation_rows.extend([[r,1,e,v] for e,v in zip(eps,vx)])
        variation_rows.extend([[r,2,e,v] for e,v in zip(eps,vt)])
        variation_rows.extend([[r,3,e,v] for e,v in zip(eps,vq)])
        finals=[]
        for final_index,zf in enumerate([.003,.001]):
            vf,_=integral(baseline,r,zf=zf,zmax=640,epsilons=eps);ff,_=extrapolate(eps,vf);finals.append(ff)
            variation_rows.extend([[r,4+final_index,e,v] for e,v in zip(eps,vf)])
        errs=[abs(fx-central),regerr,abs(ft-central),abs(fq-central),max(abs(np.array(finals)-central)),max(fixed)]
        uncertainty=sum(errs)
        summary.append([r,central,*errs,uncertainty])
        for e,v in zip(eps,vals):rows.append([r,1,e,640,v,v])
    save_table(p/'results/variation_integrals.csv','r,variation_code,epsilon,F_regulated',variation_rows)
    save_table(p/'results/cutoff_and_regulator.csv','r,test_type,epsilon,zmax,F_first,F_doubled_cutoff',rows)
    save_table(p/'results/error_budget.csv','r,F_extrapolated,initial_X,regulator,ODE_tolerance,quadrature,final_time,early_tail,total_conservative_error',summary)
    summary=np.array(summary)
    plt.figure(figsize=(6.5,3.5));plt.errorbar(summary[:,0],summary[:,1],yerr=summary[:,-1],fmt='o',capsize=4);plt.xscale('log');plt.xlabel('r = kL/kS');plt.ylabel('Extrapolated F with numerical error budget');figure(p,'bispectrum_convergence')
    # Only a convergence-qualified pilot can justify a physical frequency scan.
    envelope=max(abs(summary[:,1])/summary[:,0]**1.5)
    qualified=bool(np.all(summary[:,-1]/summary[:,0]**1.5<.01*envelope))
    save_json(p/'results/summary.json',{'settings':SETTINGS,'initial_state':'sixth-order ultraviolet Riccati expansion','pilot_relative_errors':(summary[:,-1]/np.maximum(abs(summary[:,1]),1e-30)).tolist(),'pilot_below_one_percent_of_clock_envelope':qualified,'clock_detection_gate':'The clock component must exceed 10 times the numerical uncertainty; total-signal precision alone is insufficient.','pilot_clock_envelope_error':float(max(summary[:,-1]/summary[:,0]**1.5)/envelope),'software':software()})
if __name__=='__main__':
    main()
    product_control(output(8))
    independent_quadrature(output(8))
