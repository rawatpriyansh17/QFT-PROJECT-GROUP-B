"""Initial regulated spectator calculation; convergence is assessed in Step 8."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *
from common.bispectrum import integral

SETTINGS={'M': 2, 'g': 0.1, 'lambda3': 0.01, 'X': 800, 'xf': 1e-05, 'zf': 0.01, 'epsilon': 0.025, 'zmax': 320, 'initial_state': 'adiabatic_order_6', 'rtol': 1e-12, 'atol': 1e-15}

def main():
    p=output(7);rows=[];mm=modes(X=800,xf=1e-5,rtol=1e-12,atol=1e-15,state='adiabatic')
    for r in [.1,.03,.01]:
        vals,raw=integral(mm,r,zmax=320,epsilons=[.025])
        rows.append([r,vals[0]])
        if r==.03:
            np.savez_compressed(p/'results/complex_integrand.npz',**raw,r=r,epsilon=.025,lambda3=.01,zf=.01,X=800)
            cumulative=2*.01*np.cumsum(raw['weights']*np.exp(-.025*(2+r)*raw['z'])*raw['integrand'].imag)/(raw['power_long_dimensionless']*raw['power_short_dimensionless'])
            plt.figure(figsize=(6.6,3.5));plt.semilogx(raw['z'],cumulative);plt.xlabel('z = -kS eta');plt.ylabel('Cumulative regulated F');figure(p,'cumulative_integral')
    controls=[]
    for g in [0.,.025,.05,.1]:
        model=mm if g==.1 else modes(g=g,X=800,xf=1e-5,rtol=1e-12,atol=1e-15,state='adiabatic')
        val,_=integral(model,.03,zmax=320,epsilons=[.025]);controls.append([g,val[0],val[0]/g**3 if g else 0])
    save_table(p/'results/regulated_bispectrum.csv','r,regulated_F_epsilon_0025',rows)
    save_table(p/'results/mixing_controls.csv','g,F,F_over_g_cubed',controls)
    save_modes(p/'results/complex_modes.npz',mm,np.r_[np.linspace(800,1,7991),np.geomspace(.99,1e-5,501)])
    z,_=integral(mm,.03,zmax=320,epsilons=[.025],lam=0)
    twice,_=integral(mm,.03,zmax=320,epsilons=[.025],lam=.02)
    assert controls[0][1]==0 and z[0]==0
    save_json(p/'results/summary.json',{'settings':SETTINGS,'M':2.,'g':.1,'lambda3_over_H':.01,'state':'adiabatic_order_6','X':800,'zf':.01,'zmax':320,'epsilon':.025,'zero_controls_pass':True,'doubling_lambda_ratio':float(twice[0]/controls[-1][1]),'status':'Regulated finite-initial-time result; see Step 8 for convergence','software':software()})
if __name__=='__main__':main()
