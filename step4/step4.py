"""Free Bunch-Davies modes, commutators and finite-time initial-state tests."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *

SETTINGS={'M': 2, 'g_free': 0, 'g_initial_time_test': 0.1, 'X_tests': [100, 200, 400, 800], 'product_tolerances': [[1e-10, 1e-12], [1e-11, 1e-13]], 'adiabatic_orders': [4, 6], 'adiabatic_tolerance': [2.5e-14, 1e-17]}

def main():
    p=output(4);xs=np.geomspace(100,1e-3,81)
    free=modes(g=0,X=100,xf=1e-3);u=evaluate(free,xs)
    ana=np.array([bd(float(x),2.) for x in xs]).T
    qe=np.abs(u[0,0]-ana[0])/(1+np.abs(ana[0]));se=np.abs(u[2,1]-ana[2])/(1+np.abs(ana[2]))
    fd=diagnostics(free,xs);rows=[];hist=[]
    for X in [100.,200.,400.,800.]:
        for tight in [False,True]:
            mm=modes(X=X,xf=1e-3,rtol=1e-11 if tight else 1e-10,atol=1e-13 if tight else 1e-12)
            d=diagnostics(mm,[1,.1,.01,.001])
            for z,powr,nu,S,comm,pure in d:rows.append([X,int(tight),z,powr,S,comm,pure])
            if tight:hist.append([X,*d[:,3]])
    save_table(p/'results/initial_time_checks.csv','X,tight_tolerance,x,power,entropy,commutator_residual,purity_error',rows)
    save_table(p/'results/free_mode_errors.csv','x,light_scaled_error,heavy_scaled_error,commutator_residual',np.c_[xs,qe,se,fd[:,4]])
    fig,ax=plt.subplots(1,2,figsize=(7,3));ax[0].loglog(xs,np.abs(ana[0]),label='Analytic Q');ax[0].loglog(xs,np.abs(u[0,0]),'--',label='Numerical Q');ax[0].loglog(xs,np.abs(ana[2]),label='Analytic Sigma');ax[0].loglog(xs,np.abs(u[2,1]),'--',label='Numerical Sigma');ax[0].invert_xaxis();ax[0].legend(fontsize=8);ax[0].set(xlabel='x',ylabel='Mode magnitude')
    ax[1].loglog(xs,qe,label='Light');ax[1].loglog(xs,se,label='Heavy');ax[1].invert_xaxis();ax[1].legend();ax[1].set(xlabel='x',ylabel='Scaled complex error');figure(p,'free_modes')
    hh=np.array(hist);plt.figure(figsize=(6.5,3.4))
    for j,z in enumerate([1,.1,.01,.001]):plt.plot(hh[:,0],hh[:,j+1],'-o',label=f'x={z:g}')
    plt.xlabel('Initial X');plt.ylabel('Entropy per real sector (nats)');plt.legend();figure(p,'initial_time_convergence')
    adrows=[]
    for X in [100.,200.,400.,800.]:
        for order in [4,6]:
            ad=modes(X=X,xf=1e-3,rtol=2.5e-14,atol=1e-17,state='adiabatic',order=order)
            dd=diagnostics(ad,[.01,.001])
            for zz,pw,nv,en,co,pu in dd:adrows.append([X,order,zz,pw,en,co,pu])
    save_table(p/'results/adiabatic_initial_state.csv','X,order,x,power,entropy,commutator_residual,purity_error',adrows)
    save_modes(p/'results/free_complex_modes.npz',free,xs)
    save_json(p/'results/summary.json',{'settings':SETTINGS,'max_light_error':float(max(qe)),'max_heavy_error':float(max(se)),'late_free_power':float(fd[-1,1]),'max_free_commutator_residual':float(max(fd[:,4])),'entropy_X400_vs_800':(hh[-1,1:]-hh[-2,1:]).tolist(),'software':software()})
if __name__=='__main__':main()
