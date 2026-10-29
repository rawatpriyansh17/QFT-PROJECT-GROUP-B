"""Gaussian entropy and power scan with explicit precision diagnostics."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *

SETTINGS={'M_values': [1.6, 1.8, 2, 2.4], 'g_values': [0, 0.05, 0.1, 0.2], 'x_comparisons': [1, 0.1, 0.01, 0.001], 'X': 400, 'initial_state': 'adiabatic_order_6', 'rtol': 2.5e-14, 'atol': 1e-17, 'commutator_limit': 1e-06, 'purity_limit': 1e-05}

def main():
    p=output(5);xs=np.geomspace(100,1e-4,220)
    mm=modes(X=400,xf=1e-4,rtol=2.5e-14,atol=1e-17,state='adiabatic');d=diagnostics(mm,xs)
    zero=diagnostics(modes(g=0,X=400,xf=1e-4,rtol=2.5e-14,atol=1e-17,state='adiabatic'),xs)
    good=np.logical_and.accumulate((d[:,4]<1e-6)&(d[:,5]<1e-5)&np.isfinite(d[:,3]))
    rows=[]
    for M in [1.6,1.8,2.,2.4]:
        for g in [0.,.05,.1,.2]:
            mod=modes(M,g,X=400,xf=1e-3,rtol=2.5e-14,atol=1e-17,state='adiabatic')
            ds=diagnostics(mod,[1,.1,.01,.001])
            for z,powr,nu,S,comm,pure in ds:rows.append([M,g,z,S,powr/(.5*(1+z*z))-1,comm,pure,int(comm<1e-6 and pure<1e-5 and np.isfinite(S))])
    rows=np.array(rows)
    save_table(p/'results/gaussian_scan.csv','M,g,x,entropy,relative_power_change,commutator_residual,purity_error,reliable',rows)
    save_table(p/'results/entropy_history.csv','x,power,nu,entropy,commutator_residual,purity_error',d)
    save_modes(p/'results/benchmark_complex_modes.npz',mm,xs)
    fig,ax=plt.subplots(1,2,figsize=(7.4,3.3));N=-np.log(xs)
    ax[0].plot(N[good],d[good,3],label='g=0.1');ax[0].plot(N,np.nan_to_num(zero[:,3]),'--',label='g=0');ax[0].axvline(0,color='gray',ls=':');ax[0].set(xlabel='N = -ln x',ylabel='Entropy (nats/real sector)');ax[0].legend()
    ax[1].semilogy(N,np.maximum(d[:,4],1e-16),label='Commutator');ax[1].semilogy(N,np.maximum(d[:,5],1e-16),label='Full-state purity');ax[1].axhline(1e-6,color='gray',ls=':');ax[1].set(xlabel='N',ylabel='Numerical residual');ax[1].legend(fontsize=8);figure(p,'entropy_and_purity')
    fig,ax=plt.subplots(1,2,figsize=(7.4,3.3));selected=rows[np.isclose(rows[:,2],.01)]
    for M in [1.6,1.8,2.,2.4]:
        a=selected[selected[:,0]==M];valid=a[:,7]==1
        ax[0].plot(a[valid,1],a[valid,3],'-o',label=f'M={M}')
        ax[1].plot(a[valid,1],a[valid,4],'-o')
    ax[0].legend(fontsize=8);ax[0].set(xlabel='g',ylabel='Entropy at x=0.01');ax[1].set(xlabel='g',ylabel='Fractional power change');figure(p,'mass_mixing_scan')
    neg=diagnostics(modes(g=-.1,X=400,xf=.01,rtol=2.5e-14,atol=1e-17,state='adiabatic'),[.01])[0]
    pos=diagnostics(mm,[.01])[0]
    save_json(p/'results/summary.json',{'settings':SETTINGS,'last_reliable_x':float(xs[good][-1]),'max_zero_mixing_entropy':float(np.nanmax(zero[:,3])),'g_sign_entropy_difference':float(abs(neg[3]-pos[3])),'accepted_scan_rows':int(sum(rows[:,7])),'total_scan_rows':len(rows),'benchmark_entropy_x001':float(pos[3]),'benchmark_power_change_x001':float(pos[1]/(.5*(1+.01**2))-1),'software':software()})
if __name__=='__main__':main()
