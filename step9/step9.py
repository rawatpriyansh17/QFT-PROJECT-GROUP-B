"""Mock-data calibration, no-clock controls and convergence-qualified mass fits."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *
from common.bispectrum import integral,extrapolate
from common.fitting import profile,design

SETTINGS={'M': 2, 'g': 0.1, 'r_min': 0.001, 'r_max': 0.2, 'r_count': 50, 'mock_A': 0.1, 'mock_B': 0.05, 'mock_b0': 0, 'mock_b2': 0.2, 'mock_seed': 20261029, 'held_out_null_seeds': [20261031, 20261035], 'trials_each': 100, 'profile_mu_grid': [0.1, 4, 0.01], 'fit_upper_endpoints': [0.05, 0.1, 0.2], 'fit_lower_endpoints': [0.001, 0.003], 'background_controls': ['1,r2', '1,r2,r4']}

def main():
    p=output(9);r=np.geomspace(1e-3,.2,50);mu=np.sqrt(1.76)
    coef=np.array([0.,.2,.1,.05]);signal=design(r,mu)@coef
    oscill=design(r,mu)[:,2:]@coef[2:];rms=float(np.sqrt(np.mean(oscill**2)))
    grid=np.arange(.1,4.0001,.01)
    Q=np.array([np.linalg.qr(design(r,m))[0] for m in grid])
    rng=np.random.default_rng(20261029);rows=[];nulls={};coverage={}
    for frac in [.01,.05]:
        sigma=frac*rms;cov=0;est=[];null=[]
        for j in range(100):
            y=signal+rng.normal(0,sigma,len(r))
            fit=profile(r,y,sigma)
            inside=any(a<=mu<=b for a,b in fit['intervals_95']);cov+=inside;est.append(fit['mu'])
            rows.append([frac,j,fit['mu'],np.sqrt(fit['mu']**2+2.25-.01),int(inside)])
            yn=design(r,mu)[:,:2]@coef[:2]+rng.normal(0,sigma,len(r))
            nf=profile(r,yn,sigma)
            A=design(r,mu)[:,:2]/sigma;bb=yn/sigma
            chi0=np.sum((A@np.linalg.lstsq(A,bb,rcond=None)[0]-bb)**2)
            null.append(float(chi0-nf['chi2']))
            if j==0:
                save_table(p/f'results/mock_profile_{int(frac*100)}pct.csv','mu,delta_chi2',np.c_[fit['grid'],fit['delta_chi2']])
        threshold=float(np.quantile(null,.95))
        held=np.random.default_rng(20261030+int(frac*100));false_count=0
        for j in range(100):
            yn=design(r,mu)[:,:2]@coef[:2]+held.normal(0,sigma,len(r))
            nf=profile(r,yn,sigma);A=design(r,mu)[:,:2]/sigma;bb=yn/sigma
            chi0=np.sum((A@np.linalg.lstsq(A,bb,rcond=None)[0]-bb)**2)
            false_count+=chi0-nf['chi2']>threshold
        nulls[str(frac)]={'threshold_95':threshold,'calibration_false_alarm_rate':float(np.mean(np.array(null)>threshold)),'held_out_false_alarm_rate':false_count/100,'held_out_trials':100}
        coverage[str(frac)]={'coverage':cov/100,'frequency_bias':float(np.mean(est)-mu),'frequency_RMSE':float(np.sqrt(np.mean((np.array(est)-mu)**2)))}
    save_table(p/'results/mock_recovery.csv','noise_fraction,trial,mu_estimate,mass_estimate,interval_contains_truth',rows)
    fig,ax=plt.subplots(1,2,figsize=(7.4,3.2))
    for frac in [.01,.05]:
        a=np.array(rows);v=a[a[:,0]==frac,2];ax[0].hist(v,bins=14,alpha=.5,label=f'{frac*100:g}% noise')
        a=np.loadtxt(p/f'results/mock_profile_{int(frac*100)}pct.csv',delimiter=',',skiprows=1);ax[1].plot(a[:,0],a[:,1],label=f'{frac*100:g}%')
    ax[0].axvline(mu,color='black',ls=':');ax[0].set(xlabel='Recovered frequency',ylabel='Mock realisations');ax[0].legend(fontsize=8)
    ax[1].axhline(3.841459,color='gray',ls=':');ax[1].set(xlabel='Frequency mu',ylabel='Profile delta chi squared',ylim=(0,15));ax[1].legend();figure(p,'mock_mass_recovery')
    gate=json.loads((ROOT/'step8/results/summary.json').read_text())
    physical={'status':'Not fitted: pilot convergence gate failed'}
    if gate['pilot_below_one_percent_of_clock_envelope']:
        model=modes(X=1600,xf=1e-6,rtol=1e-12,atol=1e-15,state='adiabatic')
        vals=[];errors=[]
        budget=np.loadtxt(ROOT/'step8/results/error_budget.csv',delimiter=',',skiprows=1)
        budget=budget[np.argsort(budget[:,0])]
        other_error=budget[:,-1]-budget[:,3]
        eps=[.1,.05,.025,.0125]
        for rr in r:
            vv,_=integral(model,float(rr),zf=.01,zmax=640,epsilons=eps)
            f,e=extrapolate(eps,vv);vals.append(f);extra=float(np.interp(np.log(rr),np.log(budget[:,0]),other_error));errors.append(max(e+extra,1e-14))
        vals=np.array(vals);errors=np.array(errors)
        save_table(p/'results/physical_bispectrum.csv','r,F,conservative_numerical_error',np.c_[r,vals,errors])
        fits=[]
        for upper in [.05,.1,.2]:
            for lower in [.001,.003]:
                for extra in [False,True]:
                    use=(r>=lower)&(r<=upper)
                    fit=profile(r[use],vals[use],errors[use],extra)
                    fits.append({'r_min':lower,'r_max':upper,'extra_r4':extra,**{k:v for k,v in fit.items() if k not in ['grid','delta_chi2']}})
        full=profile(r,vals,errors);amplitude=float(np.hypot(*full['coefficients'][2:4]))
        envelope_error=float(max(errors/r**1.5))
        physical={'status':'Leading bispectrum fitted; numerical profile widths are diagnostics, not statistical confidence intervals','fits':fits,'full_range_frequency':full['mu'],'clock_amplitude':amplitude,'numerical_error_over_clock_envelope':envelope_error/max(amplitude,1e-30),'visible_cycles':float(full['mu']*np.log(r[-1]/r[0])/(2*np.pi)),'frequency_best_fit_envelope':[min(f['mu'] for f in fits),max(f['mu'] for f in fits)],'mass_best_fit_envelope':[float(np.sqrt(min(f['mu'] for f in fits)**2+2.24)),float(np.sqrt(max(f['mu'] for f in fits)**2+2.24))]}
        fig,ax=plt.subplots(1,2,figsize=(7.5,3.3));ax[0].errorbar(r,vals/r**1.5,yerr=errors/r**1.5,fmt='.',ms=3);ax[0].plot(r,(design(r,full['mu'])@np.array(full['coefficients']))/r**1.5);ax[0].set_xscale('log');ax[0].set(xlabel='r',ylabel='F / r^(3/2)')
        ax[1].plot(full['grid'],full['delta_chi2']);ax[1].axvline(np.sqrt(1.76),color='gray',ls=':',label='Expected frequency');ax[1].set(xlabel='Frequency mu (minimum region)',ylabel='Numerically weighted profile',xlim=(full['mu']-.018,full['mu']+.018),ylim=(0,25));ax[1].legend(fontsize=7);figure(p,'physical_clock_fit')
    save_json(p/'results/summary.json',{'settings':SETTINGS,'seed':20261029,'trials_per_noise':100,'true_frequency':mu,'mock_rms':rms,'mock_results':coverage,'no_clock_controls':nulls,'physical':physical,'software':software()})
if __name__=='__main__':main()
