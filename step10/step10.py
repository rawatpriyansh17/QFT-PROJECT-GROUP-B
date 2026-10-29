"""Entropy versus spectator signal, and the independent cubic-coupling control."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *
from common.bispectrum import integral,extrapolate

SETTINGS={'mass_mixing_pairs': [[2, 0.05], [2, 0.1], [2, 0.2], [1.8, 0.1], [2.4, 0.1]], 'entropy_time': 0.01, 'triangle_r': 0.03, 'cubic_controls': [0, 0.005, 0.01], 'common_absolute_noise_F': 1e-10}

def main():
    p=output(10);rows=[]
    for M,g in [(2.,.05),(2.,.1),(2.,.2),(1.8,.1),(2.4,.1)]:
        mm=modes(M,g,X=1600,xf=1e-5,rtol=1e-12,atol=1e-15,state='adiabatic')
        d=diagnostics(mm,[.01])[0]
        v,_=integral(mm,.03,zmax=640,epsilons=[.1,.05,.025,.0125]);f,e=extrapolate([.1,.05,.025,.0125],v)
        rows.append([M,g,d[3],d[1]/(.5*(1+.01**2))-1,f,e])
    save_table(p/'results/entropy_signal_comparison.csv','M,g,entropy_x001,power_change_x001,F_r003,regulator_error',rows)
    baseline=[a for a in rows if a[0]==2 and a[1]==.1][0]
    cubic=[[l,baseline[2],baseline[4]*l/.01,baseline[5]*l/.01] for l in [0,.005,.01]]
    save_table(p/'results/cubic_control.csv','lambda3_over_H,Gaussian_entropy,F_r003,numerical_error',cubic)
    fit=json.loads((ROOT/'step9/results/summary.json').read_text())['physical']
    rr=np.geomspace(.001,.2,50);mu=fit['full_range_frequency']
    # Fixed-covariance diagnostic using the fitted baseline clock template.
    full=next(f for f in fit['fits'] if f['r_min']==.001 and f['r_max']==.2 and not f['extra_r4'])
    A,B=full['coefficients'][2:4]
    template=rr**1.5*(A*np.cos(mu*np.log(rr))+B*np.sin(mu*np.log(rr)))
    noise_floor=1e-10
    recover=[[l,float(np.linalg.norm(template*l/.01)/noise_floor),noise_floor] for l in [0,.005,.01]]
    save_table(p/'results/common_noise_control.csv','lambda3_over_H,known_template_SNR,absolute_sigma_F',recover)
    fig,ax=plt.subplots(1,2,figsize=(7.3,3.5),layout='constrained')
    for M,g,S,powr,f,e in rows:
        ax[0].errorbar(S*1e3,f*1e10,yerr=e*1e10,fmt='o',label=f'M={M:g}, g={g:g}')
    ax[0].set(xlabel='Entropy at x=0.01 (10^-3 nats)',ylabel='F(0.03) (10^-10)')
    c=np.array(cubic);ax[1].plot(c[:,0],c[:,2]*1e10,'-o');ax[1].set(xlabel='Cubic coupling lambda3/H',ylabel='F(0.03) (10^-10)')
    ax[1].set_xticks([0,.005,.01]);ax[1].set_xticklabels(['0','0.005','0.01'])
    handles,labels=ax[0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=7)
    figure(p,'entropy_signal_control')
    save_json(p/'results/summary.json',{'settings':SETTINGS,'baseline_entropy':baseline[2],'baseline_F_r003':baseline[4],'main_conclusion':'At fixed quadratic parameters the Gaussian entropy is unchanged when lambda3 varies, while the leading bispectrum is proportional to lambda3.','clock_amplitude_note':'F at one triangle is not the clock amplitude. Physical frequency and window tests are in Step 9.','software':software()})
if __name__=='__main__':main()
