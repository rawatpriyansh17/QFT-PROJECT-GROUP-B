"""Heavy-clock exponents and an independently integrated analytic in-in test."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *
from scipy.integrate import quad

SETTINGS={'M': 2, 'g': 0.1, 'gamma': 1, 'H': 1, 'r': [1, 0.1, 0.03, 0.01], 'epsilons': [0.1, 0.05, 0.025, 0.01, 0.001]}

def main():
    p=output(6);rows=[]
    for r in [1.,.1,.03,.01]:
        K=2+r;exact=-1/(2*r*K**3)
        for eps in [.1,.05,.025,.01,.001]:
            val,err=quad(lambda y:y*y*np.exp(-eps*y),0,np.inf,weight='sin',wvar=1.,epsabs=1e-10,limit=500)
            numerical=val/(4*r*K**3)
            reg=-np.imag(2/(eps+1j)**3)/(4*r*K**3)
            assert abs(numerical-reg)<1e-8*abs(exact)
            rows.append([r,eps,exact,numerical,abs(numerical/exact-1),err])
    rows=np.array(rows);assert np.max(rows[rows[:,1]==.001,4])<.001
    save_table(p/'results/in_in_benchmark.csv','r,epsilon,exact_bispectrum,numerical_bispectrum,relative_error,quadrature_error',rows)
    plt.figure(figsize=(6.5,3.4))
    for r in [1.,.1,.03,.01]:
        a=rows[rows[:,0]==r];plt.loglog(a[:,1],a[:,4],'-o',label=f'r={r:g}')
    plt.xlabel('Regulator epsilon');plt.ylabel('Relative error from zero-regulator answer');plt.legend();figure(p,'in_in_benchmark')
    mu=np.sqrt(2**2+.1**2-2.25)
    save_json(p/'results/summary.json',{'settings':SETTINGS,'clock_frequency_M2_g01':mu,'test_gamma':1.,'H':1.,'max_error_epsilon_0001':float(max(rows[rows[:,1]==.001,4])),'quadrature_matches_regulated_formula':True,'software':software()})
if __name__=='__main__':main()
