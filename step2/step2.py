"""Ground-state entanglement of two coupled oscillators."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *

SETTINGS={'K_diagonal': [1, 2], 'couplings': [0, 0.3, 0.6],
          'rescalings': [0.1, 0.5, 1, 2, 10],
          'curve_range': [0, 1.4], 'curve_points': 701}


def exact_entropy(c):
    """Closed-form S(c), valid for |c| < sqrt(2)."""
    c=np.asarray(c,dtype=float)
    if np.any(np.abs(c)>=np.sqrt(2)):
        raise ValueError('A stable ground state requires |c| < sqrt(2)')
    d=np.sqrt(2-c*c)
    nu=np.sqrt((2+d)*(1+d)/(4*d*(3+2*d)))
    return nu,entropy(nu)


def wavefunction_entropy(c):
    """Trace the ground-state Gaussian; use its geometric reduced spectrum."""
    w,U=np.linalg.eigh([[1.,c],[c,2.]])
    if np.min(w)<=0:
        raise ValueError('K must be positive')
    A=(U*np.sqrt(w))@U.T
    purity=float(np.sqrt(np.linalg.det(A)/(A[0,0]*A[1,1])))
    # Eigenvalues of the reduced Gaussian kernel: p_n=(1-t)*t**n.
    t=max(0.,(1-purity)/(1+purity))
    S=-np.log1p(-t)-(t/(1-t))*np.log(t) if t>0 else 0.
    return purity,float(S)

def main():
    p=output(2);rows=[];max_rescale=0.
    for c in [0,.3,.6]:
        X,P,nu,S=oscillator(c)
        for b in [.1,.5,1.,2.,10.]:
            det=(b*b*X[0,0])*(P[0,0]/b**2)
            max_rescale=max(max_rescale,abs(det-nu*nu))
            rows.append([c,b,nu,S,2*S,det])
    assert abs(oscillator(.3)[2]-.502816542)<1e-9
    assert abs(oscillator(.3)[3]-.0193599304)<1e-9
    assert oscillator(0)[3]<1e-12 and max_rescale<1e-13
    save_table(p/'results/oscillators.csv','coupling,rescaling,nu,entropy_nats,pair_entropy_nats,determinant',rows)
    cs=np.linspace(*SETTINGS['curve_range'],SETTINGS['curve_points'])
    nus,ss=exact_entropy(cs)
    covariance=np.array([oscillator(c)[3] for c in cs])
    gaussian=np.array([wavefunction_entropy(c) for c in cs])
    residual=gaussian[:,1]-ss
    max_wave_error=float(np.max(abs(residual)))
    max_cov_error=float(np.max(abs(covariance-ss)))
    max_purity_error=float(np.max(abs(gaussian[:,0]-1/(2*nus))))
    assert max_wave_error<1e-12 and max_cov_error<1e-12 and max_purity_error<1e-12
    assert np.max(abs(exact_entropy(-cs)[1]-ss))<1e-14
    save_table(p/'results/entropy_curve.csv',
               'coupling,exact_nu,exact_entropy_nats,covariance_entropy_nats,gaussian_purity,gaussian_entropy_nats,gaussian_minus_exact_nats',
               np.column_stack([cs,nus,ss,covariance,gaussian,residual]))
    plt.figure(figsize=(6.4,3.5))
    plt.plot(cs,ss,color='#2166ac',lw=2.8,label=r'Exact $S(c)=s(\nu(c))$')
    plt.plot(cs,gaussian[:,1],color='#d95f02',lw=1.7,ls='--',
             label='Gaussian wavefunction (check 2)')
    plt.scatter(SETTINGS['couplings'],[oscillator(c)[3] for c in SETTINGS['couplings']],
                color='black',s=30,zorder=4,label='Three covariance benchmarks')
    plt.xlim(-.025,1.425);plt.ylim(bottom=-.02)
    plt.xlabel('Off-diagonal coupling c');plt.ylabel('Entropy of oscillator 1 (nats)')
    plt.legend(loc='upper left',frameon=False,fontsize=9)
    plt.grid(alpha=.18);plt.tight_layout();figure(p,'oscillator_entropy')
    # Keep the vector input used by the report in sync with the plotted PNG.
    import shutil
    shutil.copyfile(p/'figures/oscillator_entropy.pdf',p/'latex/figures/oscillator_entropy.pdf')
    save_json(p/'results/summary.json',{'settings':SETTINGS,'benchmark_nu':oscillator(.3)[2],'benchmark_entropy':oscillator(.3)[3],'coupling_0_6_entropy':oscillator(.6)[3],'max_rescaling_determinant_error':max_rescale,'max_wavefunction_entropy_error_nats':max_wave_error,'max_covariance_entropy_error_nats':max_cov_error,'max_gaussian_purity_error':max_purity_error,'checks_passed':True,'software':software()})
if __name__=='__main__':main()
