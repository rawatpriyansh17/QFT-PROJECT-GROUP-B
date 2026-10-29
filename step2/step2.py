"""Ground-state entanglement of two coupled oscillators."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common.physics import *

SETTINGS={'K_diagonal': [1, 2], 'couplings': [0, 0.3, 0.6], 'rescalings': [0.1, 0.5, 1, 2, 10]}

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
    # Independent Gaussian-kernel result: Tr rho_1^2=1/(2 nu).
    for c in [.3,.6]:
        X,P,nu,S=oscillator(c);A=np.linalg.inv(2*X)
        kernel_purity=np.sqrt(np.linalg.det(A)/(A[0,0]*A[1,1]))
        assert abs(kernel_purity-1/(2*nu))<1e-13
    save_table(p/'results/oscillators.csv','coupling,rescaling,nu,entropy_nats,pair_entropy_nats,determinant',rows)
    cs=np.linspace(0,1.2,151);ss=[oscillator(c)[3] for c in cs]
    plt.figure(figsize=(6.4,3.4));plt.plot(cs,ss);plt.scatter([0,.3,.6],[oscillator(c)[3] for c in [0,.3,.6]])
    plt.xlabel('Off-diagonal coupling c');plt.ylabel('Entropy of oscillator 1 (nats)');figure(p,'oscillator_entropy')
    save_json(p/'results/summary.json',{'settings':SETTINGS,'benchmark_nu':oscillator(.3)[2],'benchmark_entropy':oscillator(.3)[3],'coupling_0_6_entropy':oscillator(.6)[3],'max_rescaling_determinant_error':max_rescale,'checks_passed':True,'software':software()})
if __name__=='__main__':main()
