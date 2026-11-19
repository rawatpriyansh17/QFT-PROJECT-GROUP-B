# QFT Project B

**Curvature-spectator entanglement and the cosmological collider**

Priyansh Rawat, Samarth Bhatnagar and Mayank Sinha  
BITS Pilani, K. K. Birla Goa Campus | Supervisor: Professor Rudranil Basu

**[Combined report PDF](final_report/Final_Report.pdf)** · **[Combined report LaTeX](final_report/latex/)**

| Step | Calculation |
| --- | --- |
| [1](step1/) | Spacetime, fields and single-field reference |
| [2](step2/) | Two-oscillator Gaussian entropy |
| [3](step3/) | Canonical equations |
| [4](step4/) | Free modes and initial states |
| [5](step5/) | Gaussian entanglement |
| [6](step6/) | Heavy clock and in-in benchmark |
| [7](step7/) | Spectator bispectrum |
| [8](step8/) | Convergence and numerical errors |
| [9](step9/) | Mass recovery and fit stability |
| [10](step10/) | Entropy and the spectator signal |

Each directory contains the calculation, figures and numerical records for that step. The combined report ends with one contribution page covering the whole project.

## Reproduce

Use Python 3.12 from the repository root:

```bash
python3 -m pip install -r requirements.txt
python3 run_steps.py
```

To run selected steps, pass their numbers, for example `python3 run_steps.py 2 3`. Steps 8, 9 and 10 have data dependencies and should be run in that order. Step 1 has its own command in its README. To rebuild the PDF, follow the [report build instructions](final_report/latex/README.md).

The final result includes a validated entropy map, a resolved baseline heavy-clock shape and explicit fit-window/background sensitivity. The physical mass range is a numerical/model sensitivity envelope, not an observational confidence interval.
