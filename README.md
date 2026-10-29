# QFT Project B

**Curvature-spectator entanglement and the cosmological collider**

Priyansh Rawat, Samarth Bhatnagar and Mayank Sinha  
BITS Pilani, K. K. Birla Goa Campus | Supervisor: Professor Rudranil Basu

**[Final report PDF](final_report/Final_Report.pdf)** · **[Final report LaTeX](final_report/latex/)**

| Step | Calculation | Report |
| --- | --- | --- |
| [1](step1/) | Spacetime, fields and single-field reference | [PDF](step1/Step_01_Report.pdf) |
| [2](step2/) | Two-oscillator Gaussian entropy | [PDF](step2/Step_02_Report.pdf) |
| [3](step3/) | Canonical equations | [PDF](step3/Step_03_Report.pdf) |
| [4](step4/) | Free modes and initial states | [PDF](step4/Step_04_Report.pdf) |
| [5](step5/) | Gaussian entanglement | [PDF](step5/Step_05_Report.pdf) |
| [6](step6/) | Heavy clock and in-in benchmark | [PDF](step6/Step_06_Report.pdf) |
| [7](step7/) | Spectator bispectrum | [PDF](step7/Step_07_Report.pdf) |
| [8](step8/) | Convergence and numerical errors | [PDF](step8/Step_08_Report.pdf) |
| [9](step9/) | Mass recovery and fit stability | [PDF](step9/Step_09_Report.pdf) |
| [10](step10/) | Entropy and the spectator signal | [PDF](step10/Step_10_Report.pdf) |

Each step explains its essential figures and data, with the calculation, LaTeX source and downloadable PDF.

## Reproduce

Use Python 3.12 from the repository root:

```bash
python3 -m pip install -r requirements.txt
python3 run_steps.py
```

To run selected steps, pass their numbers, for example `python3 run_steps.py 2 3`. Steps 8, 9 and 10 have data dependencies and should be run in that order. Step 1 has its own command in its README. Report sources compile with pdfLaTeX from each `latex/` folder.

The final result includes a validated entropy map, a resolved baseline heavy-clock shape and explicit fit-window/background sensitivity. The physical mass range is a numerical/model sensitivity envelope, not an observational confidence interval.
