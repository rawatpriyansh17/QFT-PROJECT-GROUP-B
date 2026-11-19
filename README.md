# What Quantum Trace Can a Heavy Particle Leave?

**[Read the final report (PDF)](Final_Report.pdf)**

Curvature–spectator entanglement and the cosmological collider. A project by Priyansh Rawat, Mayank Sinha and Samarth Bhatnagar, supervised by Professor Rudranil Basu at BITS Pilani, K. K. Birla Goa Campus.

| File | Purpose |
| --- | --- |
| [Final_Report.tex](Final_Report.tex) | Single LaTeX entry point for the report. |
| [assets/](assets/) | The institute logo and five plots |
| [calculations.py](calculations.py) | Numerical work behind the power, entropy, bispectrum and mass-recovery results. |

## Compile the report

From the repository root, run pdfLaTeX twice to resolve the contents and cross-references.

```sh
pdflatex Final_Report.tex
pdflatex Final_Report.tex
```

## Reproduce a result

With Python 3.12, the following commands compute the regulated-integral error budget, fit the clock on 50 triangles and evaluate the cubic-coupling controls.

```sh
python3 -m pip install numpy scipy mpmath
python3 calculations.py clock
```

The full-window fit gives $m/H=1.9930$ for an input of $2$. Across the twelve fit choices the range is $1.9930$–$2.0198$. Numerical residuals close to machine precision can vary with the solver and platform.

Use `reference` for the squeezed-limit and oscillator benchmarks, `gaussian` for the mode and covariance checks, `benchmark` for the analytic in-in integral, `pilot` for the regulated bispectrum, or `mock` for the seeded recovery trials. The `grid` calculation evaluates all 625 mass–mixing pairs in Figure 10. Running `all` performs every calculation. Generated tables go to a local `results/` directory, or to the directory passed with `--out`.
