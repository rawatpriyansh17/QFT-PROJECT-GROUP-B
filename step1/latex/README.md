# Step 1 report source

Main file: [step_01_report.tex](step_01_report.tex). The logo, vector figures and table inputs are included.

Compile with pdfLaTeX from the repository root:

```bash
cd step1/latex
mkdir -p build
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build step_01_report.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build step_01_report.tex
```

The report is written to `build/step_01_report.pdf`. Two passes resolve citations and cross-references.
