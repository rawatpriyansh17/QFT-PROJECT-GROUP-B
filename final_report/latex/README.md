# Cumulative Research Report

Compile from this folder using pdfLaTeX (two passes):

```bash
mkdir -p build
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
```

The logo and figure inputs are included.

The report has one section for each completed Step (1–10), ending with a short line naming all three members and their specific contributions. These lines are generated in `contributions.tex` from the individual Step reports. After changing a Step contribution, run `python3 sync_contributions.py --build` from the repository root to synchronise the lines and rebuild `final_report/Final_Report.pdf`. Keep `main.tex`, `contributions.tex` and the included assets together when compiling.
