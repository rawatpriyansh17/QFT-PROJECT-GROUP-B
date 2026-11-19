# Build the combined report

From the repository root, with a standard pdfLaTeX installation and Python with PyMuPDF installed, run

```bash
python3 -m pip install PyMuPDF
python3 final_report/build_report.py
```

The script compiles twice, restores the original foundation-page hyperlinks and writes `final_report/Final_Report.pdf`.

- `main.tex` supplies the cover, contents, numbering and document structure.
- `assets/foundation.pdf` contains the typeset scientific pages from Samarth’s Steps 1–6, including his analytic benchmark. The original per-step credit tables are consolidated on the closing page.
- `bispectrum.tex`, `convergence.tex`, `mass.tex` and `entropy.tex` contain the continuation. Their figures are in `figures/`.
- `contributions.tex` contains the project-wide attribution table.
