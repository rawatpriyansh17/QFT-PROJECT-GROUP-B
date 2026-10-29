# Final Research Report

Compile from this folder using pdfLaTeX (two passes):

```bash
mkdir -p build
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
```

The logo and figure inputs are included.

The contribution page is generated from the ten individual step reports. After changing a step contribution, run `python3 sync_contributions.py --build` from the repository root to synchronise it and rebuild the final PDF.
