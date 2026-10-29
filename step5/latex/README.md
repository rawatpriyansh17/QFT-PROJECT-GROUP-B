# Step 5: Gaussian Entanglement

Compile from this folder using pdfLaTeX (two passes):

```bash
mkdir -p build
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=build main.tex
```

The logo and figure inputs are included.
