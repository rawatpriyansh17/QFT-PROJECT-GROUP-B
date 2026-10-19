# Step 1: Spacetime, fields and observables

The script checks the single-field reference `B_zeta/(P_L P_S) -> 1 - ns = 0.035` at `ns = 0.965`, using direct differentiation and a finite spatial dilation. It also checks the de Sitter horizon convention `x = k/(aH)`.

## Run

From the repository root, using Python 3.12:

```bash
python3 -m pip install -r step1/requirements.txt
python3 step1/step1.py --config step1/config.json
```

Or, from this folder:

```bash
python3 -m pip install -r requirements.txt
python3 step1.py --config config.json
```

## Files

| File | Purpose |
| --- | --- |
| `step1.py` | Numerical checks and figure generation |
| `config.json` | Physical inputs, grids and tolerances |
| `requirements.txt` | Package versions used for the calculation |

Running the script creates `results/` with CSV tables, check results, software versions and source/configuration hashes, and `figures/` with the horizon map and validation plots. All 10 checks must pass; otherwise the script raises an error.

The calculation is deterministic. No coupled-field entropy or spectator bispectrum is computed in this step. The tilted single-field reference is separate from the exact de Sitter background of the two-field model and is not added to a spectator bispectrum.
