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

The checked outputs are included in [`results/`](results/) and [`figures/`](figures/). Running the script regenerates them. All 10 checks must pass; otherwise the script raises an error.

The calculation is deterministic. No coupled-field entropy or spectator bispectrum is computed in this step. The tilted single-field reference is separate from the exact de Sitter background of the two-field model and is not added to a spectator bispectrum.

## Saved results

- [Validation summary](results/summary.json): all 10 checks passed, with maximum numerical errors and source/configuration hashes.
- [Software versions and command](results/environment.json).
- [Squeezed reference table](results/squeezed_reference.csv) and [horizon map data](results/horizon_map.csv).
- [Direct-derivative convergence](results/derivative_convergence.csv), [dilation convergence](results/dilation_convergence.csv), and [tilt controls](results/tilt_controls.csv).

## Figures

![de Sitter horizon map](figures/horizon_map.png)

![Single-field reference validation](figures/reference_validation.png)

Vector versions: [horizon map PDF](figures/horizon_map.pdf) and [reference validation PDF](figures/reference_validation.pdf).
