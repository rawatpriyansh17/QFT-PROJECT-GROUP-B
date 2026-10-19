#!/usr/bin/env python3
"""Reproduce Project B, Step 1: horizon conventions and the squeezed reference.

No coupled-field evolution or bispectrum integration is performed here.
The reference spectrum has arbitrary positive amplitude; it is not a data fit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STEP = Path(__file__).resolve().parent


def power(k, ns, amplitude=1.0, pivot=1.0):
    """P_zeta(k) = A/k^3 * (k/k_star)^(n_s-1), for k>0."""
    k = np.asarray(k, dtype=float)
    if np.any(k <= 0) or amplitude <= 0 or pivot <= 0:
        raise ValueError("k, amplitude and pivot must be positive")
    return amplitude * (k / pivot) ** (ns - 1.0) / k**3


def k_derivative_response(k, ns, h, amplitude=1.0, pivot=1.0):
    """Approximate -(3+k*d/dk)P/P using a fractional central difference."""
    pk = power(k, ns, amplitude, pivot)
    kp = power(k * (1 + h), ns, amplitude, pivot)
    km = power(k * (1 - h), ns, amplitude, pivot)
    return -3.0 - (kp - km) / (2 * h * pk)


def dilation_response(k, ns, eps, amplitude=1.0, pivot=1.0):
    """Differentiate e^(-3*zeta_L) P(k*e^(-zeta_L)) at zeta_L=0."""
    plus = np.exp(-3 * eps) * power(k * np.exp(-eps), ns, amplitude, pivot)
    minus = np.exp(3 * eps) * power(k * np.exp(eps), ns, amplitude, pivot)
    return (plus - minus) / (2 * eps * power(k, ns, amplitude, pivot))


def save_csv(path, names, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(names)
        writer.writerows(rows)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=STEP / "config.json")
    args = parser.parse_args()
    cfg_path = args.config.resolve()
    cfg = json.loads(cfg_path.read_text())
    if not (cfg["H"] > 0 and cfg["k_horizon"] > 0 and
            0 < cfg["k_min"] < cfg["k_max"] and
            0 < cfg["x_min"] < 1 < cfg["x_max"] and
            cfg["k_points"] >= 3 and cfg["x_points"] >= 3 and
            all(0 < h < 1 for h in cfg["fractional_steps"]) and
            all(e > 0 for e in cfg["dilation_steps"])):
        raise ValueError("Invalid positive scales, grids, or finite-difference steps")
    results = STEP / "results"
    figures = STEP / "figures"
    results.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    k = np.geomspace(cfg["k_min"], cfg["k_max"], cfg["k_points"])
    ns, amp, pivot = cfg["ns"], cfg["amplitude"], cfg["pivot"]
    expected = 1 - ns
    convergence = []
    for h in cfg["fractional_steps"]:
        response = k_derivative_response(k, ns, h, amp, pivot)
        convergence.append((h, float(np.max(np.abs(response - expected)))))
    dilations = []
    for eps in cfg["dilation_steps"]:
        response = dilation_response(k, ns, eps, amp, pivot)
        dilations.append((eps, float(np.max(np.abs(response - expected)))))
    h = min(cfg["fractional_steps"])
    eps = min(cfg["dilation_steps"])
    derivative = k_derivative_response(k, ns, h, amp, pivot)
    dilation = dilation_response(k, ns, eps, amp, pivot)
    controls = []
    for tilt in cfg["ns_controls"]:
        target = 1 - tilt
        numerical = dilation_response(k, tilt, eps, amp, pivot)
        controls.append((tilt, target, float(np.max(np.abs(numerical - target)))))

    # These are different checks: direct k differentiation, finite coordinate
    # dilation, a scale-invariant control, and normalization invariance.
    tol = cfg["validation_abs_tolerance"]
    checks = {
        "requested_reference_is_0_035": abs(expected - 0.035) < 1e-14 if ns == 0.965 else True,
        "k_derivative_matches_reference": bool(np.max(np.abs(derivative - expected)) < tol),
        "finite_dilation_matches_reference": bool(np.max(np.abs(dilation - expected)) < tol),
        "two_numerical_routes_agree": bool(np.max(np.abs(derivative - dilation)) < tol),
        "scale_invariant_control_is_zero": bool(np.max(np.abs(dilation_response(k, 1.0, eps))) < tol),
        "red_and_blue_tilt_controls_match": all(row[2] < tol for row in controls),
        "amplitude_and_pivot_cancel": bool(np.max(np.abs(
            dilation_response(k, ns, eps, 7 * amp, 2 * pivot) - dilation)) < tol),
        "central_difference_second_order": bool(3.8 < convergence[0][1] / convergence[1][1] < 4.2),
    }
    H, kh = cfg["H"], cfg["k_horizon"]
    x = np.geomspace(cfg["x_max"], cfg["x_min"], cfg["x_points"])
    eta = -x / kh
    a = -1 / (H * eta)
    ahorizon = kh / H
    checks["horizon_identity"] = bool(np.allclose(x, kh / (a * H), rtol=1e-14))
    checks["horizon_crossing"] = bool(np.isclose(kh / (ahorizon * H), 1))
    save_csv(results / "squeezed_reference.csv",
             ["k_over_pivot", "P_zeta", "analytic_ratio", "k_derivative_ratio", "dilation_ratio"],
             zip(k / pivot, power(k, ns, amp, pivot), np.full_like(k, expected), derivative, dilation))
    save_csv(results / "derivative_convergence.csv", ["fractional_step", "max_absolute_error"], convergence)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.bbox": "tight"})
    fig, ax = plt.subplots(figsize=(7.1, 2.75), constrained_layout=True)
    ax.axvspan(1, cfg["x_max"], color="#edf3f9")
    ax.axvspan(cfg["x_min"], 1, color="#f6f0e4")
    ax.loglog(x, a / ahorizon, color="#193b60", lw=2.2)
    ax.scatter([1], [1], color="#b16f13", zorder=4)
    ax.axvline(1, color="0.5", ls="--", lw=0.8)
    ax.set(xlim=(cfg["x_max"], cfg["x_min"]), xlabel=r"$x=k/(aH)$  (time evolves to the right)",
           ylabel=r"$a/a_*$", title="de Sitter horizon convention")
    ax.text(20, 14, "Inside: x > 1", ha="center")
    ax.text(0.05, 0.05, "Outside: x < 1", ha="center")
    ax.annotate("Crossing: x = 1", xy=(1, 1), xytext=(4, 0.17),
                arrowprops={"arrowstyle": "->", "color": "0.3"})
    fig.savefig(figures / "horizon_map.pdf")
    fig.savefig(figures / "horizon_map.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.95), constrained_layout=True)
    axes[0].axhline(expected, color="#b16f13", lw=2, label="Analytic: 1 - ns")
    axes[0].semilogx(k / pivot, derivative, color="#193b60", lw=1.2, ls="--", label="k derivative")
    axes[0].semilogx(k[::8] / pivot, dilation[::8], "o", color="#38826d", ms=3, label="Finite dilation")
    axes[0].set(xlabel=r"$k/k_*$", ylabel=r"$B_\zeta/(P_L P_S)$ reference",
                ylim=(expected - 0.001, expected + 0.001), title="Reference agreement")
    axes[0].legend(fontsize=8, loc="upper right")
    hs, errors = np.array(convergence).T
    axes[1].loglog(hs, errors, "o-", color="#193b60", label="k-derivative error")
    axes[1].loglog(hs, errors[0] * (hs / hs[0])**2, "--", color="#b16f13", label=r"$h^2$ guide")
    axes[1].set(xlabel="Fractional step h", ylabel="Maximum absolute error", title="Numerical refinement")
    axes[1].legend(fontsize=8)
    fig.savefig(figures / "reference_validation.pdf")
    fig.savefig(figures / "reference_validation.png", dpi=180)
    plt.close(fig)

    environment = {"python": platform.python_version(), "numpy": np.__version__,
                   "matplotlib": matplotlib.__version__,
                   "reproduce_command": "python3 step1/step1.py --config step1/config.json"}
    summary = {
        "scope": "Step 1 conventions and single-clock reference only",
        "ns": ns, "reference_ratio": expected, "k_points": len(k),
        "k_range": [float(k[0]), float(k[-1])],
        "finest_fractional_step": h, "finest_dilation_step": eps,
        "max_derivative_error": convergence[-1][1],
        "max_dilation_error": dilations[-1][1],
        "tolerance_absolute": tol, "checks": checks,
        "all_checks_passed": all(checks.values()),
        "random_seed": None,
        "randomness": "None; all computations are deterministic",
        "source_sha256": sha256(Path(__file__)),
        "config_sha256": sha256(cfg_path),
        "software": environment,
    }
    (results / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not summary["all_checks_passed"]:
        raise AssertionError(f"Validation failed: {checks}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
