"""Numerical calculations for Curvature–Spectator Entanglement and the Cosmological Collider.

Run with --help for the individual calculations. CSV tables are written locally.
The report builds independently of this module.
"""

import argparse
from pathlib import Path
import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

OMEGA = np.kron(np.eye(2), np.array([[0.0, 1.0], [-1.0, 0.0]]))


def entropy(nu):
    nu = np.asarray(nu, dtype=float)
    if np.min(nu) < 0.5 - 1e-7:
        raise ArithmeticError("Unphysical reduced covariance")
    n = np.maximum(nu - 0.5, 0)
    return (n + 1) * np.log1p(n) - np.where(n > 0, n * np.log(np.maximum(n, 1e-300)), 0)


def oscillator(c):
    K = np.array([[1.0, c], [c, 2.0]])
    w, U = np.linalg.eigh(K)
    if np.min(w) <= 0:
        raise ValueError("K must be positive")
    X = (U * (0.5 / np.sqrt(w))) @ U.T
    P = (U * (0.5 * np.sqrt(w))) @ U.T
    nu = np.sqrt(X[0, 0] * P[0, 0])
    return X, P, float(nu), float(entropy(nu))


def bd(x, M):
    mp.mp.dps = 35
    z = mp.mpf(float(x))
    nu = mp.sqrt(mp.mpf("2.25") - M * M)
    pref = mp.sqrt(mp.pi) / 2 * mp.exp(mp.j * mp.pi * (nu / 2 + mp.mpf(".25")))
    h = mp.hankel1(nu, z)
    dh = (mp.hankel1(nu - 1, z) - mp.hankel1(nu + 1, z)) / 2
    sig = pref * z ** mp.mpf("1.5") * h
    ds = pref * (mp.mpf("1.5") * mp.sqrt(z) * h + z ** mp.mpf("1.5") * dh)
    q = (1 - 1j * x) * np.exp(1j * x) / np.sqrt(2)
    pq = -np.exp(1j * x) / (np.sqrt(2) * x)
    return q, pq, complex(sig), complex(-ds / z**2)


def initial(X, M, g, state="product", order=6):
    if state == "adiabatic":
        I = np.eye(2)
        J = np.array([[0.0, 1.0], [-1.0, 0.0]])
        D = np.array([[0.0, -3 * g], [0.0, -M * M]])
        coeff = {1: I + g * J / 2}
        for n in range(2, order + 1):
            cross = sum(
                (coeff[j] @ coeff[n - j] for j in range(1, n)),
                start=np.zeros((2, 2), complex),
            )
            coeff[n] = (
                ((n + 1) * I + g * J) @ coeff[n - 1] - cross + (D if n == 2 else 0)
            ) / (2j)
        R = 1j * I + sum(
            (coeff[n] / X**n for n in coeff), start=np.zeros((2, 2), complex)
        )
        A = -R / X**2 + g * np.array([[0.0, 1.0], [0.0, 0.0]]) / X**3
        assert np.max(abs(A - A.T)) < 1e-13
        C = 0.5 * X**2 * np.linalg.inv(R.imag)
        w, U = np.linalg.eigh(C)
        if np.min(w) <= 0:
            raise ValueError("Adiabatic covariance not positive")
        F = (U * np.sqrt(w)) @ U.T
        vel = -X * R @ F
        return np.array([F[0], vel[0], F[1], vel[1]])
    q, pq, s, ps = bd(X, M)
    u = np.zeros((4, 2), complex)
    u[:, 0] = [q, X**3 * pq, 0, 0]
    u[:, 1] = [0, -g * s, s, X**3 * ps]
    return u


def rhs_N(N, y, M, g):
    q, v, s, w = y.reshape(4, 2)
    x = np.exp(-N)
    return np.array(
        [
            v,
            -3 * v - x * x * q - g * (w + 3 * s),
            w,
            -3 * w - (x * x + M * M) * s + g * v,
        ]
    ).ravel()


def rhs_x(x, y, M, g):
    q, p, s, t = y.reshape(4, 2)
    return np.array(
        [
            g * s / x - x * x * p,
            q / x**2,
            -x * x * t,
            s / x**2 + (M * M + g * g) * s / x**4 - g * p / x,
        ]
    ).ravel()


def modes(
    M=2.0,
    g=0.1,
    X=400.0,
    xf=1e-4,
    rtol=1e-10,
    atol=1e-12,
    form="N",
    state="product",
    order=6,
):
    u = initial(X, M, g, state, order)
    if form == "N":
        sol = solve_ivp(
            lambda n, y: rhs_N(n, y, M, g),
            (-np.log(X), -np.log(xf)),
            u.ravel(),
            method="DOP853",
            rtol=rtol,
            atol=atol,
            dense_output=True,
        )
    else:
        c = u.copy()
        c[1] = (u[1] + g * u[2]) / X**3
        c[3] = u[3] / X**3
        sol = solve_ivp(
            lambda x, y: rhs_x(x, y, M, g),
            (X, xf),
            c.ravel(),
            method="DOP853",
            rtol=rtol,
            atol=atol,
            dense_output=True,
        )
    if not sol.success:
        raise RuntimeError(sol.message)
    return {
        "sol": sol,
        "M": M,
        "g": g,
        "X": X,
        "xf": xf,
        "rtol": rtol,
        "atol": atol,
        "form": form,
        "state": state,
        "adiabatic_order": order,
    }


def evaluate(m, x, canonical=True):
    x = np.atleast_1d(np.asarray(x, float))
    if np.min(x) < m["xf"] * (1 - 1e-10) or np.max(x) > m["X"] * (1 + 1e-10):
        raise ValueError("Outside integrated mode interval")
    y = m["sol"].sol(-np.log(x) if m["form"] == "N" else x).reshape(4, 2, -1)
    if m["form"] == "N" and canonical:
        y = y.copy()
        y[1] = (y[1] + m["g"] * y[2]) / x**3
        y[3] = y[3] / x**3
    return y


def diagnostics(m, x):
    x = np.atleast_1d(x)
    u = evaluate(m, x)
    out = []
    for j, z in enumerate(x):
        a = u[:, :, j].astype(np.clongdouble)
        raw = a @ a.conj().T
        comm = float(np.max(np.abs(2 * raw.imag - OMEGA)))
        a[1] += a[0] / z
        # Local canonical scalings balance the full covariance for eigensolving.
        v = (a @ a.conj().T).real
        scales = np.ones(4, dtype=np.longdouble)
        for k in [0, 2]:
            b = (v[k + 1, k + 1] / v[k, k]) ** 0.25
            scales[k] = b
            scales[k + 1] = 1 / b
        vb = np.asarray(v * scales[:, None] * scales[None, :], float)
        eig = np.sort(np.abs(np.linalg.eigvals(1j * OMEGA @ vb)))[::2]
        det = v[0, 0] * v[1, 1] - v[0, 1] ** 2
        nu = float(np.sqrt(max(det, 0)))
        S = float(entropy(nu)) if nu >= 0.5 - 1e-7 else np.nan
        out.append(
            [z, float(raw[0, 0].real), nu, S, comm, float(np.max(np.abs(eig - 0.5)))]
        )
    return np.array(out)


def mesh(zf, zmax, r, refine=1, order=12):
    t, w = leggauss(order)
    lo = np.linspace(np.log(zf), 0, int(np.ceil(-np.log(zf) / (0.08 / refine))) + 1)
    hi = np.linspace(
        1, zmax, int(np.ceil((zmax - 1) / (np.pi / (4 * (2 + r)) / refine))) + 1
    )
    z = []
    weights = []
    for edges, logarithmic in [(lo, True), (hi, False)]:
        left = edges[:-1, None]
        right = edges[1:, None]
        x = (left + right) / 2 + (right - left) / 2 * t
        ws = (right - left) / 2 * w
        if logarithmic:
            ws = ws * np.exp(x)
            x = np.exp(x)
        z.append(x.ravel())
        weights.append(ws.ravel())
    return np.concatenate(z), np.concatenate(weights)


def integral(
    m,
    r=0.03,
    zf=0.01,
    zmax=320.0,
    epsilons=(0.1, 0.05, 0.025),
    lam=0.01,
    refine=1,
    order=12,
):
    z, w = mesh(zf, zmax, r, refine, order)
    qf = evaluate(m, [zf])[0, :, 0].astype(np.clongdouble)
    ql = evaluate(m, [r * zf])[0, :, 0].astype(np.clongdouble)
    sig = evaluate(m, z)[2].astype(np.clongdouble)
    sigl = evaluate(m, r * z)[2].astype(np.clongdouble)
    cs = np.sum(qf[:, None] * sig.conj(), axis=0)
    cl = np.sum(ql[:, None] * sigl.conj(), axis=0)
    complex_integrand = cl * cs * cs / z**4
    pl = float(np.sum(np.abs(ql) ** 2))
    ps = float(np.sum(np.abs(qf) ** 2))
    values = []
    for eps in epsilons:
        damp = np.exp(-eps * (2 + r) * z)
        value = (
            2
            * lam
            * np.sum(w * damp * complex_integrand.imag, dtype=np.longdouble)
            / (pl * ps)
        )
        values.append(float(value))
    return np.array(values), {
        "z": z,
        "weights": w,
        "integrand": np.asarray(complex_integrand, complex),
        "power_long_dimensionless": pl,
        "power_short_dimensionless": ps,
    }


def extrapolate(eps, values):
    eps = np.asarray(eps)
    values = np.asarray(values)
    # Same last three points for the two extrapolation orders.
    lin = np.polyfit(eps[-3:], values[-3:], 1)[-1]
    quad = np.polyfit(eps[-3:], values[-3:], 2)[-1]
    return float(quad), float(abs(quad - lin))


def adaptive_check(m, r=0.03, zf=0.01, zmax=640, epsilon=0.025, lam=0.01):
    qf = evaluate(m, [zf])[0, :, 0].astype(np.clongdouble)
    ql = evaluate(m, [r * zf])[0, :, 0].astype(np.clongdouble)
    norm = 2 * lam / (float(sum(abs(qf) ** 2)) * float(sum(abs(ql) ** 2)))

    def f(z):
        s = evaluate(m, [z])[2, :, 0].astype(np.clongdouble)
        sl = evaluate(m, [r * z])[2, :, 0].astype(np.clongdouble)
        cs = sum(qf * s.conj())
        cl = sum(ql * sl.conj())
        return float((cl * cs * cs).imag / z**4 * np.exp(-epsilon * (2 + r) * z))

    val, err = quad(
        lambda y: np.exp(y) * f(np.exp(y)),
        np.log(zf),
        0,
        epsabs=1e-15,
        epsrel=1e-8,
        limit=200,
    )
    edges = np.linspace(1, zmax, int(np.ceil((zmax - 1) / (np.pi / (4 * (2 + r))))) + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        v, e = quad(f, a, b, epsabs=1e-16, epsrel=1e-8, limit=100)
        val += v
        err += e
    return norm * val, norm * err


def design(r, mu, extra=False):
    cols = [
        np.ones(len(r)),
        r * r,
        r**1.5 * np.cos(mu * np.log(r)),
        r**1.5 * np.sin(mu * np.log(r)),
    ]
    if extra:
        cols.append(r**4)
    return np.array(cols).T


def profile(r, y, sigma, extra=False):
    r = np.asarray(r)
    y = np.asarray(y)
    sigma = np.broadcast_to(sigma, y.shape)
    grid = np.arange(0.1, 4.0001, 0.01)

    def solve(mu):
        A = design(r, mu, extra) / sigma[:, None]
        b = y / sigma
        c = np.linalg.lstsq(A, b, rcond=None)[0]
        return float(np.sum((A @ c - b) ** 2)), c

    chi = np.array([solve(mu)[0] for mu in grid])
    candidates = []
    for j in range(len(grid)):
        if (j == 0 or chi[j] <= chi[j - 1]) and (
            j == len(grid) - 1 or chi[j] <= chi[j + 1]
        ):
            a = grid[max(0, j - 1)]
            b = grid[min(len(grid) - 1, j + 1)]
            opt = minimize_scalar(
                lambda m: solve(m)[0],
                bounds=(a, b),
                method="bounded",
                options={"xatol": 1e-10},
            )
            candidates.append((float(opt.fun), float(opt.x)))
    best, mu = min(candidates)
    coeff = solve(mu)[1]
    fine = np.unique(
        np.r_[grid, mu, np.linspace(max(0.1, mu - 0.05), min(4, mu + 0.05), 101)]
    )
    values = np.array([solve(m)[0] for m in fine])
    allowed = values <= best + 3.841459
    intervals = []
    start = None
    for j in range(len(fine)):
        if allowed[j] and start is None:
            start = (
                fine[j]
                if j == 0
                else brentq(
                    lambda m: solve(m)[0] - best - 3.841459, fine[j - 1], fine[j]
                )
            )
        if start is not None and (not allowed[j] or j == len(fine) - 1):
            end = (
                fine[j]
                if allowed[j]
                else brentq(
                    lambda m: solve(m)[0] - best - 3.841459, fine[j - 1], fine[j]
                )
            )
            intervals.append([float(start), float(end)])
            start = None
    return {
        "mu": mu,
        "chi2": best,
        "coefficients": coeff.tolist(),
        "intervals_95": intervals,
        "competitive_minima": [
            [m, c - best] for c, m in sorted(candidates) if c - best < 10
        ],
        "grid": fine,
        "delta_chi2": values - best,
    }


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


def exact_entropy(c):
    """Closed-form S(c), valid for |c| < sqrt(2)."""
    c = np.asarray(c, dtype=float)
    if np.any(np.abs(c) >= np.sqrt(2)):
        raise ValueError("A stable ground state requires |c| < sqrt(2)")
    d = np.sqrt(2 - c * c)
    nu = np.sqrt((2 + d) * (1 + d) / (4 * d * (3 + 2 * d)))
    return nu, entropy(nu)


def wavefunction_entropy(c):
    """Trace the ground-state Gaussian; use its geometric reduced spectrum."""
    w, U = np.linalg.eigh([[1.0, c], [c, 2.0]])
    if np.min(w) <= 0:
        raise ValueError("K must be positive")
    A = (U * np.sqrt(w)) @ U.T
    purity = float(np.sqrt(np.linalg.det(A) / (A[0, 0] * A[1, 1])))
    # Eigenvalues of the reduced Gaussian kernel: p_n=(1-t)*t**n.
    t = max(0.0, (1 - purity) / (1 + purity))
    S = -np.log1p(-t) - (t / (1 - t)) * np.log(t) if t > 0 else 0.0
    return purity, float(S)


def write_csv(path, header, rows):
    np.savetxt(
        path, np.asarray(rows), delimiter=",", header=header, comments="", fmt="%.14g"
    )


def reference(out):
    """Sections 2–3: dilation reference and the two-oscillator covariance benchmark."""
    k = np.geomspace(0.001, 1000, 121)
    ns = 0.965
    rows = []
    for h in [0.01, 0.005, 0.0025, 0.00125, 0.0001, 0.00001]:
        slope = k_derivative_response(k, ns, h)
        rows.append((h, abs(slope - (1 - ns)).max()))
    write_csv(
        out / "derivative_convergence.csv", "fractional_step,max_absolute_error", rows
    )
    dilations = dilation_response(k, ns, 0.0001)
    write_csv(
        out / "squeezed_reference.csv",
        "k,analytic,derivative,dilation",
        np.column_stack((k, np.full_like(k, 1 - ns), slope, dilations)),
    )
    tol = 5e-8
    assert np.max(abs(slope - (1 - ns))) < tol
    assert np.max(abs(dilations - (1 - ns))) < tol
    assert np.max(abs(dilations - slope)) < tol
    assert 3.8 < rows[0][1] / rows[1][1] < 4.2
    for tilt in [0.94, 0.965, 1.0, 1.02]:
        assert np.max(abs(dilation_response(k, tilt, 0.0001) - (1 - tilt))) < tol
    assert np.max(abs(dilation_response(k, ns, 0.0001, 7, 2) - dilations)) < tol
    x = np.geomspace(100, 0.01, 161)
    eta = -x  # k=H=1
    scale_factor = -1 / eta
    assert np.allclose(x, 1 / scale_factor) and np.isclose(scale_factor[len(x) // 2], 1)
    curve = []
    for c in np.linspace(0, 1.4, 701):
        nu, s_exact = exact_entropy(c)
        purity, s_kernel = wavefunction_entropy(c)
        s_cov = oscillator(c)[3]
        curve.append((c, nu, s_exact, s_cov, purity, s_kernel))
    write_csv(
        out / "oscillator_entropy.csv",
        "c,nu,exact_entropy,covariance_entropy,purity,kernel_entropy",
        curve,
    )
    samples = np.array(curve)
    assert np.max(abs(samples[:, 2] - samples[:, 3])) < 1e-12
    assert np.max(abs(samples[:, 2] - samples[:, 5])) < 1e-12
    for coupling in [0, 0.3, 0.6]:
        position, momentum, nu, _ = oscillator(coupling)
        for scale in [0.1, 0.5, 1, 2, 10]:
            determinant = (position[0, 0] * scale**2) * (momentum[0, 0] / scale**2)
            assert abs(determinant - nu**2) < 1e-13
    assert abs(oscillator(0.3)[3] - 0.0193599304) < 1e-9
    print(f"Single-field reference {1 - ns:.6f}; S(c=0.3) = {oscillator(0.3)[3]:.10f}")


def gaussian(out):
    """Sections 4–6: compare the evolution equations, free light power and entropy."""
    grid = np.geomspace(100, 0.001, 100)
    comparison = []
    for rtol, atol in [(1e-10, 1e-12), (1e-11, 1e-13)]:
        vel = modes(X=100, xf=0.001, rtol=rtol, atol=atol)
        can = modes(X=100, xf=0.001, rtol=rtol, atol=atol, form="x")
        a, b = evaluate(vel, grid), evaluate(can, grid)
        points = [1, 0.1, 0.01, 0.001]
        comparison.append(
            (
                rtol,
                atol,
                np.max(abs(a - b) / (1 + abs(a))),
                np.max(
                    abs(diagnostics(vel, points)[:, 3] - diagnostics(can, points)[:, 3])
                ),
            )
        )
    write_csv(
        out / "formulation_check.csv",
        "rtol,atol,max_scaled_difference,max_entropy_difference",
        comparison,
    )
    free = []
    for X, rtol, atol in [
        (100, 1e-10, 1e-12),
        (200, 1e-10, 1e-12),
        (400, 1e-10, 1e-12),
        (100, 1e-11, 1e-13),
    ]:
        model = modes(g=0, X=X, xf=0.001, rtol=rtol, atol=atol, form="x")
        d = diagnostics(model, np.geomspace(X, 0.001, 600))
        free.append((X, rtol, atol, d[-1, 1], d[-1, 2], d[-1, 3], d[:, 4].max()))
    write_csv(
        out / "free_modes.csv",
        "X,rtol,atol,power,nu,entropy,max_commutator_residual",
        free,
    )
    times = np.geomspace(100, 1e-4, 2000)
    history = modes(X=100, xf=1e-4, rtol=2.3e-14, atol=1e-15, form="x")
    baseline = diagnostics(history, times)
    write_csv(
        out / "entropy_history.csv",
        "x,power,nu,entropy,commutator,purity_error",
        baseline,
    )
    shear = []
    for x in [1, 0.1, 0.01, 0.001, 0.0001]:
        q, p = evaluate(history, [x])[:2, :, 0].astype(np.clongdouble)

        def determinant(momentum):
            return (
                sum(abs(q) ** 2) * sum(abs(momentum) ** 2)
                - np.real(q @ momentum.conj()) ** 2
            )

        shear.append((x, determinant(p), determinant(p + q / x)))
    write_csv(
        out / "canonical_shear.csv", "x,original_determinant,sheared_determinant", shear
    )
    sign = []
    for g in [-0.1, 0, 0.1]:
        d = diagnostics(
            modes(g=g, X=100, xf=1e-4, rtol=2.3e-14, atol=1e-15, form="x"), times
        )
        sign.append((g, np.nanmax(d[:, 3]), d[-1, 3], np.max(d[:, 4])))
    write_csv(
        out / "mixing_symmetry.csv",
        "g,max_entropy,final_entropy,max_commutator_residual",
        sign,
    )


def mass_mixing_grid(out):
    """Figure 10: the full 25 by 25 scan at common x=0.01. This is the longest run."""
    rows = []
    for mass in np.linspace(1.6, 2.4, 25):
        for mixing in np.linspace(0, 0.2, 25):
            model = modes(
                mass,
                mixing,
                X=400,
                xf=0.01,
                rtol=2.5e-14,
                atol=1e-17,
                state="adiabatic",
            )
            x, power, nu, s, comm, pure = diagnostics(model, [0.01])[0]
            good = comm < 1e-6 and pure < 1e-5 and np.isfinite(s)
            rows.append(
                (
                    mass,
                    mixing,
                    x,
                    s,
                    power / (0.5 * (1 + x * x)) - 1,
                    comm,
                    pure,
                    int(good),
                )
            )
        print(f"Completed entropy scan at M={mass:.4f}", flush=True)
    write_csv(
        out / "mass_mixing_grid.csv",
        "M,g,x,entropy,power_change,commutator,purity_error,reliable",
        rows,
    )


def benchmark(out):
    """Section 7: independent sine quadrature of the regulated analytic benchmark."""
    rows = []
    for r in [1.0, 0.1, 0.03, 0.01]:
        K = 2 + r
        exact = -1 / (2 * r * K**3)
        for eps in [0.1, 0.05, 0.025, 0.01, 0.001]:
            val, err = quad(
                lambda y: y * y * np.exp(-eps * y),
                0,
                np.inf,
                weight="sin",
                wvar=1.0,
                epsabs=1e-10,
                limit=500,
            )
            numerical = val / (4 * r * K**3)
            reg = -np.imag(2 / (eps + 1j) ** 3) / (4 * r * K**3)
            assert abs(numerical - reg) < 1e-8 * abs(exact)
            rows.append([r, eps, exact, numerical, abs(numerical / exact - 1), err])
    rows = np.array(rows)
    assert np.max(rows[rows[:, 1] == 0.001, 4]) < 0.001
    write_csv(
        out / "in_in_benchmark.csv",
        "r,epsilon,exact_bispectrum,numerical_bispectrum,relative_error,quadrature_error",
        rows,
    )


def pilot(out):
    """Section 8: regulated pilot, zero-coupling and weak-mixing tests."""
    rows = []
    mm = modes(X=800, xf=1e-5, rtol=1e-12, atol=1e-15, state="adiabatic")
    for r in [0.1, 0.03, 0.01]:
        vals, raw = integral(mm, r, zmax=320, epsilons=[0.025])
        rows.append([r, vals[0]])
        if r == 0.03:
            cumulative = (
                2
                * 0.01
                * np.cumsum(
                    raw["weights"]
                    * np.exp(-0.025 * (2 + r) * raw["z"])
                    * raw["integrand"].imag
                )
                / (raw["power_long_dimensionless"] * raw["power_short_dimensionless"])
            )
            write_csv(
                out / "cumulative_integral.csv",
                "z,cumulative_F",
                np.c_[raw["z"], cumulative],
            )
    controls = []
    for g in [0.0, 0.025, 0.05, 0.1]:
        model = (
            mm
            if g == 0.1
            else modes(g=g, X=800, xf=1e-5, rtol=1e-12, atol=1e-15, state="adiabatic")
        )
        val, _ = integral(model, 0.03, zmax=320, epsilons=[0.025])
        controls.append([g, val[0], val[0] / g**3 if g else 0])
    write_csv(out / "regulated_bispectrum.csv", "r,regulated_F_epsilon_0025", rows)
    write_csv(out / "mixing_controls.csv", "g,F,F_over_g_cubed", controls)
    z, _ = integral(mm, 0.03, zmax=320, epsilons=[0.025], lam=0)
    twice, _ = integral(mm, 0.03, zmax=320, epsilons=[0.025], lam=0.02)
    assert np.isclose(twice[0] / controls[-1][1], 2.0)
    assert controls[0][1] == 0 and z[0] == 0


def product_control(p):
    rows = []
    for X in [800, 1600]:
        for state in ["product", "adiabatic"]:
            m = modes(X=X, xf=1e-5, rtol=1e-12, atol=1e-15, state=state)
            vals, _ = integral(m, 0.03, zmax=640, epsilons=[0.1, 0.05, 0.025, 0.0125])
            f, e = extrapolate([0.1, 0.05, 0.025, 0.0125], vals)
            rows.append([X, int(state == "adiabatic"), f, e])
    write_csv(
        p / "initial_state_control.csv",
        "X,adiabatic,F_extrapolated,regulator_error",
        rows,
    )


def independent_quadrature(p):
    m = modes(X=1600, xf=1e-5, rtol=1e-12, atol=1e-15, state="adiabatic")
    v, _ = integral(m, 0.03, zmax=640, epsilons=[0.025])
    a, e = adaptive_check(m)
    write_csv(
        p / "adaptive_quadrature.csv",
        "Gauss_Legendre,adaptive_Gauss_Kronrod,difference,adaptive_error_estimate",
        [[v[0], a, abs(v[0] - a), e]],
    )


def convergence(out):
    """Section 9: build the six-component error budget for the physical fits."""
    eps = np.array([0.1, 0.05, 0.025, 0.0125])
    rs = [0.2, 0.1, 0.03, 0.01, 0.003, 0.001]
    baseline = modes(X=1600, xf=1e-6, rtol=1e-12, atol=1e-15, state="adiabatic")
    other = modes(X=800, xf=1e-6, rtol=1e-12, atol=1e-15, state="adiabatic")
    loose = modes(X=1600, xf=1e-6, rtol=1e-11, atol=1e-14, state="adiabatic")
    rows = []
    summary = []
    variation_rows = []
    for r in rs:
        vals, _ = integral(baseline, r, zmax=640, epsilons=eps)
        variation_rows.extend([[r, 0, e, v] for e, v in zip(eps, vals)])
        fixed = []
        for e, L in zip([0.1, 0.05, 0.025], [80, 160, 320]):
            short, _ = integral(baseline, r, zmax=L, epsilons=[e])
            long, _ = integral(baseline, r, zmax=2 * L, epsilons=[e])
            fixed.append(abs(short[0] - long[0]))
            rows.append([r, 0, e, L, short[0], long[0]])
        central, regerr = extrapolate(eps, vals)
        vx, _ = integral(other, r, zmax=640, epsilons=eps)
        fx, _ = extrapolate(eps, vx)
        vt, _ = integral(loose, r, zmax=640, epsilons=eps)
        ft, _ = extrapolate(eps, vt)
        vq, _ = integral(baseline, r, zmax=640, epsilons=eps, refine=2, order=16)
        fq, _ = extrapolate(eps, vq)
        variation_rows.extend([[r, 1, e, v] for e, v in zip(eps, vx)])
        variation_rows.extend([[r, 2, e, v] for e, v in zip(eps, vt)])
        variation_rows.extend([[r, 3, e, v] for e, v in zip(eps, vq)])
        finals = []
        for final_index, zf in enumerate([0.003, 0.001]):
            vf, _ = integral(baseline, r, zf=zf, zmax=640, epsilons=eps)
            ff, _ = extrapolate(eps, vf)
            finals.append(ff)
            variation_rows.extend([[r, 4 + final_index, e, v] for e, v in zip(eps, vf)])
        errs = [
            abs(fx - central),
            regerr,
            abs(ft - central),
            abs(fq - central),
            max(abs(np.array(finals) - central)),
            max(fixed),
        ]
        uncertainty = sum(errs)
        summary.append([r, central, *errs, uncertainty])
        for e, v in zip(eps, vals):
            rows.append([r, 1, e, 640, v, v])
    write_csv(
        out / "variation_integrals.csv",
        "r,variation_code,epsilon,F_regulated",
        variation_rows,
    )
    write_csv(
        out / "cutoff_and_regulator.csv",
        "r,test_type,epsilon,zmax,F_first,F_doubled_cutoff",
        rows,
    )
    write_csv(
        out / "error_budget.csv",
        "r,F_extrapolated,initial_X,regulator,ODE_tolerance,quadrature,final_time,early_tail,total_conservative_error",
        summary,
    )
    summary = np.array(summary)

    ratio = max(summary[:, -1] / summary[:, 0] ** 1.5) / max(
        abs(summary[:, 1]) / summary[:, 0] ** 1.5
    )
    print(f"Pilot envelope error = {ratio:.6g}")
    if ratio >= 0.01:
        raise RuntimeError("Physical frequency scan requires pilot error below 1%.")
    product_control(out)
    independent_quadrature(out)
    return summary


def clock_fit(out, budget):
    """Section 10: compute the 50-triangle physical shape and all 12 mass fits."""
    r = np.geomspace(0.001, 0.2, 50)
    model = modes(X=1600, xf=1e-6, rtol=1e-12, atol=1e-15, state="adiabatic")
    vals = []
    errors = []
    regulated = []
    budget = budget[np.argsort(budget[:, 0])]
    other_error = budget[:, -1] - budget[:, 3]
    eps = [0.1, 0.05, 0.025, 0.0125]
    for rr in r:
        vv, _ = integral(model, float(rr), zf=0.01, zmax=640, epsilons=eps)
        regulated.extend([[rr, ee, v0] for ee, v0 in zip(eps, vv)])
        f, e = extrapolate(eps, vv)
        vals.append(f)
        extra = float(np.interp(np.log(rr), np.log(budget[:, 0]), other_error))
        errors.append(max(e + extra, 1e-14))
    vals = np.array(vals)
    errors = np.array(errors)
    write_csv(
        out / "physical_bispectrum.csv",
        "r,F,conservative_numerical_error",
        np.c_[r, vals, errors],
    )
    write_csv(out / "physical_regulated.csv", "r,epsilon,F_regulated", regulated)
    fits = []
    for upper in [0.05, 0.1, 0.2]:
        for lower in [0.001, 0.003]:
            for extra in [False, True]:
                use = (r >= lower) & (r <= upper)
                fit = profile(r[use], vals[use], errors[use], extra)
                mu = fit["mu"]
                mass = np.sqrt(mu**2 + 2.25 - 0.1**2)
                fits.append([lower, upper, int(extra), mu, mass, fit["chi2"]])
                if lower == r[0] and upper == r[-1] and not extra:
                    full = fit
    write_csv(
        out / "clock_fits.csv", "r_min,r_max,extra_r4,mu,mass,weighted_residual", fits
    )
    masses = np.array(fits)[:, 4]
    amplitude = np.hypot(*full["coefficients"][2:4])
    error_fraction = np.max(errors / r**1.5) / amplitude
    cycles = full["mu"] * np.log(r[-1] / r[0]) / (2 * np.pi)
    print(f"Best-fit mass envelope: {masses.min():.8f} to {masses.max():.8f}")
    print(
        f"Clock amplitude {amplitude:.8g}, error/envelope {error_fraction:.6g}, cycles {cycles:.6f}"
    )
    return full


def mock_trials(out):
    """Section 10: seeded mock-frequency and held-out no-clock trials."""
    r = np.geomspace(1e-3, 0.2, 50)
    mu = np.sqrt(1.76)
    coef = np.array([0.0, 0.2, 0.1, 0.05])
    signal = design(r, mu) @ coef
    oscill = design(r, mu)[:, 2:] @ coef[2:]
    rms = float(np.sqrt(np.mean(oscill**2)))
    rng = np.random.default_rng(20261029)
    rows = []
    nulls = {}
    coverage = {}
    for frac in [0.01, 0.05]:
        sigma = frac * rms
        cov = 0
        est = []
        null = []
        for j in range(100):
            y = signal + rng.normal(0, sigma, len(r))
            fit = profile(r, y, sigma)
            inside = any(a <= mu <= b for a, b in fit["intervals_95"])
            cov += inside
            est.append(fit["mu"])
            rows.append(
                [frac, j, fit["mu"], np.sqrt(fit["mu"] ** 2 + 2.25 - 0.01), int(inside)]
            )
            yn = design(r, mu)[:, :2] @ coef[:2] + rng.normal(0, sigma, len(r))
            nf = profile(r, yn, sigma)
            A = design(r, mu)[:, :2] / sigma
            bb = yn / sigma
            chi0 = np.sum((A @ np.linalg.lstsq(A, bb, rcond=None)[0] - bb) ** 2)
            null.append(float(chi0 - nf["chi2"]))
            if j == 0:
                write_csv(
                    out / f"mock_profile_{int(frac * 100)}pct.csv",
                    "mu,delta_chi2",
                    np.c_[fit["grid"], fit["delta_chi2"]],
                )
        threshold = float(np.quantile(null, 0.95))
        held = np.random.default_rng(20261030 + int(frac * 100))
        false_count = 0
        for j in range(100):
            yn = design(r, mu)[:, :2] @ coef[:2] + held.normal(0, sigma, len(r))
            nf = profile(r, yn, sigma)
            A = design(r, mu)[:, :2] / sigma
            bb = yn / sigma
            chi0 = np.sum((A @ np.linalg.lstsq(A, bb, rcond=None)[0] - bb) ** 2)
            false_count += chi0 - nf["chi2"] > threshold
        nulls[str(frac)] = {
            "threshold_95": threshold,
            "calibration_false_alarm_rate": float(np.mean(np.array(null) > threshold)),
            "held_out_false_alarm_rate": false_count / 100,
            "held_out_trials": 100,
        }
        coverage[str(frac)] = {
            "coverage": cov / 100,
            "frequency_bias": float(np.mean(est) - mu),
            "frequency_RMSE": float(np.sqrt(np.mean((np.array(est) - mu) ** 2))),
        }
    write_csv(
        out / "mock_recovery.csv",
        "noise_fraction,trial,mu_estimate,mass_estimate,interval_contains_truth",
        rows,
    )

    for noise in coverage:
        print(f"Mock noise {noise}: {coverage[noise]}, null tests {nulls[noise]}")


def cubic_controls(out, full):
    """Section 11: fixed-entropy controls and known-template signal-to-noise."""
    rows = []
    for M, g in [(2.0, 0.05), (2.0, 0.1), (2.0, 0.2), (1.8, 0.1), (2.4, 0.1)]:
        mm = modes(M, g, X=1600, xf=1e-5, rtol=1e-12, atol=1e-15, state="adiabatic")
        d = diagnostics(mm, [0.01])[0]
        v, _ = integral(mm, 0.03, zmax=640, epsilons=[0.1, 0.05, 0.025, 0.0125])
        f, e = extrapolate([0.1, 0.05, 0.025, 0.0125], v)
        rows.append([M, g, d[3], d[1] / (0.5 * (1 + 0.01**2)) - 1, f, e])
    write_csv(
        out / "entropy_signal_comparison.csv",
        "M,g,entropy_x001,power_change_x001,F_r003,regulator_error",
        rows,
    )
    baseline = [a for a in rows if a[0] == 2 and a[1] == 0.1][0]
    cubic = [
        [l, baseline[2], baseline[4] * l / 0.01, baseline[5] * l / 0.01]
        for l in [0, 0.005, 0.01]
    ]
    write_csv(
        out / "cubic_control.csv",
        "lambda3_over_H,Gaussian_entropy,F_r003,numerical_error",
        cubic,
    )
    rr = np.geomspace(0.001, 0.2, 50)
    mu = full["mu"]
    # Fixed-covariance diagnostic using the fitted baseline clock template.
    A, B = full["coefficients"][2:4]
    template = rr**1.5 * (A * np.cos(mu * np.log(rr)) + B * np.sin(mu * np.log(rr)))
    noise_floor = 1e-10
    recover = [
        [l, float(np.linalg.norm(template * l / 0.01) / noise_floor), noise_floor]
        for l in [0, 0.005, 0.01]
    ]
    write_csv(
        out / "common_noise_control.csv",
        "lambda3_over_H,known_template_SNR,absolute_sigma_F",
        recover,
    )

    print(f"Baseline entropy {baseline[2]:.10g}, F(0.03) {baseline[4]:.10g}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "calculation",
        nargs="?",
        default="all",
        choices=[
            "reference",
            "gaussian",
            "benchmark",
            "pilot",
            "mock",
            "clock",
            "grid",
            "all",
        ],
        help="clock runs convergence, the physical frequency fits and the cubic controls",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results"),
        help="directory for generated tables",
    )
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    tasks = {
        "reference": reference,
        "gaussian": gaussian,
        "benchmark": benchmark,
        "pilot": pilot,
        "mock": mock_trials,
    }
    for name, run in tasks.items():
        if args.calculation in (name, "all"):
            print(f"Running {name}", flush=True)
            run(args.out)
    if args.calculation in ("clock", "all"):
        budget = convergence(args.out)
        fitted = clock_fit(args.out, budget)
        cubic_controls(args.out, fitted)
    if args.calculation in ("grid", "all"):
        mass_mixing_grid(args.out)


if __name__ == "__main__":
    main()
