"""Reference maximum likelihood fits of the generalized logistic (GLO) distribution, and the
comparison with xclim's log-logistic maximum likelihood fit.

For every basin, window (1, 3 and 12 months) and calendar month, the water balance of the
calibration period 1991-2020 is fitted by maximum likelihood in Hosking's GLO form (location xi,
scale alpha, shape kappa), written here independently of the Java code. Nelder-Mead runs from the
PWM fit of the R package SPEI (r_glo_params.csv) and is restarted until the log-likelihood stops
improving.

Writes src/test/resources/glo_ml_params.csv with the parameters and the log-likelihood, and
src/test/resources/xclim_fisk_ml.csv with the fit xclim makes for dist="fisk", method="ML"
(its starting values from xclim.indices.stats._fit_start, a private function, then
scipy.stats.fisk.fit), converted to GLO form, and its log-likelihood.

Run from the repository root:  python validation/glo_ml_reference.py
Tested with numpy 2.4, scipy 1.18.0 and xclim 0.61.1.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.optimize as so
import scipy.stats as st
import xclim
from xclim.indices.stats import _fit_start

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "src" / "test" / "resources"
BASINS = ("1391", "132", "4994")
WINDOWS = (1, 3, 12)


def glo_loglik(x, xi, alpha, kappa):
    if alpha <= 0 or abs(kappa) >= 1:
        return -np.inf
    u = (x - xi) / alpha
    if abs(kappa) < 1e-12:
        y = u
    else:
        z = 1 - kappa * u
        if np.any(z <= 0):
            return -np.inf
        y = -np.log1p(-kappa * u) / kappa
    # log of 1 + exp(-y), written so that it does not overflow for large negative y
    log1pexp = np.where(y > 0, np.log1p(np.exp(-y)), -y + np.log1p(np.exp(y)))
    return float(np.sum(-np.log(alpha) - (1 - kappa) * y - 2 * log1pexp))


def fit_ml(x, start):
    def nll(q):
        v = glo_loglik(x, q[0], np.exp(q[1]), np.tanh(q[2]))
        return np.inf if not np.isfinite(v) else -v

    q = np.array([start[0], np.log(start[1]), np.arctanh(start[2])])
    if not np.isfinite(nll(q)):
        q = np.array([x.mean(), np.log(x.std() * np.sqrt(3) / np.pi), 0.0])  # logistic start
    best = nll(q)
    for _ in range(50):
        r = so.minimize(nll, q, method="Nelder-Mead",
                        options={"xatol": 1e-13, "fatol": 1e-15, "maxiter": 100000, "maxfev": 200000})
        q = r.x
        if best - r.fun < 1e-13:
            break
        best = r.fun
    return q[0], np.exp(q[1]), np.tanh(q[2]), -nll(q)


def window_mean(values, window):
    out = np.full(values.shape, np.nan)
    for i in range(window - 1, len(values)):
        out[i] = values[i - window + 1 : i + 1].sum() / window
    return out


ppt = pd.read_csv(RES / "basin_ppt.csv", index_col="month")
pet = pd.read_csv(RES / "basin_pet.csv", index_col="month")
pwm = pd.read_csv(RES / "r_glo_params.csv", dtype={"basin": str})
months = pd.to_datetime(ppt.index + "-01")
cal = (months.year >= 1991) & (months.year <= 2020)

ml_rows, xclim_rows = [], []
for basin in BASINS:
    wb = (ppt[f"ppt_{basin}"] - pet[f"pet_{basin}"]).to_numpy()
    for window in WINDOWS:
        acc = window_mean(wb, window)
        for month in range(1, 13):
            x = acc[(months.month == month) & cal]
            x = x[np.isfinite(x)]
            p = pwm[(pwm.basin == basin) & (pwm.window == window) & (pwm.month == month)].iloc[0]
            # R fits sums, this script means: location and scale scale with the window
            start = (p.xi / window, p.alpha / window, p.kappa)
            xi, alpha, kappa, ll = fit_ml(x, start)
            ml_rows.append((basin, window, month, xi, alpha, kappa, ll))

            args, kwargs = _fit_start(x, "fisk")
            c, loc, scale = st.fisk.fit(x, *args, **kwargs)
            k, a, xi_f = -1 / c, scale / c, loc + scale
            xclim_rows.append((basin, window, month, c, xi_f, a, k, glo_loglik(x, xi_f, a, k)))

cols = ["basin", "window", "month", "xi", "alpha", "kappa", "loglik"]
pd.DataFrame(ml_rows, columns=cols).to_csv(RES / "glo_ml_params.csv", index=False, float_format="%.17g")
pd.DataFrame(xclim_rows, columns=["basin", "window", "month", "fisk_c"] + cols[3:]).to_csv(
    RES / "xclim_fisk_ml.csv", index=False, float_format="%.17g")
print(f"xclim {xclim.__version__}: {len(ml_rows)} fits written")
