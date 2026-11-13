"""Compare the maximum likelihood and the L-moments (PWM) fits for SPI and SPEI, by basin, window and month.

The four indices come from reference files that the Java code matches (see validation/README.md):

  SPI, maximum likelihood    xclim_spi.csv (xclim 0.61.1, gamma, method="ML")
  SPI, PWM                   computed below with numpy and scipy, written independently of the Java
                             code: gamma fitted to the non-zero values with Hosking's PELGAM, the
                             share of zeros added as in xclim and in the Java SPI. It is written to
                             src/test/resources/spi_pwm.csv for the Java tests.
  SPEI, maximum likelihood   computed below from the generalized logistic parameters in
                             glo_ml_params.csv (independent of the Java fit)
  SPEI, PWM                  r_spei.csv (R package SPEI 1.8.1)

Values beyond the fitted support are set to +/-8.21, as in the Java code. The calibration period is
1991-2020. Writes validation/fit_comparison.csv with, per index, basin, window and calendar month:
the largest and the mean absolute difference between the two fits inside and outside the calibration
period, how many values differ by more than 0.1, and, inside the calibration period, how many months
change class (boundaries at -2, -1.5, -1, 1, 1.5 and 2), how many of those lie on the dry side (mean
of the two values below 0) and how many cross -1, the usual start of a drought. It also prints, per
basin and window, the mean difference where both fits give at least 1.5 in absolute value over the
mean difference where both stay within 1.

Run from the repository root:  python validation/fit_comparison.py
Tested with numpy 2.5.4, scipy 1.18.1 and pandas 3.1.1.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as st

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "src" / "test" / "resources"
BASINS = ("1391", "132", "4994")
WINDOWS = (1, 3, 12)
LIMIT = 8.21
HOSKING_PI = 3.1415927  # as written in Hosking's PELGAM, which lmom uses

ppt = pd.read_csv(RES / "basin_ppt.csv", index_col="month")
pet = pd.read_csv(RES / "basin_pet.csv", index_col="month")
months = pd.to_datetime(ppt.index + "-01")
cal = np.asarray((months.year >= 1991) & (months.year <= 2020))
month_of = np.asarray(months.month)


def window_mean(values, window):
    out = np.full(values.shape, np.nan)
    for i in range(window - 1, len(values)):
        out[i] = values[i - window + 1 : i + 1].sum() / window
    return out


def to_normal(p):
    return np.clip(st.norm.ppf(p), -LIMIT, LIMIT)


def drought_class(v):
    return np.digitize(v, [-2, -1.5, -1, 1, 1.5, 2])


def pelgam(x):
    """Gamma shape and scale from the unbiased PWMs of positive values (Hosking's rational approximation)."""
    x = np.sort(x)
    n = len(x)
    if n < 4:
        return np.nan, np.nan
    j = np.arange(n)
    b0 = x.mean()
    b1 = (x * j / (n - 1)).sum() / n
    cv = (2 * b1 - b0) / b0
    if cv >= 0.5:
        t = 1 - cv
        shape = t * (0.7213 + t * -0.5947) / (1 + t * (-2.1817 + t * 1.2113))
    else:
        t = HOSKING_PI * cv * cv
        shape = (1 + -0.3080 * t) / (t * (1 + t * (-0.05812 + t * 0.01765)))
    return shape, b0 / shape


def spi_pwm(pr, window):
    acc = window_mean(pr, window)
    out = np.full(acc.shape, np.nan)
    for m in range(1, 13):
        sample = acc[(month_of == m) & cal]
        sample = sample[np.isfinite(sample)]
        p0 = np.mean(sample == 0)
        shape, scale = pelgam(sample[sample > 0])
        idx = (month_of == m) & np.isfinite(acc)
        x = acc[idx]
        p = np.where(x > 0, p0 + (1 - p0) * st.gamma.cdf(x, shape, scale=scale), p0)
        out[idx] = to_normal(p)
    return out


def glo_cdf(x, xi, alpha, kappa):
    u = (x - xi) / alpha
    if kappa == 0:
        y = u
    else:
        z = 1 - kappa * u
        with np.errstate(invalid="ignore", divide="ignore"):
            y = np.where(z > 0, -np.log(np.where(z > 0, z, 1)) / kappa, np.inf if kappa > 0 else -np.inf)
    return 1 / (1 + np.exp(-y))


def spei_ml(wb, window, params):
    acc = window_mean(wb, window)
    out = np.full(acc.shape, np.nan)
    for m in range(1, 13):
        p = params[(params.window == window) & (params.month == m)].iloc[0]
        idx = (month_of == m) & np.isfinite(acc)
        out[idx] = to_normal(glo_cdf(acc[idx], p.xi, p.alpha, p.kappa))
    return out


xclim_spi = pd.read_csv(RES / "xclim_spi.csv", index_col="month")
r_spei = pd.read_csv(RES / "r_spei.csv", index_col="month").clip(-LIMIT, LIMIT)
glo_ml = pd.read_csv(RES / "glo_ml_params.csv", dtype={"basin": str})

spi_ref = pd.DataFrame(index=ppt.index)
rows, ratios = [], []
for basin in BASINS:
    pr = ppt[f"ppt_{basin}"].to_numpy()
    wb = pr - pet[f"pet_{basin}"].to_numpy()
    for window in WINDOWS:
        pwm = spi_pwm(pr, window)
        spi_ref[f"spi_{basin}_w{window}"] = pwm
        pairs = {
            "spi": (xclim_spi[f"spi_{basin}_w{window}_ML"].to_numpy(), pwm),
            "spei": (spei_ml(wb, window, glo_ml[glo_ml.basin == basin]), r_spei[f"spei_{basin}_w{window}"].to_numpy()),
        }
        for index, (ml, pw) in pairs.items():
            d = np.abs(ml - pw)
            ok = cal & np.isfinite(d)
            tails = ok & (np.minimum(np.abs(ml), np.abs(pw)) >= 1.5)
            centre = ok & (np.maximum(np.abs(ml), np.abs(pw)) < 1)
            ratios.append((index, basin, window, d[tails].mean() / d[centre].mean()))
            for m in range(1, 13):
                row = {"index": index, "basin": basin, "window": window, "month": m}
                for period, sel in (("cal", cal), ("out", ~cal)):
                    v = d[(month_of == m) & sel & np.isfinite(d)]
                    row[f"max_{period}"] = v.max()
                    row[f"mean_{period}"] = v.mean()
                    row[f"over01_{period}"] = int((v > 0.1).sum())
                    row[f"n_{period}"] = len(v)
                c = (month_of == m) & cal & np.isfinite(d)
                changed = drought_class(ml[c]) != drought_class(pw[c])
                row["class_changes_cal"] = int(changed.sum())
                row["class_changes_dry_cal"] = int((changed & ((ml[c] + pw[c]) < 0)).sum())
                row["drought_flips_cal"] = int(((ml[c] <= -1) != (pw[c] <= -1)).sum())
                rows.append(row)

spi_ref.to_csv(RES / "spi_pwm.csv", float_format="%.17g", na_rep="NaN", lineterminator="\n")
table = pd.DataFrame(rows)
table.to_csv(ROOT / "validation" / "fit_comparison.csv", index=False, float_format="%.4f", lineterminator="\n")

for index in ("spi", "spei"):
    t = table[table["index"] == index]
    print(index.upper(), "calibration: max", round(t.max_cal.max(), 3), "| >0.1:", t.over01_cal.sum(), "of", t.n_cal.sum(),
          "| outside: max", round(t.max_out.max(), 3), ">0.1:", t.over01_out.sum(), "of", t.n_out.sum(),
          "| class changes:", t.class_changes_cal.sum(), "(dry side", t.class_changes_dry_cal.sum(), ")",
          "| crossing -1:", t.drought_flips_cal.sum())
    r = [x[3] for x in ratios if x[0] == index]
    print(f"  tails over centre, mean difference: {min(r):.1f} to {max(r):.1f}")
