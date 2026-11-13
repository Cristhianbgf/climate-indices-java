"""Reference R95pTOT and R99pTOT for the Java tests, from an implementation written here and from icclim and xclim.

Reads src/test/resources/era5_daily_pr.csv and writes two files:

  etccdi_percentile_thresholds.csv  per point, the 95th and 99th percentiles of wet days (>= 1 mm) in the
                                     base period 1961-1990, by two quantile methods:
                                       median_unbiased  Hyndman and Fan type 8, icclim's default
                                       linear           type 7, the numpy and R default
  etccdi_percentiles.csv            per calendar year (1960-2025), point and variant:
    own_<index>_mm_<method>    written below with numpy: total precipitation of the days above the threshold
    own_<index>_pct_<method>   the same as a percentage of PRCPTOT (days >= 1 mm)
    icclim_<index>_pct         icclim's ECA&D R95pTOT and R99pTOT (a percentage), base 1961-1990
    icclim_<index>_mm          icclim's generic sum with the same threshold
    xclim_<index>_pct_gt       100 x xclim fraction_over_precip_thresh, condition ">"
    xclim_<index>_pct_ge       the same with condition ">="
  where <index> is r95ptot or r99ptot. A day counts when it is above the threshold (strictly).

xclim takes the threshold from the caller; it gets the type 8 percentiles computed here. Its condition
applies to both the wet-day test and the percentile test, so neither choice matches the ETCCDI definition
exactly (wet: >= 1 mm; extreme: > percentile).

Run from the repository root, in an environment with xclim and icclim:
  python validation/etccdi_percentiles_reference.py
Tested with icclim 7.3.7, xclim 1.0.0 and xarray 2026.9.
"""
from pathlib import Path

import icclim
import numpy as np
import pandas as pd
import xarray as xr
import xclim
from xclim import compute as xi

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "src" / "test" / "resources"
BASE = ("1961-01-01", "1990-12-31")
WET = 1.0
LEVELS = {"r95ptot": 95, "r99ptot": 99}
METHODS = ("median_unbiased", "linear")

daily = pd.read_csv(RES / "era5_daily_pr.csv", index_col="date")
dates = pd.to_datetime(daily.index)
years = dates.year.to_numpy()
in_base = (dates >= BASE[0]) & (dates <= BASE[1])
time = pd.DatetimeIndex(dates.to_numpy(), name="time")

annual, thresholds = {}, {}
for col in daily.columns:
    point = col.removeprefix("pr_")
    pr = daily[col].to_numpy()
    wet_base = pr[in_base & (pr >= WET)]
    prcptot = pd.Series(np.where(pr >= WET, pr, 0.0)).groupby(years).sum()

    row = {}
    for index, p in LEVELS.items():
        for method in METHODS:
            q = float(np.quantile(wet_base, p / 100, method=method))
            row[f"q{p}_{method}"] = q
            over = pd.Series(np.where(pr > q, pr, 0.0)).groupby(years).sum()
            annual[f"own_{index}_mm_{method}_{point}"] = over
            annual[f"own_{index}_pct_{method}_{point}"] = 100 * over / prcptot
    thresholds[point] = row

    da = xr.DataArray(pr, coords={"time": time}, dims="time", attrs={"units": "mm/d", "standard_name": "precipitation_flux"})
    ds = xr.Dataset({"pr": da})
    for index, p in LEVELS.items():
        v = icclim.index(index_name=f"R{p}pTOT", in_files=ds, slice_mode="year", base_period_time_range=list(BASE)).to_array().squeeze()
        annual[f"icclim_{index}_pct_{point}"] = pd.Series(v.to_numpy(), index=v.time.dt.year.to_numpy())

        thr = icclim.build_threshold(f"> {p} period_per", threshold_min_value="1 mm/day", reference_period=list(BASE))
        v = icclim.sum(in_files=ds, threshold=thr, slice_mode="year").to_array().squeeze()
        annual[f"icclim_{index}_mm_{point}"] = pd.Series(v.to_numpy(), index=v.time.dt.year.to_numpy())

        per = xr.DataArray(row[f"q{p}_median_unbiased"], attrs={"units": "mm/d"})
        for name, cond in (("gt", ">"), ("ge", ">=")):
            v = xi.fraction_over_precip_thresh(da, per, condition=cond) * 100
            annual[f"xclim_{index}_pct_{name}_{point}"] = pd.Series(v.to_numpy(), index=v.time.dt.year.to_numpy())

df = pd.DataFrame(annual)
df.index.name = "year"
df.to_csv(RES / "etccdi_percentiles.csv", float_format="%.12g", na_rep="NaN", lineterminator="\n")
th = pd.DataFrame(thresholds).T
th.index.name = "point"
th.to_csv(RES / "etccdi_percentile_thresholds.csv", float_format="%.12g", lineterminator="\n")
print(f"xclim {xclim.__version__}, icclim {icclim.__version__}: {df.shape[1]} series, {df.index[0]}-{df.index[-1]}")
print(th)
