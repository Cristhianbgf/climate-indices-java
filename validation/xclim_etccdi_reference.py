"""Reference ETCCDI precipitation indices from xclim for the Java tests.

Reads src/test/resources/era5_daily_pr.csv and writes src/test/resources/xclim_etccdi.csv,
with one value per calendar year (1960-2025) for each point and index:

  PRCPTOT  total precipitation on wet days (>= 1 mm)    xclim.indices.prcptot(thresh="1 mm/d")
  RX1DAY   highest 1-day precipitation                   max_1day_precipitation_amount
  RX5DAY   highest 5-day precipitation                   max_n_day_precipitation_amount(window=5)
  SDII     mean precipitation on wet days (>= 1 mm)      daily_pr_intensity
  R10MM    days with at least 10 mm                      wetdays(thresh="10 mm/d")
  R20MM    days with at least 20 mm                      wetdays(thresh="20 mm/d")

Run from the repository root:  python validation/xclim_etccdi_reference.py
Tested with xclim 0.61.1 and xarray 2026.7.
"""
from pathlib import Path

import pandas as pd
import xarray as xr
import xclim
from xclim import indices as xi

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "src" / "test" / "resources"

daily = pd.read_csv(RES / "era5_daily_pr.csv", index_col="date")
time = pd.DatetimeIndex(pd.to_datetime(daily.index).to_numpy(), name="time")

INDICES = {
    "prcptot": lambda pr: xi.prcptot(pr, thresh="1 mm/d"),
    "rx1day": xi.max_1day_precipitation_amount,
    "rx5day": lambda pr: xi.max_n_day_precipitation_amount(pr, window=5),
    "sdii": lambda pr: xi.daily_pr_intensity(pr, thresh="1 mm/d"),
    "r10mm": lambda pr: xi.wetdays(pr, thresh="10 mm/d"),
    "r20mm": lambda pr: xi.wetdays(pr, thresh="20 mm/d"),
}

out = {}
for col in daily.columns:
    pr = xr.DataArray(daily[col].to_numpy(), coords={"time": time}, dims="time", attrs={"units": "mm/d"})
    for name, f in INDICES.items():
        v = f(pr)
        out[f"{name}_{col.removeprefix('pr_')}"] = pd.Series(v.to_numpy(), index=v.time.dt.year.to_numpy())

df = pd.DataFrame(out)
df.index.name = "year"
df.to_csv(RES / "xclim_etccdi.csv", float_format="%.10g", na_rep="NaN", lineterminator="\n")
print(f"xclim {xclim.__version__}: {df.shape[1]} series, {df.index[0]}-{df.index[-1]}")
