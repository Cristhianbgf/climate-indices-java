"""Reference CDD and CWD for the Java tests, from an implementation written here and from xclim and icclim.

Reads src/test/resources/era5_daily_pr.csv and writes src/test/resources/etccdi_spells.csv with one
value per calendar year (1960-2025), point and variant:

  own_<index>_<rule>   written below with numpy, from the ETCCDI definitions: a dry day has less than
                       1 mm and a wet day at least 1 mm. A spell that crosses 31 December is
                         split  cut into its two calendar years,
                         start  counted whole in the year it starts,
                         end    counted whole in the year it ends.
  xclim_<index>_split  xclim maximum_consecutive_dry_days / _wet_days with their defaults
                       (xclim 1.0 counts a wet day at >= 1 mm; 0.61 used > 1 mm)
  xclim_<index>_start  the same with resample_before_rl=False
  icclim_<index>       icclim's ECA&D index with its defaults

Run from the repository root, in an environment with xclim and icclim:
  python validation/etccdi_spells_reference.py
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


def spells(flags):
    """(start, end) positions of the runs of True, end inclusive."""
    out, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        elif not f and start is not None:
            out.append((start, i - 1))
            start = None
    if start is not None:
        out.append((start, len(flags) - 1))
    return out


def longest_per_year(flags, years, rule):
    result = {y: 0 for y in np.unique(years)}
    for s, e in spells(flags):
        if rule == "split":
            # cut the spell at every change of year
            i = s
            while i <= e:
                y = years[i]
                j = i
                while j + 1 <= e and years[j + 1] == y:
                    j += 1
                result[y] = max(result[y], j - i + 1)
                i = j + 1
        else:
            y = years[s] if rule == "start" else years[e]
            result[y] = max(result[y], e - s + 1)
    return result


daily = pd.read_csv(RES / "era5_daily_pr.csv", index_col="date")
dates = pd.to_datetime(daily.index)
years = dates.year.to_numpy()
time = pd.DatetimeIndex(dates.to_numpy(), name="time")

out = {}
for col in daily.columns:
    point = col.removeprefix("pr_")
    pr = daily[col].to_numpy()
    for index, flags in (("cdd", pr < 1.0), ("cwd", pr >= 1.0)):
        for rule in ("split", "start", "end"):
            out[f"own_{index}_{rule}_{point}"] = pd.Series(longest_per_year(flags, years, rule))

    da = xr.DataArray(pr, coords={"time": time}, dims="time", attrs={"units": "mm/d", "standard_name": "precipitation_flux"})
    for index, f in (("cdd", xi.maximum_consecutive_dry_days), ("cwd", xi.maximum_consecutive_wet_days)):
        for rule, before in (("split", True), ("start", False)):
            v = f(da, resample_before_rl=before)
            out[f"xclim_{index}_{rule}_{point}"] = pd.Series(v.to_numpy(), index=v.time.dt.year.to_numpy())
    ds = xr.Dataset({"pr": da})
    for index in ("CDD", "CWD"):
        v = icclim.index(index_name=index, in_files=ds, slice_mode="year").to_array().squeeze()
        out[f"icclim_{index.lower()}_{point}"] = pd.Series(v.to_numpy(), index=v.time.dt.year.to_numpy())

df = pd.DataFrame(out)
df.index.name = "year"
df.to_csv(RES / "etccdi_spells.csv", float_format="%.10g", na_rep="NaN", lineterminator="\n")
print(f"xclim {xclim.__version__}, icclim {icclim.__version__}: {df.shape[1]} series, {df.index[0]}-{df.index[-1]}")
