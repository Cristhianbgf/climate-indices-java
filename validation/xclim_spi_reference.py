"""Reference SPI values from xclim for the Java tests.

Reads src/test/resources/basin_ppt.csv and writes src/test/resources/xclim_spi.csv, with one
column per basin, window (1, 3 and 12 months) and fit method. Calibration period 1991-2020,
two-parameter gamma (floc=0), zeros handled as xclim does by default.

The accumulation over the window is done here with exact sums and passed to xclim with window=1.
xclim's own rolling mean (xarray, computed with moving sums) leaves residues of about
1e-15 where all months in the window are dry, so those months are no longer zeros: they enter
the gamma fit instead of the probability of zero. That changes SPI-3 in basin 1391 by up to
2.8. See validation/README.md.

Run from the repository root:  python validation/xclim_spi_reference.py
Tested with xclim 0.61.1, scipy 1.18.0 and xarray 2026.7.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import xclim
from xclim.indices import standardized_precipitation_index

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "src" / "test" / "resources"
WINDOWS = (1, 3, 12)
METHODS = ("APP", "ML")


def window_mean(values: np.ndarray, window: int) -> np.ndarray:
    """Mean of the last `window` months, summed term by term so that dry spells stay at 0."""
    out = np.full(values.shape, np.nan)
    for i in range(window - 1, len(values)):
        out[i] = values[i - window + 1 : i + 1].sum() / window
    return out


ppt = pd.read_csv(RES / "basin_ppt.csv", index_col="month")
time = pd.DatetimeIndex(pd.to_datetime(ppt.index + "-01").to_numpy(), name="time")

out = pd.DataFrame(index=ppt.index)
for col in ppt.columns:
    for window in WINDOWS:
        acc = window_mean(ppt[col].to_numpy(), window)
        pr = xr.DataArray(acc, coords={"time": time}, dims="time", attrs={"units": "mm/month"})
        for method in METHODS:
            spi = standardized_precipitation_index(
                pr, freq=None, window=1, dist="gamma", method=method,
                fitkwargs={"floc": 0}, cal_start="1991-01-01", cal_end="2020-12-31",
            )
            out[f"{col.replace('ppt_', 'spi_')}_w{window}_{method}"] = spi.to_numpy()

out.round(10).to_csv(RES / "xclim_spi.csv", na_rep="NaN")
print(f"xclim {xclim.__version__}: {out.shape[1]} series written to {RES / 'xclim_spi.csv'}")
