"""Daily precipitation for the ETCCDI tests: ERA5 at one point in each of the three test basins.

Downloads 1960-2025 daily precipitation sums from the Open-Meteo historical weather API with
models=era5 and writes src/test/resources/era5_daily_pr.csv (mm/day).

The points are the centroids of the hydrographic units, all inside their unit:
  1391 Intercuenca 1391 (Pacific coast)  -4.683, -81.166
  132  Quilca-Vítor-Chili               -16.218, -71.593
  4994 Urubamba                         -12.244, -72.548

Data: ERA5 (Hersbach et al., 2020, doi:10.1002/qj.3803), Copernicus Climate Change Service,
served by Open-Meteo (Zippenfenig, 2023, doi:10.5281/zenodo.7970649) under CC BY 4.0. See the
attribution notes in validation/README.md.

Run from the repository root:  python validation/fetch_era5_daily.py
"""
import json
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

try:
    # use the operating system's certificate store when the Python one is out of date
    import truststore

    truststore.inject_into_ssl()
except ImportError:
    pass

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "src" / "test" / "resources" / "era5_daily_pr.csv"
POINTS = {"1391": (-4.683, -81.166), "132": (-16.218, -71.593), "4994": (-12.244, -72.548)}

series = {}
for unit, (lat, lon) in POINTS.items():
    query = urllib.parse.urlencode({
        "latitude": lat, "longitude": lon, "start_date": "1960-01-01", "end_date": "2025-12-31",
        "daily": "precipitation_sum", "models": "era5", "timezone": "UTC",
    })
    with urllib.request.urlopen("https://archive-api.open-meteo.com/v1/archive?" + query, timeout=120) as r:
        daily = json.load(r)["daily"]
    series[f"pr_{unit}"] = pd.Series(daily["precipitation_sum"], index=daily["time"], dtype=float)

df = pd.DataFrame(series)
df.index.name = "date"
df.round(2).to_csv(OUT, na_rep="NaN", lineterminator="\n")
print(f"{len(df)} days written to {OUT}; missing values: {int(df.isna().sum().sum())}")
