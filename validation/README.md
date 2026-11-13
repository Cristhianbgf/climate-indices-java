# Validation

## Data

`src/test/resources/basin_ppt.csv` and `basin_pet.csv` hold monthly precipitation and reference evapotranspiration (mm) from January 1958 to December 2025 for three hydrographic units of Peru. They are averages of TerraClimate v1.1 (Abatzoglou et al., licence CC0 1.0) over each unit, computed with the pipeline of [peru-basin-drought](https://github.com/Cristhianbgf/peru-basin-drought), which publishes SPEI per unit but not these series:

| Unit | Name | Why it is here |
|---|---|---|
| 1391 | Intercuenca 1391, Pacific coast | Hyper-arid: 108 of 816 months without rain, so the probability of zero matters |
| 132 | Quilca-Vítor-Chili | Semi-arid Andean basin with a long dry season |
| 4994 | Urubamba | Humid Andean and Amazonian basin |

`src/test/resources/era5_daily_pr.csv` holds daily precipitation (mm) from 1 January 1960 to 31 December 2025 at the centroid of each of the same three units, made by `fetch_era5_daily.py`. It comes from ERA5 (Hersbach et al., 2020, [doi:10.1002/qj.3803](https://doi.org/10.1002/qj.3803)), served by the Open-Meteo historical weather API (Zippenfenig, 2023, [doi:10.5281/zenodo.7970649](https://doi.org/10.5281/zenodo.7970649)) under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Changes: Open-Meteo extracts the series at the nearest grid point, and the values are rounded to 0.01 mm. Contains modified Copernicus Climate Change Service information (2026); neither the European Commission nor ECMWF is responsible for any use that may be made of the Copernicus information or data it contains. Weather data by [Open-Meteo.com](https://open-meteo.com/).

## SPI against xclim

`xclim_spi_reference.py` computes SPI-1, SPI-3 and SPI-12 with xclim 0.61.1 (two-parameter gamma, approximate and maximum likelihood fits, calibration 1991-2020) and writes `src/test/resources/xclim_spi.csv`. `SpiXclimTest` checks every value with a tolerance of 1e-8.

The 18 series (3 basins, 3 windows, 2 fit methods) agree with xclim to within 1e-8 in all 14 610 months with a value, including the extreme tail (SPI-3 of -6.95 in November 1972 at Urubamba). The other 78 are the first months of each window, missing in both.

Two things had to be handled to get there:

- **Tail precision.** Commons Math inverts the normal CDF through `erfInv(2p - 1)`, which loses precision for probabilities near 1e-12. `StandardNormal` uses Acklam's approximation refined with Halley's method instead.
- **Digamma.** Commons Math's digamma is accurate only to about 1e-8 near 1, which moved the maximum likelihood shape by the same amount. `Polygamma` uses the recurrence and the asymptotic series.

### A difference in xclim

xclim accumulates precipitation with xarray's rolling mean, which uses moving sums. Where every month in the window is dry, the result is a residue of about 1e-15 instead of 0. Those windows are then fitted as rain instead of counted in the probability of zero. For SPI-3 in unit 1391 this changes 68 of the 816 months by more than 0.01, by up to 2.8 with the approximate fit and 2.5 with maximum likelihood. For example, September 1958 is 0.15 in xclim and -0.46 here with the approximate fit. The other two basins have no dry windows and are not affected.

This library sums the window term by term, so dry windows stay at zero. The reference script does the same and passes the accumulated series to xclim with `window=1`, so that the comparison tests the fit and the transformation, not the residue.

## SPEI against the R package SPEI

`r_spei_reference.R` computes SPEI-1, SPEI-3 and SPEI-12 with the R package SPEI 1.8.1, which fits a log-logistic distribution (Hosking's generalized logistic) by unbiased probability weighted moments, with 1991-2020 as reference period. It writes `src/test/resources/r_spei.csv` and the fitted parameters to `r_glo_params.csv`.

- **Parameters.** All 108 fits (3 basins, 3 windows, 12 months) agree with R: the shape to within 1e-12 and the location and scale to within 1e-9 in relative terms (`LogLogisticFitTest`). R accumulates sums and this library means; the fit scales with the data, so location and scale are compared after multiplying by the window.
- **Index.** The 9 series agree to within 1e-8 in all 7 278 months with a finite value (`SpeiRTest`).
- **Out of support.** In 27 months, all outside the reference period, the water balance falls beyond the bound of the fitted distribution (it has one whenever the shape is not zero), and R returns an infinite SPEI. This library returns ±8.21 there, the same limit as for the SPI. The other 39 values are the first months of each window, missing in both.

## SPEI by maximum likelihood

`glo_ml_reference.py` fits the same 108 samples by maximum likelihood twice: in the generalized logistic form, with code written independently of the Java one (Nelder-Mead from the PWM fit, restarted until the log-likelihood stops improving), and as xclim 0.61.1 does it for `dist="fisk"`, `method="ML"` (xclim's starting values from `_fit_start`, then `scipy.stats.fisk.fit`). It writes `src/test/resources/glo_ml_params.csv` and `xclim_fisk_ml.csv`.

- **Against the independent fit.** All 108 log-likelihoods agree to within 1e-7, and the parameters to within 1e-4 of the scale (`LogLogisticMaximumLikelihoodTest`). The optimum is flat, so the log-likelihood is the meaningful comparison.
- **Against xclim.** The fit here is never worse, and it is better in 44 of the 108 samples, by up to 1.8 in log-likelihood.

### Why xclim's fit falls short

xclim fits scipy's `fisk`, the log-logistic with shape `c > 0`, which in the generalized logistic form is `kappa = -1/c < 0`: a distribution bounded below. In 44 of the 108 samples the best fit has `kappa > 0` instead, bounded above, because the water balance of that month is skewed to the left. `fisk` cannot reach it; the closest it gets is the logistic, at `c` going to infinity. The optimizer heads there and stops wherever its tolerance says, with shapes between 1e5 and 1e9. These are exactly the 44 samples where the fit here does better. In one more sample (Urubamba, SPEI-3, March) xclim also returns a shape above 1000, but there the optimum is practically the logistic (`kappa` = -0.0008) and both fits reach the same log-likelihood. In the generalized logistic form both signs of `kappa` are ordinary values, so the fit converges.

### Effect on the index, and why PWM stays the default

In the 3 240 monthly values of the calibration period, the SPEI from the maximum likelihood fit here differs:

- from xclim's by up to 0.72, and by more than 0.1 in 191 values;
- from the PWM fit by up to 0.81, and by more than 0.1 in 290 values.

Outside the calibration period, values that fall outside the support of the fitted distribution get ±8.21. That happens in 98 months with the maximum likelihood fit (80 dry months below the lower bound of a fit with `kappa < 0`, and 18 wet months at the upper bound of one with `kappa > 0`: 17 beyond it and one so close that the probability rounds to 1), in 80 with xclim's (all dry) and in 27 with PWM (23 dry and 4 wet). The maximum likelihood fit follows the calibration sample more closely, and its bounds sit closer to the sample; PWM places them further out and extrapolates better. That is why `Spei.compute` uses PWM unless asked otherwise.

## ETCCDI precipitation indices against xclim

The daily series of ERA5 (see Data above) cover 24 107 days with none missing. They are point series, not basin averages. The point in Intercuenca 1391 has at least 1 mm on 8 % of the days, the one in Quilca-Vítor-Chili on 43 % and the one in Urubamba on 86 %.

`xclim_etccdi_reference.py` computes the annual indices with xclim 0.61.1 and writes `src/test/resources/xclim_etccdi.csv`:

| Index | Definition | xclim function |
|---|---|---|
| PRCPTOT | total precipitation on wet days (at least 1 mm) | `prcptot(thresh="1 mm/d")` |
| Rx1day | highest daily precipitation | `max_1day_precipitation_amount` |
| Rx5day | highest precipitation over five consecutive days | `max_n_day_precipitation_amount(window=5)` |
| SDII | mean precipitation on wet days | `daily_pr_intensity` |
| R10mm, R20mm | days with at least 10 and 20 mm | `wetdays(thresh=...)` |

All 18 series (6 indices, 3 points, 66 years) agree with xclim to within 1e-9 relative (`PrecipitationAmountsTest`).

A five-day window belongs to the year of its last day, so the one ending on 1 to 4 January takes days from December, as in xclim.

## CDD and CWD

`etccdi_spells_reference.py` computes the longest dry spell (CDD, days with less than 1 mm) and wet spell (CWD, days with at least 1 mm) per calendar year on the same ERA5 series, in three ways, and writes `src/test/resources/etccdi_spells.csv`:

- with numpy, written from the ETCCDI definitions and independently of xclim, under each of three rules for a spell that crosses 31 December: cut at the change of year, counted whole in the year it starts, or counted whole in the year it ends;
- with xclim 1.0.0, with its default (`resample_before_rl=True`, which cuts the spell) and with `resample_before_rl=False` (the whole spell goes to the year it starts);
- with icclim 7.3.7 and its ECA&D defaults, which count the whole spell in the year it starts.

icclim runs on xclim, so the numpy version is the independent check. All three agree wherever they follow the same rule: `Spells` matches the numpy version under all three rules, xclim under its two, and icclim, in every year at the three points (`SpellsTest`).

xclim 1.0 counts a wet day for CWD at 1 mm or more, as the ETCCDI list does. xclim 0.61, used above for the amount indices, counted only more than 1 mm, which differs on days of exactly 1.00 mm.

### The rule for the new year matters more than the definition

| Point | Index | Years where the rules disagree | Largest difference (year of start vs year of end) |
|---|---|---|---|
| 1391, Pacific coast | CDD | 12 of 66 | 247 days |
| Quilca-Vítor-Chili | CWD | 25 of 66 | 108 days |
| Urubamba | CWD | 65 of 66 | 168 days |

Where the dry or the wet season spans the new year, the rule changes the index in most years and by months. `Spells` takes the rule as an argument and cuts at the change of year unless told otherwise, as xclim does by default; the choice should be stated with any result.

ERA5 gives at least 1 mm on 86 % of the days at the Urubamba point, so its wet spells reach 233 days. These series test the code, not the climate of the basin; station data would likely give shorter spells.

## R95pTOT and R99pTOT

The threshold is the 95th or 99th percentile of the wet days (at least 1 mm) of 1961-1990, and the index adds up the precipitation of the days above it, strictly. The ETCCDI list gives the result as a total in mm; ECA&D, and so icclim, as a percentage of PRCPTOT. `VeryWetDays` gives both.

`etccdi_percentiles_reference.py` computes the thresholds and the annual indices on the same ERA5 series and writes `src/test/resources/etccdi_percentile_thresholds.csv` and `etccdi_percentiles.csv`:

- with numpy, written from the definitions, with the percentile by two methods: type 8 of Hyndman and Fan (median-unbiased, icclim's default) and type 7 (linear, the default of numpy and R);
- with icclim 7.3.7: its ECA&D R95pTOT and R99pTOT for the percentage, and its generic sum with the same threshold for the total in mm;
- with xclim 1.0.0, `fraction_over_precip_thresh`, given the type 8 threshold computed in numpy.

`VeryWetDays` matches the numpy version under both methods and icclim in every year at the three points, to within 1e-9 relative, and the thresholds to within 1e-12 (`VeryWetDaysTest`). The base period holds 759 wet days at the point in Intercuenca 1391, 4 925 in Quilca-Vítor-Chili and 9 499 in Urubamba.

climdex.pcic, the R implementation of the ETCCDI indices, could not be used: it was archived from CRAN in 2023, so there is no binary for R 4.6, and building it needs Rtools. icclim runs on xclim, so the numpy version is the independent check here.

### Where xclim differs

xclim applies one comparison to both tests: whether a day is wet and whether it is above the percentile. The ETCCDI definition needs at least 1 mm for the first and more than the percentile for the second, so neither setting reproduces it.

| Point | Index | Years that differ with `">"` (largest) | Years that differ with `">="` (largest) |
|---|---|---|---|
| 1391 | R95pTOT | 20 (1.09 points) | 0 |
| 1391 | R99pTOT | 5 (0.67) | 0 |
| Quilca-Vítor-Chili | R95pTOT | 62 (0.16) | 19 (2.83) |
| Quilca-Vítor-Chili | R99pTOT | 51 (0.08) | 5 (4.30) |
| Urubamba | R95pTOT | 66 (0.07) | 12 (1.78) |
| Urubamba | R99pTOT | 60 (0.03) | 0 |

With `">"`, days of exactly 1.0 mm drop out of PRCPTOT; the series has 112, 204 and 274 of them. With `">="`, the days equal to the percentile count as extreme. Open-Meteo gives ERA5 in steps of 0.1 mm, so the percentile can land on a recorded value: 12.8 and 17.5 mm in Quilca-Vítor-Chili, 25.8 mm in Urubamba. Where it does not, as at 1391, `">="` gives the same result as the definition.

### The percentile method

Type 8 and type 7 give the same 95th percentile at Quilca-Vítor-Chili and Urubamba, and thresholds 0.12 mm apart at 1391 that change no year. The 99th percentile sits among few days and the methods part: at 1391 the threshold drops from 51.2 to 48.3 mm with type 7, which moves one year by 50.5 mm (2.7 points); at Quilca-Vítor-Chili five years move, by up to 35.0 mm (4.3 points); at Urubamba one year moves by 48.6 mm (1.5 points). `wetDayPercentile` takes the method as an argument; `r95pTot` and the other index methods use type 8.
