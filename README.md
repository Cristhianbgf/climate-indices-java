# climind

[![tests](https://github.com/Cristhianbgf/climate-indices-java/actions/workflows/tests.yml/badge.svg)](https://github.com/Cristhianbgf/climate-indices-java/actions/workflows/tests.yml)

Drought and climate extreme indices in Java, checked against [xclim](https://github.com/Ouranosinc/xclim), [icclim](https://github.com/cerfacs-globc/icclim) and the R package [SPEI](https://cran.r-project.org/package=SPEI). SPI, SPEI and the ETCCDI precipitation indices of amount, spell length and very wet days are done.

I could not find a maintained Java library for these indices. People who work in Java (THREDDS, ERDDAP, GeoServer or institutional systems) usually call Python or R to get them. The aim here is to compute them in Java and to document how the results compare with the reference implementations.

Work in progress. See the open issues for what comes next.

## What is there

- **SPI** (McKee et al., 1993) from monthly precipitation, at any accumulation window, with a two-parameter gamma fitted per calendar month on a chosen calibration period. Zeros get their own probability, as in xclim.
- Fits by Thom's approximation, by maximum likelihood or by L-moments (probability weighted moments), as in the R package SPEI.
- **SPEI** (Vicente-Serrano et al., 2010) from monthly precipitation and reference evapotranspiration, with a log-logistic distribution fitted by unbiased probability weighted moments, as in the R package SPEI, or by maximum likelihood.
- **ETCCDI precipitation amounts** (Zhang et al., 2011) from daily precipitation: PRCPTOT, Rx1day, Rx5day, SDII, R10mm and R20mm, per calendar year.
- **ETCCDI spells**: CDD and CWD, with a choice of how a spell that crosses 31 December is counted.
- **ETCCDI very wet days**: R95pTOT and R99pTOT from a percentile of base-period wet days, as a total in mm (the ETCCDI list) or as a percentage of PRCPTOT (ECA&D and icclim).

```java
YearMonth start = YearMonth.of(1958, 1);
MonthlySeries pr = new MonthlySeries(start, monthlyPrecipitation);
MonthlySeries pet = new MonthlySeries(start, monthlyEvapotranspiration);
YearMonth calStart = YearMonth.of(1991, 1);
YearMonth calEnd = YearMonth.of(2020, 12);

MonthlySeries spi3 = Spi.compute(pr, 3, FitMethod.MAXIMUM_LIKELIHOOD, calStart, calEnd);
MonthlySeries spei12 = Spei.compute(pr, pet, 12, calStart, calEnd);

DailySeries daily = new DailySeries(LocalDate.of(1960, 1, 1), dailyPrecipitation);
AnnualSeries rx5day = PrecipitationAmounts.rx5day(daily);
AnnualSeries r95pTot = VeryWetDays.r95pTot(daily, 1961, 1990);
```

## Validation

All of it is checked on Peruvian data. Details in [validation/README.md](validation/README.md).

- SPI-1, SPI-3 and SPI-12 agree with xclim 0.61.1 to within 1e-8 in every month. Getting there showed that xclim's rolling mean turns dry windows into tiny positive values, which shifts SPI-3 in a hyper-arid basin by up to 2.8.
- SPEI-1, SPEI-3 and SPEI-12 agree with the R package SPEI 1.8.1 to within 1e-8, and so do its 108 fitted distributions. Where R returns an infinite SPEI (27 months outside the reference period), this library returns ±8.21.
- PRCPTOT, Rx1day, Rx5day, SDII, R10mm and R20mm agree with xclim 0.61.1 to within 1e-9 relative on ERA5 daily precipitation at three points, 1960-2025.
- CDD and CWD agree in every year with an independent implementation under three rules for spells that cross the new year, and with xclim 1.0.0 and icclim 7.3.7 under the rules each of them uses. At the Urubamba point the rule changes CWD in 65 of 66 years, by up to 168 days.
- R95pTOT and R99pTOT agree in every year with an independent implementation and with icclim 7.3.7, both in mm and as a percentage. xclim 1.0.0 differs, because it uses one comparison for whether a day is wet and whether it is above the percentile; at the three points that changes the percentage by up to 4.3 points.
- SPI by L-moments agrees with the R package SPEI 1.8.1 to within 1e-8 in every calendar month whose reference period has no month without rain, and so do its 108 fitted distributions. SPEI 1.8.1 leaves the probability of zero out, so it returns -Inf for every month without rain and lowers the other values of those calendar months; this library adds it, as xclim does.
- Maximum likelihood and L-moments give SPI and SPEI values up to 0.8 apart in the calibration period, and the gap is largest in the tails and in the hyper-arid basin, where SPI-12 crosses -1 under one fit and not the other in 29 months. The comparison by basin, window and month is in [validation/fit_comparison.csv](validation/fit_comparison.csv).
- The maximum likelihood fit is never worse than xclim's and better in 44 of the 108 samples. In those 44, the best fit is bounded above, a shape xclim's log-logistic cannot take, so xclim's optimizer runs off towards infinite parameters.

## Build

Requires Java 21 and Maven.

```bash
mvn test
```

## Licence

The code is under the MIT licence; see [LICENSE](LICENSE).

The test data are not. `src/test/resources/basin_*.csv` derive from TerraClimate (CC0 1.0). `src/test/resources/era5_daily_pr.csv`, and the reference values computed from it, contain modified Copernicus Climate Change Service information (2026), obtained through [Open-Meteo.com](https://open-meteo.com/) under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Neither the European Commission nor ECMWF is responsible for any use that may be made of the Copernicus information or data it contains. Details in [validation/README.md](validation/README.md).
