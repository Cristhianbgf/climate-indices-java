# climind

Drought and climate extreme indices in Java, checked against [xclim](https://github.com/Ouranosinc/xclim) and the R package [SPEI](https://cran.r-project.org/package=SPEI). SPI and SPEI are done; the ETCCDI indices are planned.

I could not find a maintained Java library for these indices. People who work in Java (THREDDS, ERDDAP, GeoServer or institutional systems) usually call Python or R to get them. The aim here is to compute them in Java and to document how the results compare with the reference implementations.

Work in progress. See the open issues for what comes next.

## What is there

- **SPI** (McKee et al., 1993) from monthly precipitation, at any accumulation window, with a two-parameter gamma fitted per calendar month on a chosen calibration period. Zeros get their own probability, as in xclim.
- Fits by Thom's approximation or by maximum likelihood.

- **SPEI** (Vicente-Serrano et al., 2010) from monthly precipitation and reference evapotranspiration, with a log-logistic distribution fitted by unbiased probability weighted moments, as in the R package SPEI.

```java
YearMonth start = YearMonth.of(1958, 1);
MonthlySeries pr = new MonthlySeries(start, monthlyPrecipitation);
MonthlySeries pet = new MonthlySeries(start, monthlyEvapotranspiration);
YearMonth calStart = YearMonth.of(1991, 1);
YearMonth calEnd = YearMonth.of(2020, 12);

MonthlySeries spi3 = Spi.compute(pr, 3, FitMethod.MAXIMUM_LIKELIHOOD, calStart, calEnd);
MonthlySeries spei12 = Spei.compute(pr, pet, 12, calStart, calEnd);
```

## Validation

Both indices are checked on three Peruvian basins, 1958-2025. Details in [validation/README.md](validation/README.md).

- SPI-1, SPI-3 and SPI-12 agree with xclim 0.61.1 to within 1e-8 in every month. Getting there showed that xclim's rolling mean turns dry windows into tiny positive values, which shifts SPI-3 in a hyper-arid basin by up to 2.8.
- SPEI-1, SPEI-3 and SPEI-12 agree with the R package SPEI 1.8.1 to within 1e-8, and so do its 108 fitted distributions. Where R returns an infinite SPEI (27 months outside the reference period), this library returns ±8.21.

## Build

Requires Java 21 and Maven.

```bash
mvn test
```

## Licence

MIT. See [LICENSE](LICENSE).
