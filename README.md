# climind

Drought and climate extreme indices in Java, checked against [xclim](https://github.com/Ouranosinc/xclim). The SPI comes first; SPEI and the ETCCDI indices are planned.

I could not find a maintained Java library for these indices. People who work in Java (THREDDS, ERDDAP, GeoServer or institutional systems) usually call Python or R to get them. The aim here is to compute them in Java and to document how the results compare with the reference implementations.

Work in progress. See the open issues for what comes next.

## What is there

- **SPI** (McKee et al., 1993) from monthly precipitation, at any accumulation window, with a two-parameter gamma fitted per calendar month on a chosen calibration period. Zeros get their own probability, as in xclim.
- Fits by Thom's approximation or by maximum likelihood.

```java
MonthlySeries pr = new MonthlySeries(YearMonth.of(1958, 1), monthlyTotals);
MonthlySeries spi3 = Spi.compute(pr, 3, FitMethod.MAXIMUM_LIKELIHOOD,
        YearMonth.of(1991, 1), YearMonth.of(2020, 12));
```

## Validation

SPI-1, SPI-3 and SPI-12 for three Peruvian basins, 1958-2025, agree with xclim 0.61.1 to within 1e-8 in every month. Getting there showed that xclim's rolling mean turns dry windows into tiny positive values, which shifts SPI-3 in a hyper-arid basin by up to 2.8. Details in [validation/README.md](validation/README.md).

## Build

Requires Java 21 and Maven.

```bash
mvn test
```

## Licence

MIT. See [LICENSE](LICENSE).
