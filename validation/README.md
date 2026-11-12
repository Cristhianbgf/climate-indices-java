# Validation against xclim

## Data

`src/test/resources/basin_ppt.csv` holds monthly precipitation (mm) from January 1958 to December 2025 for three hydrographic units of Peru. They are averages of TerraClimate v1.1 (Abatzoglou et al., licence CC0 1.0) over each unit, computed with the pipeline of [peru-basin-drought](https://github.com/Cristhianbgf/peru-basin-drought), which publishes SPEI per unit but not these precipitation series:

| Column | Unit | Why it is here |
|---|---|---|
| `ppt_1391` | Intercuenca 1391, Pacific coast | Hyper-arid: 108 of 816 months without rain, so the probability of zero matters |
| `ppt_132` | Quilca-Vítor-Chili | Semi-arid Andean basin with a long dry season |
| `ppt_4994` | Urubamba | Humid Andean and Amazonian basin |

## Reference values

`xclim_spi_reference.py` computes SPI-1, SPI-3 and SPI-12 with xclim 0.61.1 (two-parameter gamma, approximate and maximum likelihood fits, calibration 1991-2020) and writes `src/test/resources/xclim_spi.csv`. `SpiXclimTest` checks every value with a tolerance of 1e-8.

## Result

The 18 series (3 basins, 3 windows, 2 fit methods) agree with xclim to within 1e-8 in all 14 610 months with a value, including the extreme tail (SPI-3 of -6.95 in November 1972 at Urubamba). The other 78 are the first months of each window, missing in both.

Two things had to be handled to get there:

- **Tail precision.** Commons Math inverts the normal CDF through `erfInv(2p - 1)`, which loses precision for probabilities near 1e-12. `StandardNormal` uses Acklam's approximation refined with Halley's method instead.
- **Digamma.** Commons Math's digamma is accurate only to about 1e-8 near 1, which moved the maximum likelihood shape by the same amount. `Polygamma` uses the recurrence and the asymptotic series.

## A difference in xclim

xclim accumulates precipitation with xarray's rolling mean, which uses moving sums. Where every month in the window is dry, the result is a residue of about 1e-15 instead of 0. Those windows are then fitted as rain instead of counted in the probability of zero. For SPI-3 in unit 1391 this changes 68 of the 816 months by more than 0.01, by up to 2.8 with the approximate fit and 2.5 with maximum likelihood. For example, September 1958 is 0.15 in xclim and -0.46 here with the approximate fit. The other two basins have no dry windows and are not affected.

This library sums the window term by term, so dry windows stay at zero. The reference script does the same and passes the accumulated series to xclim with `window=1`, so that the comparison tests the fit and the transformation, not the residue.
