package io.github.cristhianbgf.climind;

import io.github.cristhianbgf.climind.dist.FitMethod;
import io.github.cristhianbgf.climind.dist.GammaFit;
import io.github.cristhianbgf.climind.dist.StandardNormal;
import java.time.YearMonth;
import java.util.Arrays;

/**
 * Standardized Precipitation Index (McKee et al., 1993) from monthly precipitation.
 *
 * <p>For each calendar month, a two-parameter gamma distribution is fitted to the non-zero values of
 * the calibration period, and the probability of zero is the share of zeros in that period (the
 * "ecdf" plotting position with the "upper" choice in xclim). Values are then converted to standard
 * normal quantiles and limited to ±8.21, as in xclim.
 */
public final class Spi {

    /** The largest absolute value xclim returns: the normal quantile of the float64 neighbours of 0 and 1. */
    public static final double LIMIT = 8.21;

    private Spi() {}

    /** SPI calibrated on the whole series. */
    public static MonthlySeries compute(MonthlySeries precipitation, int window, FitMethod method) {
        return compute(precipitation, window, method, precipitation.start(), precipitation.end());
    }

    /**
     * @param precipitation monthly totals or means; the units do not matter
     * @param window        accumulation period in months (SPI-1, SPI-3, SPI-12...)
     * @param calibrationStart first month used to fit the distributions
     * @param calibrationEnd   last month used to fit the distributions
     */
    public static MonthlySeries compute(MonthlySeries precipitation, int window, FitMethod method,
                                        YearMonth calibrationStart, YearMonth calibrationEnd) {
        if (calibrationEnd.isBefore(calibrationStart)) {
            throw new IllegalArgumentException("calibration ends before it starts");
        }
        MonthlySeries acc = precipitation.rollingMean(window);

        GammaFit[] fits = new GammaFit[13];
        double[] probZero = new double[13];
        for (int m = 1; m <= 12; m++) {
            double[] cal = acc.calendarMonth(m, calibrationStart, calibrationEnd);
            long valid = Arrays.stream(cal).filter(v -> !Double.isNaN(v)).count();
            long zeros = Arrays.stream(cal).filter(v -> v == 0).count();
            fits[m] = GammaFit.fit(cal, method);
            probZero[m] = valid == 0 ? Double.NaN : (double) zeros / valid;
        }

        double[] spi = new double[acc.size()];
        for (int i = 0; i < spi.length; i++) {
            int m = acc.monthAt(i).getMonthValue();
            double x = acc.get(i);
            double p = x == 0 ? probZero[m] : probZero[m] + (1 - probZero[m]) * fits[m].cdf(x);
            spi[i] = toStandardNormal(p);
        }
        return new MonthlySeries(acc.start(), spi);
    }

    static double toStandardNormal(double p) {
        if (Double.isNaN(p)) {
            return Double.NaN;
        }
        double z = StandardNormal.quantile(p);
        return Math.max(-LIMIT, Math.min(LIMIT, z));
    }
}
