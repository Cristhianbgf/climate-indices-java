package io.github.cristhianbgf.climind;

import io.github.cristhianbgf.climind.dist.LogLogisticFit;
import java.time.YearMonth;

/**
 * Standardized Precipitation-Evapotranspiration Index (Vicente-Serrano et al., 2010) from monthly
 * precipitation and reference evapotranspiration.
 *
 * <p>The climatic water balance (precipitation minus evapotranspiration) is accumulated over the
 * window and, for each calendar month, a log-logistic distribution is fitted by unbiased probability
 * weighted moments on the calibration period, as in the R package SPEI. Values beyond the support of
 * the fitted distribution, which R returns as infinite, are limited to ±8.21 as in {@link Spi}.
 */
public final class Spei {

    private Spei() {}

    /** SPEI calibrated on the whole series. */
    public static MonthlySeries compute(MonthlySeries precipitation, MonthlySeries evapotranspiration, int window) {
        return compute(precipitation, evapotranspiration, window, precipitation.start(), precipitation.end());
    }

    /**
     * @param precipitation      monthly precipitation
     * @param evapotranspiration monthly reference evapotranspiration, same units and months
     * @param window             accumulation period in months
     * @param calibrationStart   first month used to fit the distributions
     * @param calibrationEnd     last month used to fit the distributions
     */
    public static MonthlySeries compute(MonthlySeries precipitation, MonthlySeries evapotranspiration, int window,
                                        YearMonth calibrationStart, YearMonth calibrationEnd) {
        if (!precipitation.start().equals(evapotranspiration.start())
                || precipitation.size() != evapotranspiration.size()) {
            throw new IllegalArgumentException("precipitation and evapotranspiration must cover the same months");
        }
        if (calibrationEnd.isBefore(calibrationStart)) {
            throw new IllegalArgumentException("calibration ends before it starts");
        }
        double[] balance = new double[precipitation.size()];
        for (int i = 0; i < balance.length; i++) {
            balance[i] = precipitation.get(i) - evapotranspiration.get(i);
        }
        MonthlySeries acc = new MonthlySeries(precipitation.start(), balance).rollingMean(window);

        LogLogisticFit[] fits = new LogLogisticFit[13];
        for (int m = 1; m <= 12; m++) {
            fits[m] = LogLogisticFit.fitPwm(acc.calendarMonth(m, calibrationStart, calibrationEnd));
        }
        double[] spei = new double[acc.size()];
        for (int i = 0; i < spei.length; i++) {
            spei[i] = Spi.toStandardNormal(fits[acc.monthAt(i).getMonthValue()].cdf(acc.get(i)));
        }
        return new MonthlySeries(acc.start(), spei);
    }
}
