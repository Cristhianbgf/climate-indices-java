package io.github.cristhianbgf.climind;

import io.github.cristhianbgf.climind.dist.FitMethod;
import io.github.cristhianbgf.climind.dist.LogLogisticFit;
import java.time.YearMonth;

/**
 * Standardized Precipitation-Evapotranspiration Index (Vicente-Serrano et al., 2010) from monthly
 * precipitation and reference evapotranspiration.
 *
 * <p>The climatic water balance (precipitation minus evapotranspiration) is accumulated over the
 * window and, for each calendar month, a log-logistic distribution is fitted on the calibration
 * period, by default with unbiased probability weighted moments as in the R package SPEI, or by
 * maximum likelihood. Values beyond the support of the fitted distribution, which R returns as
 * infinite, are limited to ±8.21 as in {@link Spi}.
 */
public final class Spei {

    private Spei() {}

    /** SPEI with the PWM fit, calibrated on the whole series. */
    public static MonthlySeries compute(MonthlySeries precipitation, MonthlySeries evapotranspiration, int window) {
        return compute(precipitation, evapotranspiration, window, precipitation.start(), precipitation.end());
    }

    /** SPEI with the PWM fit. */
    public static MonthlySeries compute(MonthlySeries precipitation, MonthlySeries evapotranspiration, int window,
                                        YearMonth calibrationStart, YearMonth calibrationEnd) {
        return compute(precipitation, evapotranspiration, window, FitMethod.PROBABILITY_WEIGHTED_MOMENTS,
                calibrationStart, calibrationEnd);
    }

    /**
     * @param precipitation      monthly precipitation
     * @param evapotranspiration monthly reference evapotranspiration, same units and months
     * @param window             accumulation period in months
     * @param method             {@code PROBABILITY_WEIGHTED_MOMENTS} or {@code MAXIMUM_LIKELIHOOD}
     * @param calibrationStart   first month used to fit the distributions
     * @param calibrationEnd     last month used to fit the distributions
     */
    public static MonthlySeries compute(MonthlySeries precipitation, MonthlySeries evapotranspiration, int window,
                                        FitMethod method, YearMonth calibrationStart, YearMonth calibrationEnd) {
        if (method == FitMethod.APPROXIMATE) {
            throw new IllegalArgumentException("the SPEI supports PROBABILITY_WEIGHTED_MOMENTS and MAXIMUM_LIKELIHOOD");
        }
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
            double[] cal = acc.calendarMonth(m, calibrationStart, calibrationEnd);
            fits[m] = method == FitMethod.MAXIMUM_LIKELIHOOD
                    ? LogLogisticFit.fitMaximumLikelihood(cal)
                    : LogLogisticFit.fitPwm(cal);
        }
        double[] spei = new double[acc.size()];
        for (int i = 0; i < spei.length; i++) {
            spei[i] = Spi.toStandardNormal(fits[acc.monthAt(i).getMonthValue()].cdf(acc.get(i)));
        }
        return new MonthlySeries(acc.start(), spei);
    }
}
