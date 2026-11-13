package io.github.cristhianbgf.climind.etccdi;

import io.github.cristhianbgf.climind.AnnualSeries;
import io.github.cristhianbgf.climind.DailySeries;
import java.util.Arrays;

/**
 * Precipitation on very wet and extremely wet days (ETCCDI R95pTOT and R99pTOT), from daily
 * precipitation in mm.
 *
 * <p>The threshold is a percentile of the wet days (at least 1 mm) of a base period, 1961-1990 in the
 * ETCCDI list; days with missing values are left out of it. A day counts when it is above the threshold,
 * strictly. The base period is used as given; the ETCCDI list asks for a bootstrap only in the
 * temperature indices.
 *
 * <p>Implementations disagree on the unit. The ETCCDI list gives R95pTOT as the total in mm;
 * ECA&amp;D, and so icclim, give it as a percentage of PRCPTOT. Both are here: {@link #r95pTot} in mm
 * and {@link #r95pTotPercent} in percent.
 *
 * <p>Only calendar years fully covered by the series get a value, and a year with any missing day is
 * {@code NaN}.
 */
public final class VeryWetDays {

    /** How a percentile is interpolated between order statistics (Hyndman and Fan, 1996). */
    public enum QuantileMethod {
        /** Type 8, approximately median-unbiased. icclim's default, and the one used by the index methods here. */
        MEDIAN_UNBIASED,
        /** Type 7, the default of numpy and R. */
        LINEAR
    }

    private VeryWetDays() {}

    /** R95pTOT in mm: precipitation on days above the 95th percentile of base-period wet days. */
    public static AnnualSeries r95pTot(DailySeries pr, int baseStart, int baseEnd) {
        return amountAbove(pr, wetDayPercentile(pr, 95, baseStart, baseEnd, QuantileMethod.MEDIAN_UNBIASED));
    }

    /** R99pTOT in mm: precipitation on days above the 99th percentile of base-period wet days. */
    public static AnnualSeries r99pTot(DailySeries pr, int baseStart, int baseEnd) {
        return amountAbove(pr, wetDayPercentile(pr, 99, baseStart, baseEnd, QuantileMethod.MEDIAN_UNBIASED));
    }

    /** R95pTOT as ECA&amp;D defines it: percentage of PRCPTOT from days above the 95th percentile. */
    public static AnnualSeries r95pTotPercent(DailySeries pr, int baseStart, int baseEnd) {
        return percentAbove(pr, wetDayPercentile(pr, 95, baseStart, baseEnd, QuantileMethod.MEDIAN_UNBIASED));
    }

    /** R99pTOT as ECA&amp;D defines it: percentage of PRCPTOT from days above the 99th percentile. */
    public static AnnualSeries r99pTotPercent(DailySeries pr, int baseStart, int baseEnd) {
        return percentAbove(pr, wetDayPercentile(pr, 99, baseStart, baseEnd, QuantileMethod.MEDIAN_UNBIASED));
    }

    /**
     * Percentile of the wet days (at least 1 mm) from 1 January {@code baseStart} to 31 December
     * {@code baseEnd}; {@code NaN} if that period has no wet day.
     *
     * @param percent between 0 and 100
     * @throws IllegalArgumentException if the series does not cover the base period
     */
    public static double wetDayPercentile(DailySeries pr, double percent, int baseStart, int baseEnd,
            QuantileMethod method) {
        if (!(percent >= 0 && percent <= 100)) {
            throw new IllegalArgumentException("percent must be between 0 and 100, got " + percent);
        }
        if (baseStart > baseEnd || !pr.coversYear(baseStart) || !pr.coversYear(baseEnd)) {
            throw new IllegalArgumentException(
                    "the series does not cover the base period " + baseStart + "-" + baseEnd);
        }
        double[] wet = new double[0];
        int n = 0;
        for (int y = baseStart; y <= baseEnd; y++) {
            for (double v : pr.year(y)) {
                if (v >= PrecipitationAmounts.WET_DAY) { // false for NaN
                    if (n == wet.length) {
                        wet = Arrays.copyOf(wet, Math.max(16, 2 * n));
                    }
                    wet[n++] = v;
                }
            }
        }
        return quantile(Arrays.copyOf(wet, n), percent / 100, method);
    }

    /**
     * Precipitation in mm on the days above {@code threshold} that are also wet days, per calendar
     * year. Every year is {@code NaN} if the threshold is.
     */
    public static AnnualSeries amountAbove(DailySeries pr, double threshold) {
        return PrecipitationAmounts.annual(pr, days -> Double.isNaN(threshold) ? Double.NaN : over(days, threshold));
    }

    /**
     * Percentage of PRCPTOT that falls on the days above {@code threshold}, per calendar year;
     * {@code NaN} in a year without wet days.
     */
    public static AnnualSeries percentAbove(DailySeries pr, double threshold) {
        return PrecipitationAmounts.annual(pr, days -> {
            double total = 0;
            for (double v : days) {
                if (v >= PrecipitationAmounts.WET_DAY) {
                    total += v;
                }
            }
            return total == 0 || Double.isNaN(threshold) ? Double.NaN : 100 * over(days, threshold) / total;
        });
    }

    private static double over(double[] days, double threshold) {
        double sum = 0;
        for (double v : days) {
            if (v > threshold && v >= PrecipitationAmounts.WET_DAY) {
                sum += v;
            }
        }
        return sum;
    }

    /**
     * Quantile {@code q} of {@code values}, computed as numpy does for the two methods, so that the
     * results match to the last bits; {@code NaN} if there are no values.
     */
    static double quantile(double[] values, double q, QuantileMethod method) {
        int n = values.length;
        if (n == 0) {
            return Double.NaN;
        }
        double[] x = values.clone();
        Arrays.sort(x);
        double alpha = method == QuantileMethod.LINEAR ? 1 : 1.0 / 3;
        double beta = alpha;
        double index = n * q + (alpha + q * (1 - alpha - beta)) - 1;
        index = Math.min(Math.max(index, 0), n - 1);
        int below = (int) Math.floor(index);
        int above = Math.min(below + 1, n - 1);
        double t = index - below;
        double a = x[below];
        double b = x[above];
        double diff = b - a;
        return t >= 0.5 ? b - diff * (1 - t) : a + diff * t;
    }
}
