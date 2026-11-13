package io.github.cristhianbgf.climind.etccdi;

import io.github.cristhianbgf.climind.AnnualSeries;
import io.github.cristhianbgf.climind.DailySeries;
import java.time.LocalDate;
import java.util.function.ToDoubleFunction;

/**
 * Annual ETCCDI indices of precipitation amount and intensity, from daily precipitation in mm.
 *
 * <p>Definitions follow the ETCCDI list (Zhang et al., 2011): a wet day has at least 1 mm. Only
 * calendar years fully covered by the series get a value, and a year with any missing day is
 * {@code NaN}, as xclim's indicators do by default. The thresholds are inclusive ({@code >=}).
 */
public final class PrecipitationAmounts {

    /** Threshold of a wet day, in mm. */
    public static final double WET_DAY = 1.0;

    private PrecipitationAmounts() {}

    /** PRCPTOT: total precipitation on wet days. */
    public static AnnualSeries prcptot(DailySeries pr) {
        return annual(pr, days -> {
            double sum = 0;
            for (double v : days) {
                if (v >= WET_DAY) {
                    sum += v;
                }
            }
            return sum;
        });
    }

    /** Rx1day: highest daily precipitation. */
    public static AnnualSeries rx1day(DailySeries pr) {
        return annual(pr, days -> {
            double max = Double.NEGATIVE_INFINITY;
            for (double v : days) {
                max = Math.max(max, v);
            }
            return max;
        });
    }

    /**
     * Rx5day: highest precipitation over five consecutive days. A window belongs to the year of its
     * last day, so one ending on 1 to 4 January takes days from December, as in xclim; windows that
     * would start before the series are skipped.
     */
    public static AnnualSeries rx5day(DailySeries pr) {
        return rxNday(pr, 5);
    }

    /** Highest precipitation over {@code n} consecutive days, with the same convention as {@link #rx5day}. */
    public static AnnualSeries rxNday(DailySeries pr, int n) {
        if (n < 1) {
            throw new IllegalArgumentException("n must be at least 1, got " + n);
        }
        int first = firstFullYear(pr);
        int last = lastFullYear(pr);
        double[] out = new double[Math.max(0, last - first + 1)];
        for (int y = first; y <= last; y++) {
            int from = pr.indexOf(LocalDate.of(y, 1, 1));
            int to = pr.indexOf(LocalDate.of(y, 12, 31));
            if (hasMissing(pr.year(y))) {
                out[y - first] = Double.NaN;
                continue;
            }
            double max = Double.NEGATIVE_INFINITY;
            for (int end = Math.max(from, n - 1); end <= to; end++) {
                double sum = 0;
                for (int k = end - n + 1; k <= end; k++) {
                    sum += pr.get(k);
                }
                if (!Double.isNaN(sum)) { // a window reaching a missing December day is ignored
                    max = Math.max(max, sum);
                }
            }
            out[y - first] = max == Double.NEGATIVE_INFINITY ? Double.NaN : max;
        }
        return new AnnualSeries(first, out);
    }

    /** SDII: mean precipitation on wet days; {@code NaN} in a year without wet days. */
    public static AnnualSeries sdii(DailySeries pr) {
        return annual(pr, days -> {
            double sum = 0;
            int wet = 0;
            for (double v : days) {
                if (v >= WET_DAY) {
                    sum += v;
                    wet++;
                }
            }
            return wet == 0 ? Double.NaN : sum / wet;
        });
    }

    /** R10mm: days with at least 10 mm. */
    public static AnnualSeries r10mm(DailySeries pr) {
        return daysAtLeast(pr, 10);
    }

    /** R20mm: days with at least 20 mm. */
    public static AnnualSeries r20mm(DailySeries pr) {
        return daysAtLeast(pr, 20);
    }

    /** Rnnmm: days with at least {@code threshold} mm. */
    public static AnnualSeries daysAtLeast(DailySeries pr, double threshold) {
        return annual(pr, days -> {
            int count = 0;
            for (double v : days) {
                if (v >= threshold) {
                    count++;
                }
            }
            return count;
        });
    }

    static AnnualSeries annual(DailySeries pr, ToDoubleFunction<double[]> f) {
        int first = firstFullYear(pr);
        int last = lastFullYear(pr);
        double[] out = new double[Math.max(0, last - first + 1)];
        for (int y = first; y <= last; y++) {
            double[] days = pr.year(y);
            out[y - first] = hasMissing(days) ? Double.NaN : f.applyAsDouble(days);
        }
        return new AnnualSeries(first, out);
    }

    private static boolean hasMissing(double[] days) {
        for (double v : days) {
            if (Double.isNaN(v)) {
                return true;
            }
        }
        return false;
    }

    static int firstFullYear(DailySeries pr) {
        return pr.coversYear(pr.firstYear()) ? pr.firstYear() : pr.firstYear() + 1;
    }

    static int lastFullYear(DailySeries pr) {
        return pr.coversYear(pr.lastYear()) ? pr.lastYear() : pr.lastYear() - 1;
    }
}
