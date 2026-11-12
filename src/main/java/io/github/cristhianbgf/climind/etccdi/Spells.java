package io.github.cristhianbgf.climind.etccdi;

import io.github.cristhianbgf.climind.AnnualSeries;
import io.github.cristhianbgf.climind.DailySeries;
import java.time.LocalDate;
import java.util.function.DoublePredicate;

/**
 * Longest dry and wet spells per calendar year (ETCCDI CDD and CWD), from daily precipitation in mm.
 *
 * <p>A dry day has less than 1 mm and a wet day at least 1 mm. The ETCCDI list does not settle how a
 * spell that crosses 31 December is counted, and implementations differ (see {@link YearRule}); the
 * choice can change CWD by months in places whose wet season spans the new year.
 *
 * <p>Only calendar years fully covered by the series get a value, and a year with any missing day is
 * {@code NaN}. A missing day ends a spell.
 */
public final class Spells {

    /** How a spell that crosses 31 December is counted. */
    public enum YearRule {
        /** Cut at the change of year, each part counted in its own year. xclim's default. */
        SPLIT,
        /** Counted whole in the year it starts. icclim's default, and xclim with {@code resample_before_rl=False}. */
        YEAR_OF_START,
        /** Counted whole in the year it ends. */
        YEAR_OF_END
    }

    private Spells() {}

    /** CDD: longest run of days with less than 1 mm. */
    public static AnnualSeries cdd(DailySeries pr, YearRule rule) {
        return longest(pr, v -> v < PrecipitationAmounts.WET_DAY, rule);
    }

    /** CWD: longest run of days with at least 1 mm. */
    public static AnnualSeries cwd(DailySeries pr, YearRule rule) {
        return longest(pr, v -> v >= PrecipitationAmounts.WET_DAY, rule);
    }

    /** CDD with spells cut at the change of year. */
    public static AnnualSeries cdd(DailySeries pr) {
        return cdd(pr, YearRule.SPLIT);
    }

    /** CWD with spells cut at the change of year. */
    public static AnnualSeries cwd(DailySeries pr) {
        return cwd(pr, YearRule.SPLIT);
    }

    /** Longest run of days that meet {@code condition}, per calendar year. */
    public static AnnualSeries longest(DailySeries pr, DoublePredicate condition, YearRule rule) {
        int first = PrecipitationAmounts.firstFullYear(pr);
        int last = PrecipitationAmounts.lastFullYear(pr);
        int n = Math.max(0, last - first + 1);
        double[] out = new double[n];
        for (int y = first; y <= last; y++) {
            out[y - first] = hasMissing(pr.year(y)) ? Double.NaN : 0;
        }

        int i = 0;
        while (i < pr.size()) {
            double v = pr.get(i);
            if (Double.isNaN(v) || !condition.test(v)) {
                i++;
                continue;
            }
            int start = i;
            while (i + 1 < pr.size() && !Double.isNaN(pr.get(i + 1)) && condition.test(pr.get(i + 1))) {
                i++;
            }
            credit(pr, start, i, rule, first, out);
            i++;
        }
        return new AnnualSeries(first, out);
    }

    /** Records the spell from {@code start} to {@code end} (inclusive) in the year or years it counts for. */
    private static void credit(DailySeries pr, int start, int end, YearRule rule, int first, double[] out) {
        switch (rule) {
            case YEAR_OF_START -> record(pr.dateAt(start).getYear(), end - start + 1, first, out);
            case YEAR_OF_END -> record(pr.dateAt(end).getYear(), end - start + 1, first, out);
            case SPLIT -> {
                int from = start;
                while (from <= end) {
                    int year = pr.dateAt(from).getYear();
                    int lastOfYear = pr.indexOf(LocalDate.of(year, 12, 31)); // -1 if the series ends earlier
                    int to = lastOfYear < 0 ? end : Math.min(end, lastOfYear);
                    record(year, to - from + 1, first, out);
                    from = to + 1;
                }
            }
        }
    }

    private static void record(int year, int length, int first, double[] out) {
        int k = year - first;
        if (k >= 0 && k < out.length && !Double.isNaN(out[k])) {
            out[k] = Math.max(out[k], length);
        }
    }

    private static boolean hasMissing(double[] days) {
        for (double v : days) {
            if (Double.isNaN(v)) {
                return true;
            }
        }
        return false;
    }
}
