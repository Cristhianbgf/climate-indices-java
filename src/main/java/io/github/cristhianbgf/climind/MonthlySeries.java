package io.github.cristhianbgf.climind;

import java.time.YearMonth;
import java.time.temporal.ChronoUnit;
import java.util.Arrays;
import java.util.Objects;
import java.util.stream.IntStream;

/**
 * A regular monthly time series. Missing months are stored as {@code NaN}.
 *
 * <p>Instances are immutable: every operation returns a new series.
 */
public final class MonthlySeries {

    private final YearMonth start;
    private final double[] values;

    public MonthlySeries(YearMonth start, double[] values) {
        this.start = Objects.requireNonNull(start, "start");
        this.values = values.clone();
    }

    public YearMonth start() {
        return start;
    }

    public YearMonth end() {
        return start.plusMonths(values.length - 1L);
    }

    public int size() {
        return values.length;
    }

    public double get(int i) {
        return values[i];
    }

    public YearMonth monthAt(int i) {
        return start.plusMonths(i);
    }

    /** Position of a month in the series, or -1 if it falls outside. */
    public int indexOf(YearMonth month) {
        long i = start.until(month, ChronoUnit.MONTHS);
        return (i < 0 || i >= values.length) ? -1 : (int) i;
    }

    public double[] toArray() {
        return values.clone();
    }

    /**
     * Mean over the current and the previous {@code window - 1} months. The first
     * {@code window - 1} months and any window with a missing month are {@code NaN}.
     */
    public MonthlySeries rollingMean(int window) {
        if (window < 1) {
            throw new IllegalArgumentException("window must be at least 1, got " + window);
        }
        double[] out = new double[values.length];
        Arrays.fill(out, Double.NaN);
        for (int i = window - 1; i < values.length; i++) {
            double sum = 0;
            for (int k = i - window + 1; k <= i; k++) {
                sum += values[k];
            }
            out[i] = sum / window; // stays NaN if any month in the window is NaN
        }
        return new MonthlySeries(start, out);
    }

    /** Values of one calendar month (1 = January) between {@code from} and {@code to}, inclusive. */
    public double[] calendarMonth(int month, YearMonth from, YearMonth to) {
        return IntStream.range(0, values.length)
                .filter(i -> {
                    YearMonth ym = monthAt(i);
                    return ym.getMonthValue() == month && !ym.isBefore(from) && !ym.isAfter(to);
                })
                .mapToDouble(i -> values[i])
                .toArray();
    }
}
