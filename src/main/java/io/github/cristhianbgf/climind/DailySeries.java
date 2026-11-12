package io.github.cristhianbgf.climind;

import java.time.LocalDate;
import java.time.temporal.ChronoUnit;
import java.util.Arrays;
import java.util.Objects;

/**
 * A regular daily time series on the proleptic Gregorian calendar. Missing days are stored as
 * {@code NaN}.
 *
 * <p>Instances are immutable: every operation returns a new array or series.
 */
public final class DailySeries {

    private final LocalDate start;
    private final double[] values;

    public DailySeries(LocalDate start, double[] values) {
        this.start = Objects.requireNonNull(start, "start");
        this.values = values.clone();
    }

    public LocalDate start() {
        return start;
    }

    public LocalDate end() {
        return start.plusDays(values.length - 1L);
    }

    public int size() {
        return values.length;
    }

    public double get(int i) {
        return values[i];
    }

    public LocalDate dateAt(int i) {
        return start.plusDays(i);
    }

    /** Position of a date in the series, or -1 if it falls outside. */
    public int indexOf(LocalDate date) {
        long i = start.until(date, ChronoUnit.DAYS);
        return (i < 0 || i >= values.length) ? -1 : (int) i;
    }

    public double[] toArray() {
        return values.clone();
    }

    /** First and last calendar years with at least one day in the series. */
    public int firstYear() {
        return start.getYear();
    }

    public int lastYear() {
        return end().getYear();
    }

    /** Whether the series covers every day of the year, from 1 January to 31 December. */
    public boolean coversYear(int year) {
        return !start.isAfter(LocalDate.of(year, 1, 1)) && !end().isBefore(LocalDate.of(year, 12, 31));
    }

    /** Days of one calendar year, in order. The year must be covered entirely. */
    public double[] year(int year) {
        if (!coversYear(year)) {
            throw new IllegalArgumentException("the series does not cover all of " + year);
        }
        int from = indexOf(LocalDate.of(year, 1, 1));
        int to = indexOf(LocalDate.of(year, 12, 31));
        return Arrays.copyOfRange(values, from, to + 1);
    }
}
