package io.github.cristhianbgf.climind.etccdi;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.AnnualSeries;
import io.github.cristhianbgf.climind.CsvTable;
import io.github.cristhianbgf.climind.DailySeries;
import io.github.cristhianbgf.climind.etccdi.VeryWetDays.QuantileMethod;
import java.io.IOException;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * R95pTOT and R99pTOT on ERA5 daily precipitation at three points, 1960-2025, base period 1961-1990,
 * against an implementation in numpy and against icclim 7.3.7. The reference files are made by
 * {@code validation/etccdi_percentiles_reference.py}.
 */
class VeryWetDaysTest {

    private static final int BASE_START = 1961;
    private static final int BASE_END = 1990;

    private static DailySeries era5(String point) throws IOException {
        CsvTable daily = CsvTable.read("/era5_daily_pr.csv");
        return new DailySeries(LocalDate.parse(daily.key(0)), daily.column("pr_" + point));
    }

    static Stream<Arguments> thresholds() {
        List<Arguments> out = new ArrayList<>();
        for (String point : List.of("1391", "132", "4994")) {
            for (int p : List.of(95, 99)) {
                out.add(Arguments.of(point, p, QuantileMethod.MEDIAN_UNBIASED, "median_unbiased"));
                out.add(Arguments.of(point, p, QuantileMethod.LINEAR, "linear"));
            }
        }
        return out.stream();
    }

    @ParameterizedTest(name = "{1}th percentile at {0}, {3}")
    @MethodSource("thresholds")
    void percentileMatchesNumpy(String point, int p, QuantileMethod method, String name) throws IOException {
        CsvTable ref = CsvTable.read("/etccdi_percentile_thresholds.csv");
        int row = ref.keys().indexOf(point);
        double expected = ref.value(row, "q" + p + "_" + name);
        assertEquals(expected, VeryWetDays.wetDayPercentile(era5(point), p, BASE_START, BASE_END, method), 1e-12);
    }

    static Stream<Arguments> indices() {
        List<Arguments> out = new ArrayList<>();
        for (String point : List.of("1391", "132", "4994")) {
            for (String index : List.of("r95ptot", "r99ptot")) {
                out.add(Arguments.of(point, index));
            }
        }
        return out.stream();
    }

    @ParameterizedTest(name = "{1} at {0}")
    @MethodSource("indices")
    void matchesReferences(String point, String index) throws IOException {
        CsvTable ref = CsvTable.read("/etccdi_percentiles.csv");
        DailySeries pr = era5(point);
        boolean r95 = index.equals("r95ptot");

        AnnualSeries mm = r95 ? VeryWetDays.r95pTot(pr, BASE_START, BASE_END) : VeryWetDays.r99pTot(pr, BASE_START, BASE_END);
        AnnualSeries pct = r95 ? VeryWetDays.r95pTotPercent(pr, BASE_START, BASE_END)
                : VeryWetDays.r99pTotPercent(pr, BASE_START, BASE_END);
        compare(ref, mm, List.of("own_%s_mm_median_unbiased", "icclim_%s_mm"), index, point);
        compare(ref, pct, List.of("own_%s_pct_median_unbiased", "icclim_%s_pct"), index, point);

        double linear = VeryWetDays.wetDayPercentile(pr, r95 ? 95 : 99, BASE_START, BASE_END, QuantileMethod.LINEAR);
        compare(ref, VeryWetDays.amountAbove(pr, linear), List.of("own_%s_mm_linear"), index, point);
        compare(ref, VeryWetDays.percentAbove(pr, linear), List.of("own_%s_pct_linear"), index, point);
    }

    private static void compare(CsvTable ref, AnnualSeries actual, List<String> patterns, String index, String point) {
        assertEquals(Integer.parseInt(ref.key(0)), actual.firstYear());
        for (String pattern : patterns) {
            String column = String.format(pattern, index) + "_" + point;
            double[] expected = ref.column(column);
            assertEquals(expected.length, actual.size(), column);
            for (int i = 0; i < expected.length; i++) {
                double tolerance = 1e-9 * Math.max(1, Math.abs(expected[i]));
                assertEquals(expected[i], actual.values()[i], tolerance, column + ", " + ref.key(i));
            }
        }
    }

    @Test
    void quantileMethods() {
        double[] x = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
        assertEquals(9.1, VeryWetDays.quantile(x, 0.9, QuantileMethod.LINEAR), 1e-12);
        assertEquals(9 + 19.0 / 30, VeryWetDays.quantile(x, 0.9, QuantileMethod.MEDIAN_UNBIASED), 1e-12);
        assertEquals(10, VeryWetDays.quantile(x, 0.99, QuantileMethod.MEDIAN_UNBIASED), 0);
        assertEquals(1, VeryWetDays.quantile(x, 0, QuantileMethod.MEDIAN_UNBIASED), 0);
        assertTrue(Double.isNaN(VeryWetDays.quantile(new double[0], 0.95, QuantileMethod.LINEAR)));
    }

    @Test
    void dayEqualToTheThresholdDoesNotCount() {
        double[] v = new double[365];
        v[10] = 20;
        v[11] = 30;
        v[12] = 0.5; // not a wet day, left out of PRCPTOT
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), v);
        assertEquals(30, VeryWetDays.amountAbove(pr, 20).get(2001), 0);
        assertEquals(60, VeryWetDays.percentAbove(pr, 20).get(2001), 1e-12);
        assertEquals(50, VeryWetDays.amountAbove(pr, 0.2).get(2001), 0); // still only wet days
    }

    @Test
    void yearWithoutWetDaysHasNoPercentage() {
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), new double[365]);
        assertEquals(0, VeryWetDays.amountAbove(pr, 10).get(2001), 0);
        assertTrue(Double.isNaN(VeryWetDays.percentAbove(pr, 10).get(2001)));
        assertTrue(Double.isNaN(VeryWetDays.wetDayPercentile(pr, 95, 2001, 2001, QuantileMethod.MEDIAN_UNBIASED)));
        assertTrue(Double.isNaN(VeryWetDays.r95pTot(pr, 2001, 2001).get(2001)));
    }

    @Test
    void missingDaysAreLeftOutOfThePercentileAndBlankTheirYear() {
        double[] v = new double[730];
        for (int i = 0; i < 730; i += 7) {
            v[i] = 1 + i % 50;
        }
        v[14] = Double.NaN;
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), v);
        double q = VeryWetDays.wetDayPercentile(pr, 95, 2001, 2002, QuantileMethod.MEDIAN_UNBIASED);
        assertTrue(!Double.isNaN(q));
        assertTrue(Double.isNaN(VeryWetDays.amountAbove(pr, q).get(2001)));
        assertTrue(!Double.isNaN(VeryWetDays.amountAbove(pr, q).get(2002)));
    }

    @Test
    void basePeriodMustBeCovered() {
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), new double[730]);
        assertThrows(IllegalArgumentException.class, () -> VeryWetDays.r95pTot(pr, 2000, 2002));
        assertThrows(IllegalArgumentException.class,
                () -> VeryWetDays.wetDayPercentile(pr, 101, 2001, 2002, QuantileMethod.LINEAR));
    }
}
