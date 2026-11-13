package io.github.cristhianbgf.climind;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.dist.FitMethod;
import io.github.cristhianbgf.climind.dist.GammaFit;
import java.io.IOException;
import java.time.YearMonth;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * SPI with the gamma fitted by probability weighted moments, against the R package SPEI 1.8.1
 * (reference 1991-2020) on three Peruvian basins, 1958-2025. The reference files are made by
 * {@code validation/r_spi_reference.R}.
 *
 * <p>The parameters are compared in all 108 fits. The index is compared only in the calendar months
 * whose reference period has no zeros: SPEI 1.8.1 leaves the share of zeros out of the transformation,
 * while this library, like xclim, adds it. The whole index, zeros included, is compared with an
 * implementation in numpy.
 */
class SpiRTest {

    private static final YearMonth CAL_START = YearMonth.of(1991, 1);
    private static final YearMonth CAL_END = YearMonth.of(2020, 12);

    static Stream<Arguments> cases() {
        List<Arguments> out = new ArrayList<>();
        for (String basin : List.of("1391", "132", "4994")) {
            for (int window : List.of(1, 3, 12)) {
                out.add(Arguments.of(basin, window));
            }
        }
        return out.stream();
    }

    @ParameterizedTest(name = "basin {0}, SPI-{1}")
    @MethodSource("cases")
    void matchesR(String basin, int window) throws IOException {
        CsvTable ppt = CsvTable.read("/basin_ppt.csv");
        CsvTable params = CsvTable.read("/r_gamma_params.csv");
        CsvTable ref = CsvTable.read("/r_spi.csv");
        MonthlySeries pr = new MonthlySeries(YearMonth.parse(ppt.key(0)), ppt.column("ppt_" + basin));
        MonthlySeries acc = pr.rollingMean(window);

        boolean[] noZeros = new boolean[13];
        for (int row = 0; row < params.size(); row++) {
            if (!params.key(row).equals(basin) || params.value(row, "window") != window) {
                continue;
            }
            int m = (int) params.value(row, "month");
            GammaFit fit = GammaFit.fitPwm(acc.calendarMonth(m, CAL_START, CAL_END));
            // R fits sums and this library means: the scale grows with the window, the shape does not
            assertEquals(params.value(row, "shape"), fit.shape(), 1e-12 * params.value(row, "shape"), "shape, month " + m);
            assertEquals(params.value(row, "scale"), fit.scale() * window, 1e-12 * params.value(row, "scale"),
                    "scale, month " + m);
            noZeros[m] = params.value(row, "zeros") == 0;
        }

        double[] expected = ref.column("spi_" + basin + "_w" + window);
        double[] actual = Spi.compute(pr, window, FitMethod.PROBABILITY_WEIGHTED_MOMENTS, CAL_START, CAL_END).toArray();
        assertEquals(expected.length, actual.length);
        int compared = 0;
        for (int i = 0; i < expected.length; i++) {
            String at = ppt.key(i);
            if (!noZeros[acc.monthAt(i).getMonthValue()]) {
                continue;
            }
            if (Double.isNaN(expected[i])) {
                assertTrue(Double.isNaN(actual[i]), "expected NaN at " + at);
            } else if (Double.isInfinite(expected[i])) {
                assertEquals(Math.signum(expected[i]) * Spi.LIMIT, actual[i], 0, "at " + at);
            } else {
                assertEquals(expected[i], actual[i], 1e-8, "at " + at);
            }
            compared++;
        }
        assertTrue(compared > 0);
    }

    /** The whole SPI, zeros included, against the numpy version in {@code validation/fit_comparison.py}. */
    @ParameterizedTest(name = "basin {0}, SPI-{1}")
    @MethodSource("cases")
    void matchesIndependentImplementation(String basin, int window) throws IOException {
        CsvTable ppt = CsvTable.read("/basin_ppt.csv");
        CsvTable ref = CsvTable.read("/spi_pwm.csv");
        MonthlySeries pr = new MonthlySeries(YearMonth.parse(ppt.key(0)), ppt.column("ppt_" + basin));
        double[] expected = ref.column("spi_" + basin + "_w" + window);
        double[] actual = Spi.compute(pr, window, FitMethod.PROBABILITY_WEIGHTED_MOMENTS, CAL_START, CAL_END).toArray();
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            if (Double.isNaN(expected[i])) {
                assertTrue(Double.isNaN(actual[i]), "expected NaN at " + ppt.key(i));
            } else {
                assertEquals(expected[i], actual[i], 1e-8, "at " + ppt.key(i));
            }
        }
    }
}
