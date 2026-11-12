package io.github.cristhianbgf.climind;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.dist.FitMethod;
import java.io.IOException;
import java.time.YearMonth;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * SPEI against the R package SPEI 1.8.1 (log-logistic, unbiased PWM, reference 1991-2020) on three
 * Peruvian basins, 1958-2025. The reference file is made by {@code validation/r_spei_reference.R}.
 */
class SpeiRTest {

    static Stream<Arguments> cases() {
        List<Arguments> out = new ArrayList<>();
        for (String basin : List.of("1391", "132", "4994")) {
            for (int window : List.of(1, 3, 12)) {
                out.add(Arguments.of(basin, window));
            }
        }
        return out.stream();
    }

    @ParameterizedTest(name = "basin {0}, SPEI-{1}")
    @MethodSource("cases")
    void matchesR(String basin, int window) throws IOException {
        CsvTable ppt = CsvTable.read("/basin_ppt.csv");
        CsvTable pet = CsvTable.read("/basin_pet.csv");
        CsvTable ref = CsvTable.read("/r_spei.csv");
        YearMonth start = YearMonth.parse(ppt.key(0));
        double[] expected = ref.column("spei_" + basin + "_w" + window);

        double[] actual = Spei.compute(new MonthlySeries(start, ppt.column("ppt_" + basin)),
                new MonthlySeries(start, pet.column("pet_" + basin)), window,
                YearMonth.of(1991, 1), YearMonth.of(2020, 12)).toArray();

        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            String at = ppt.key(i);
            if (Double.isNaN(expected[i])) {
                assertTrue(Double.isNaN(actual[i]), "expected NaN at " + at);
            } else if (Double.isInfinite(expected[i])) {
                // outside the support of the fitted distribution: R gives ±Inf, this library ±8.21
                assertEquals(Math.signum(expected[i]) * Spi.LIMIT, actual[i], 0, "at " + at);
            } else {
                assertEquals(expected[i], actual[i], 1e-8, "at " + at);
            }
        }
    }

    @Test
    void maximumLikelihoodRunsAndApproximateIsRejected() throws IOException {
        CsvTable ppt = CsvTable.read("/basin_ppt.csv");
        CsvTable pet = CsvTable.read("/basin_pet.csv");
        YearMonth start = YearMonth.parse(ppt.key(0));
        MonthlySeries p = new MonthlySeries(start, ppt.column("ppt_132"));
        MonthlySeries e = new MonthlySeries(start, pet.column("pet_132"));
        YearMonth from = YearMonth.of(1991, 1);
        YearMonth to = YearMonth.of(2020, 12);

        double[] ml = Spei.compute(p, e, 3, FitMethod.MAXIMUM_LIKELIHOOD, from, to).toArray();
        assertTrue(Double.isNaN(ml[1]));
        for (int i = 2; i < ml.length; i++) {
            assertTrue(Math.abs(ml[i]) <= Spi.LIMIT, "at " + ppt.key(i));
        }
        assertThrows(IllegalArgumentException.class, () -> Spei.compute(p, e, 3, FitMethod.APPROXIMATE, from, to));
    }

    @Test
    void rejectsMismatchedSeries() {
        MonthlySeries p = new MonthlySeries(YearMonth.of(2000, 1), new double[24]);
        MonthlySeries e = new MonthlySeries(YearMonth.of(2000, 2), new double[24]);
        assertThrows(IllegalArgumentException.class, () -> Spei.compute(p, e, 1));
    }
}
