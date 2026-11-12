package io.github.cristhianbgf.climind.dist;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.CsvTable;
import io.github.cristhianbgf.climind.MonthlySeries;
import java.io.IOException;
import java.time.YearMonth;
import org.junit.jupiter.api.Test;

class LogLogisticFitTest {

    /**
     * Parameters fitted by the R package SPEI 1.8.1 (unbiased PWM, reference 1991-2020) to the
     * monthly water balance of three basins, for every calendar month and window 1, 3 and 12.
     * Made by {@code validation/r_spei_reference.R}.
     */
    @Test
    void matchesRPackageSpei() throws IOException {
        CsvTable ppt = CsvTable.read("/basin_ppt.csv");
        CsvTable pet = CsvTable.read("/basin_pet.csv");
        CsvTable ref = CsvTable.read("/r_glo_params.csv");
        YearMonth start = YearMonth.parse(ppt.key(0));

        for (int row = 0; row < ref.size(); row++) {
            String basin = ref.key(row);
            int window = (int) ref.value(row, "window");
            int month = (int) ref.value(row, "month");

            double[] p = ppt.column("ppt_" + basin);
            double[] e = pet.column("pet_" + basin);
            double[] wb = new double[p.length];
            for (int i = 0; i < wb.length; i++) {
                wb[i] = p[i] - e[i];
            }
            // R accumulates sums and this library means; the fit scales with the data, so xi and
            // alpha are compared after multiplying by the window
            double[] cal = new MonthlySeries(start, wb).rollingMean(window)
                    .calendarMonth(month, YearMonth.of(1991, 1), YearMonth.of(2020, 12));
            LogLogisticFit f = LogLogisticFit.fitPwm(cal);

            String at = "basin " + basin + ", window " + window + ", month " + month;
            assertEquals(ref.value(row, "xi"), f.xi() * window, 1e-9 * Math.abs(ref.value(row, "xi")) + 1e-9, at);
            assertEquals(ref.value(row, "alpha"), f.alpha() * window, 1e-9 * ref.value(row, "alpha"), at);
            assertEquals(ref.value(row, "kappa"), f.kappa(), 1e-12, at);
        }
    }

    @Test
    void cdfIsZeroOrOneOutsideTheSupport() {
        LogLogisticFit lowerBound = new LogLogisticFit(0, 1, -0.5); // support x > -2
        assertEquals(0, lowerBound.cdf(-3), 0);
        assertTrue(lowerBound.cdf(-1.9) > 0);
        LogLogisticFit upperBound = new LogLogisticFit(0, 1, 0.5); // support x < 2
        assertEquals(1, upperBound.cdf(3), 0);
    }

    @Test
    void zeroKappaIsTheLogisticDistribution() {
        LogLogisticFit logistic = new LogLogisticFit(10, 2, 0);
        assertEquals(0.5, logistic.cdf(10), 1e-15);
        assertEquals(1 / (1 + Math.exp(-1.5)), logistic.cdf(13), 1e-15);
    }

    @Test
    void tooFewOrConstantValuesGiveNoFit() {
        assertSame(LogLogisticFit.UNDEFINED, LogLogisticFit.fitPwm(new double[] {1, 2, 3}));
        assertSame(LogLogisticFit.UNDEFINED, LogLogisticFit.fitPwm(new double[] {4, 4, 4, 4, 4}));
        assertTrue(Double.isNaN(LogLogisticFit.UNDEFINED.cdf(0)));
    }
}
