package io.github.cristhianbgf.climind.dist;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.CsvTable;
import io.github.cristhianbgf.climind.MonthlySeries;
import java.io.IOException;
import java.time.YearMonth;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

/**
 * Maximum likelihood GLO fits of the monthly water balance of the three test basins (windows 1, 3
 * and 12, every calendar month, calibration 1991-2020), against an independent fit in Python and
 * against xclim's log-logistic fit. Reference files made by {@code validation/glo_ml_reference.py}.
 */
class LogLogisticMaximumLikelihoodTest {

    private static CsvTable ppt;
    private static CsvTable pet;

    @BeforeAll
    static void load() throws IOException {
        ppt = CsvTable.read("/basin_ppt.csv");
        pet = CsvTable.read("/basin_pet.csv");
    }

    private static double[] calibration(String basin, int window, int month) {
        double[] p = ppt.column("ppt_" + basin);
        double[] e = pet.column("pet_" + basin);
        double[] wb = new double[p.length];
        for (int i = 0; i < wb.length; i++) {
            wb[i] = p[i] - e[i];
        }
        return new MonthlySeries(YearMonth.parse(ppt.key(0)), wb).rollingMean(window)
                .calendarMonth(month, YearMonth.of(1991, 1), YearMonth.of(2020, 12));
    }

    @Test
    void matchesIndependentFit() throws IOException {
        CsvTable ref = CsvTable.read("/glo_ml_params.csv");
        for (int row = 0; row < ref.size(); row++) {
            String basin = ref.key(row);
            int window = (int) ref.value(row, "window");
            int month = (int) ref.value(row, "month");
            double[] x = calibration(basin, window, month);
            LogLogisticFit f = LogLogisticFit.fitMaximumLikelihood(x);
            String at = "basin " + basin + ", window " + window + ", month " + month;

            // the optimum is flat, so the log-likelihood is compared tightly and the parameters loosely
            assertEquals(ref.value(row, "loglik"), f.logLikelihood(x), 1e-7, at);
            assertEquals(ref.value(row, "kappa"), f.kappa(), 1e-4, at);
            assertEquals(ref.value(row, "alpha"), f.alpha(), 1e-4 * ref.value(row, "alpha"), at);
            assertEquals(ref.value(row, "xi"), f.xi(), 1e-4 * ref.value(row, "alpha"), at);
        }
    }

    /** xclim's fit is never better than this one, and is worse where its shape runs off to infinity. */
    @Test
    void neverWorseThanXclim() throws IOException {
        CsvTable xclim = CsvTable.read("/xclim_fisk_ml.csv");
        for (int row = 0; row < xclim.size(); row++) {
            double[] x = calibration(xclim.key(row), (int) xclim.value(row, "window"), (int) xclim.value(row, "month"));
            double ours = LogLogisticFit.fitMaximumLikelihood(x).logLikelihood(x);
            assertTrue(ours >= xclim.value(row, "loglik") - 1e-7,
                    "row " + row + ": " + ours + " < " + xclim.value(row, "loglik"));
        }
    }

    @Test
    void logLikelihoodOutsideTheSupportIsMinusInfinity() {
        LogLogisticFit lowerBound = new LogLogisticFit(0, 1, -0.5); // support x > -2
        assertEquals(Double.NEGATIVE_INFINITY, lowerBound.logLikelihood(new double[] {-3, 0, 1}));
        assertTrue(Double.isFinite(lowerBound.logLikelihood(new double[] {-1.9, 0, 1})));
    }

    @Test
    void tooFewValuesGiveNoFit() {
        assertSame(LogLogisticFit.UNDEFINED, LogLogisticFit.fitMaximumLikelihood(new double[] {1, 2, 3}));
    }
}
