package io.github.cristhianbgf.climind;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.dist.FitMethod;
import io.github.cristhianbgf.climind.dist.StandardNormal;
import java.time.YearMonth;
import org.junit.jupiter.api.Test;

class SpiTest {

    /** Ten years of a single calendar month repeated; every other month is filled the same way. */
    private static MonthlySeries tenYears(double[] perYear) {
        double[] v = new double[perYear.length * 12];
        for (int y = 0; y < perYear.length; y++) {
            for (int m = 0; m < 12; m++) {
                v[y * 12 + m] = perYear[y];
            }
        }
        return new MonthlySeries(YearMonth.of(2000, 1), v);
    }

    @Test
    void dryWindowsTakeTheProbabilityOfZero() {
        // 3 of 10 years are dry, so a dry month gets p = 0.3
        MonthlySeries pr = tenYears(new double[] {0, 12, 30, 0, 8, 22, 0, 15, 40, 18});
        double[] spi = Spi.compute(pr, 1, FitMethod.MAXIMUM_LIKELIHOOD).toArray();
        assertEquals(StandardNormal.quantile(0.3), spi[0], 1e-12);
        assertEquals(spi[0], spi[36], 0);
        assertTrue(spi[12] > spi[0]);
    }

    @Test
    void leadingMonthsOfTheWindowAreMissing() {
        MonthlySeries pr = tenYears(new double[] {5, 12, 30, 7, 8, 22, 3, 15, 40, 18});
        double[] spi = Spi.compute(pr, 3, FitMethod.APPROXIMATE).toArray();
        assertTrue(Double.isNaN(spi[0]));
        assertTrue(Double.isNaN(spi[1]));
        assertTrue(Double.isFinite(spi[2]));
    }

    @Test
    void valuesAreLimited() {
        assertEquals(-Spi.LIMIT, Spi.toStandardNormal(0), 0);
        assertEquals(Spi.LIMIT, Spi.toStandardNormal(1), 0);
    }

    @Test
    void rejectsReversedCalibration() {
        MonthlySeries pr = tenYears(new double[] {5, 12, 30, 7, 8, 22, 3, 15, 40, 18});
        assertThrows(IllegalArgumentException.class, () ->
                Spi.compute(pr, 1, FitMethod.APPROXIMATE, YearMonth.of(2005, 1), YearMonth.of(2004, 12)));
    }
}
