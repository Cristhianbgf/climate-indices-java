package io.github.cristhianbgf.climind;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.time.YearMonth;
import org.junit.jupiter.api.Test;

class MonthlySeriesTest {

    private final MonthlySeries s =
            new MonthlySeries(YearMonth.of(2000, 11), new double[] {1, 2, 3, 4, Double.NaN, 6});

    @Test
    void endAndIndexing() {
        assertEquals(YearMonth.of(2001, 4), s.end());
        assertEquals(2, s.indexOf(YearMonth.of(2001, 1)));
        assertEquals(-1, s.indexOf(YearMonth.of(2000, 10)));
        assertEquals(-1, s.indexOf(YearMonth.of(2001, 5)));
    }

    @Test
    void rollingMeanLeavesLeadingAndGapWindowsMissing() {
        double[] r = s.rollingMean(3).toArray();
        assertTrue(Double.isNaN(r[0]));
        assertTrue(Double.isNaN(r[1]));
        assertEquals(2.0, r[2], 1e-12);
        assertEquals(3.0, r[3], 1e-12);
        assertTrue(Double.isNaN(r[4]));
        assertTrue(Double.isNaN(r[5]));
    }

    @Test
    void windowOfOneIsTheSeriesItself() {
        assertArrayEquals(s.toArray(), s.rollingMean(1).toArray());
        assertThrows(IllegalArgumentException.class, () -> s.rollingMean(0));
    }

    @Test
    void calendarMonthRespectsPeriod() {
        MonthlySeries years = new MonthlySeries(YearMonth.of(1990, 1), new double[36]);
        assertEquals(3, years.calendarMonth(5, YearMonth.of(1990, 1), YearMonth.of(1992, 12)).length);
        assertEquals(1, years.calendarMonth(5, YearMonth.of(1991, 1), YearMonth.of(1991, 12)).length);
    }

    @Test
    void inputArrayIsCopied() {
        double[] v = {1, 2};
        MonthlySeries m = new MonthlySeries(YearMonth.of(2020, 1), v);
        v[0] = 99;
        assertEquals(1, m.get(0), 0);
    }
}
