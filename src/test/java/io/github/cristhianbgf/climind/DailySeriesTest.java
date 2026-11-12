package io.github.cristhianbgf.climind;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.time.LocalDate;
import org.junit.jupiter.api.Test;

class DailySeriesTest {

    /** 1999-12-30 to 2001-01-02: two days of 1999, all of 2000 (leap) and two days of 2001. */
    private static DailySeries series() {
        double[] v = new double[2 + 366 + 2];
        for (int i = 0; i < v.length; i++) {
            v[i] = i;
        }
        return new DailySeries(LocalDate.of(1999, 12, 30), v);
    }

    @Test
    void datesAndIndexes() {
        DailySeries s = series();
        assertEquals(LocalDate.of(2001, 1, 2), s.end());
        assertEquals(2, s.indexOf(LocalDate.of(2000, 1, 1)));
        assertEquals(-1, s.indexOf(LocalDate.of(1999, 12, 29)));
        assertEquals(LocalDate.of(2000, 2, 29), s.dateAt(2 + 59));
        assertEquals(1999, s.firstYear());
        assertEquals(2001, s.lastYear());
    }

    @Test
    void onlyCompleteYearsCanBeTaken() {
        DailySeries s = series();
        assertTrue(s.coversYear(2000));
        assertFalse(s.coversYear(1999));
        assertFalse(s.coversYear(2001));
        double[] y = s.year(2000);
        assertEquals(366, y.length);
        assertEquals(2, y[0], 0);
        assertEquals(367, y[365], 0);
        assertThrows(IllegalArgumentException.class, () -> s.year(2001));
    }

    @Test
    void inputArrayIsCopied() {
        double[] v = {1, 2, 3};
        DailySeries s = new DailySeries(LocalDate.of(2020, 1, 1), v);
        v[0] = 99;
        assertArrayEquals(new double[] {1, 2, 3}, s.toArray());
    }

    @Test
    void annualSeriesLooksUpByYear() {
        AnnualSeries a = new AnnualSeries(1990, new double[] {1, 2, 3});
        assertEquals(1992, a.lastYear());
        assertEquals(2, a.get(1991), 0);
        assertThrows(IllegalArgumentException.class, () -> a.get(1993));
    }
}
