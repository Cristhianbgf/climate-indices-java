package io.github.cristhianbgf.climind.etccdi;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.AnnualSeries;
import io.github.cristhianbgf.climind.CsvTable;
import io.github.cristhianbgf.climind.DailySeries;
import java.io.IOException;
import java.time.LocalDate;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Stream;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * Amount indices against xclim 0.61.1 on ERA5 daily precipitation at three points, 1960-2025.
 * The reference file is made by {@code validation/xclim_etccdi_reference.py}.
 */
class PrecipitationAmountsTest {

    private static final List<String> POINTS = List.of("1391", "132", "4994");

    static Stream<Arguments> cases() {
        List<Arguments> out = new ArrayList<>();
        List<Object[]> indices = List.of(
                new Object[] {"prcptot", (Function<DailySeries, AnnualSeries>) PrecipitationAmounts::prcptot},
                new Object[] {"rx1day", (Function<DailySeries, AnnualSeries>) PrecipitationAmounts::rx1day},
                new Object[] {"rx5day", (Function<DailySeries, AnnualSeries>) PrecipitationAmounts::rx5day},
                new Object[] {"sdii", (Function<DailySeries, AnnualSeries>) PrecipitationAmounts::sdii},
                new Object[] {"r10mm", (Function<DailySeries, AnnualSeries>) PrecipitationAmounts::r10mm},
                new Object[] {"r20mm", (Function<DailySeries, AnnualSeries>) PrecipitationAmounts::r20mm});
        for (String point : POINTS) {
            for (Object[] idx : indices) {
                out.add(Arguments.of(idx[0], point, idx[1]));
            }
        }
        return out.stream();
    }

    @ParameterizedTest(name = "{0} at {1}")
    @MethodSource("cases")
    void matchesXclim(String name, String point, Function<DailySeries, AnnualSeries> index) throws IOException {
        CsvTable daily = CsvTable.read("/era5_daily_pr.csv");
        CsvTable ref = CsvTable.read("/xclim_etccdi.csv");
        DailySeries pr = new DailySeries(LocalDate.parse(daily.key(0)), daily.column("pr_" + point));
        double[] expected = ref.column(name + "_" + point);

        AnnualSeries actual = index.apply(pr);

        assertEquals(Integer.parseInt(ref.key(0)), actual.firstYear());
        assertEquals(expected.length, actual.size());
        for (int i = 0; i < expected.length; i++) {
            // the reference is written with 10 significant digits
            assertEquals(expected[i], actual.values()[i], 1e-9 * Math.max(1, Math.abs(expected[i])),
                    name + " at " + point + ", " + ref.key(i));
        }
    }

    @Test
    void thresholdsAreInclusive() {
        double[] v = new double[365];
        v[10] = 1.0;
        v[11] = 0.99;
        v[12] = 10.0;
        v[13] = 20.0;
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), v);
        assertEquals(31.0, PrecipitationAmounts.prcptot(pr).get(2001), 0);
        assertEquals(31.0 / 3, PrecipitationAmounts.sdii(pr).get(2001), 1e-15);
        assertEquals(2, PrecipitationAmounts.r10mm(pr).get(2001), 0);
        assertEquals(1, PrecipitationAmounts.r20mm(pr).get(2001), 0);
    }

    @Test
    void partialYearsAreLeftOutAndMissingDaysGiveNaN() {
        double[] v = new double[365 + 365 + 10];
        v[400] = Double.NaN;
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), v); // ends on 10 January 2003
        AnnualSeries tot = PrecipitationAmounts.prcptot(pr);
        assertEquals(2001, tot.firstYear());
        assertEquals(2002, tot.lastYear());
        assertEquals(0, tot.get(2001), 0);
        assertTrue(Double.isNaN(tot.get(2002)));
        assertTrue(Double.isNaN(PrecipitationAmounts.sdii(pr).get(2001)));
    }

    @Test
    void rx5dayWindowCanStartInDecember() {
        double[] v = new double[365 + 365];
        v[363] = 30; // 30 December 2001
        v[364] = 30; // 31 December 2001
        v[365] = 30; // 1 January 2002
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), v);
        assertEquals(60, PrecipitationAmounts.rx5day(pr).get(2001), 0);
        assertEquals(90, PrecipitationAmounts.rx5day(pr).get(2002), 0);
    }
}
