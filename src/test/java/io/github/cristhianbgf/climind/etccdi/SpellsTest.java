package io.github.cristhianbgf.climind.etccdi;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.AnnualSeries;
import io.github.cristhianbgf.climind.CsvTable;
import io.github.cristhianbgf.climind.DailySeries;
import io.github.cristhianbgf.climind.etccdi.Spells.YearRule;
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
 * CDD and CWD on ERA5 daily precipitation at three points, 1960-2025, against an implementation in
 * Python written from the ETCCDI definitions, and against xclim 1.0.0 and icclim 7.3.7. The reference
 * file is made by {@code validation/etccdi_spells_reference.py}.
 */
class SpellsTest {

    static Stream<Arguments> cases() {
        List<Arguments> out = new ArrayList<>();
        for (String point : List.of("1391", "132", "4994")) {
            for (String index : List.of("cdd", "cwd")) {
                out.add(Arguments.of(point, index, YearRule.SPLIT, List.of("own_%s_split", "xclim_%s_split")));
                out.add(Arguments.of(point, index, YearRule.YEAR_OF_START,
                        List.of("own_%s_start", "xclim_%s_start", "icclim_%s")));
                out.add(Arguments.of(point, index, YearRule.YEAR_OF_END, List.of("own_%s_end")));
            }
        }
        return out.stream();
    }

    @ParameterizedTest(name = "{1} at {0}, {2}")
    @MethodSource("cases")
    void matchesReferences(String point, String index, YearRule rule, List<String> columns) throws IOException {
        CsvTable daily = CsvTable.read("/era5_daily_pr.csv");
        CsvTable ref = CsvTable.read("/etccdi_spells.csv");
        DailySeries pr = new DailySeries(LocalDate.parse(daily.key(0)), daily.column("pr_" + point));

        AnnualSeries actual = index.equals("cdd") ? Spells.cdd(pr, rule) : Spells.cwd(pr, rule);

        assertEquals(Integer.parseInt(ref.key(0)), actual.firstYear());
        for (String pattern : columns) {
            String column = String.format(pattern, index) + "_" + point;
            double[] expected = ref.column(column);
            assertEquals(expected.length, actual.size(), column);
            for (int i = 0; i < expected.length; i++) {
                assertEquals(expected[i], actual.values()[i], 0, column + ", " + ref.key(i));
            }
        }
    }

    /** Ten wet days from 27 December 2001 to 5 January 2002, inside two otherwise dry years. */
    private static DailySeries spellAcrossNewYear() {
        double[] v = new double[730];
        for (int i = 360; i < 370; i++) {
            v[i] = 5;
        }
        return new DailySeries(LocalDate.of(2001, 1, 1), v);
    }

    @Test
    void yearRules() {
        DailySeries pr = spellAcrossNewYear();
        assertEquals(5, Spells.cwd(pr, YearRule.SPLIT).get(2001), 0);
        assertEquals(5, Spells.cwd(pr, YearRule.SPLIT).get(2002), 0);
        assertEquals(10, Spells.cwd(pr, YearRule.YEAR_OF_START).get(2001), 0);
        assertEquals(0, Spells.cwd(pr, YearRule.YEAR_OF_START).get(2002), 0);
        assertEquals(0, Spells.cwd(pr, YearRule.YEAR_OF_END).get(2001), 0);
        assertEquals(10, Spells.cwd(pr, YearRule.YEAR_OF_END).get(2002), 0);
    }

    @Test
    void oneMillimetreIsWet() {
        double[] v = new double[365];
        v[100] = 1.0;
        v[101] = 1.0;
        v[102] = 0.99;
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), v);
        assertEquals(2, Spells.cwd(pr).get(2001), 0);
        assertEquals(263, Spells.cdd(pr).get(2001), 0); // days 102 to 364
    }

    @Test
    void missingDayEndsTheSpellAndBlanksItsYear() {
        double[] v = new double[730];
        v[500] = Double.NaN;
        DailySeries pr = new DailySeries(LocalDate.of(2001, 1, 1), v);
        assertEquals(365, Spells.cdd(pr).get(2001), 0);
        assertTrue(Double.isNaN(Spells.cdd(pr).get(2002)));
    }
}
