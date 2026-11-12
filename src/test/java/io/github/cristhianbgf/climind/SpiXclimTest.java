package io.github.cristhianbgf.climind;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.dist.FitMethod;
import java.io.IOException;
import java.time.YearMonth;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

/**
 * SPI against xclim 0.61.1 on three Peruvian basins (TerraClimate, 1958-2025), calibration
 * 1991-2020. The reference file is made by {@code validation/xclim_spi_reference.py}.
 */
class SpiXclimTest {

    private static final YearMonth CAL_START = YearMonth.of(1991, 1);
    private static final YearMonth CAL_END = YearMonth.of(2020, 12);

    static Stream<Arguments> cases() {
        List<Arguments> out = new ArrayList<>();
        for (String basin : List.of("1391", "132", "4994")) {
            for (int window : List.of(1, 3, 12)) {
                out.add(Arguments.of(basin, window, FitMethod.APPROXIMATE, "APP"));
                out.add(Arguments.of(basin, window, FitMethod.MAXIMUM_LIKELIHOOD, "ML"));
            }
        }
        return out.stream();
    }

    @ParameterizedTest(name = "basin {0}, SPI-{1}, {3}")
    @MethodSource("cases")
    void matchesXclim(String basin, int window, FitMethod method, String tag) throws IOException {
        CsvTable ppt = CsvTable.read("/basin_ppt.csv");
        CsvTable ref = CsvTable.read("/xclim_spi.csv");
        MonthlySeries pr = new MonthlySeries(YearMonth.parse(ppt.key(0)), ppt.column("ppt_" + basin));
        double[] expected = ref.column("spi_" + basin + "_w" + window + "_" + tag);

        double[] actual = Spi.compute(pr, window, method, CAL_START, CAL_END).toArray();

        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            String at = ppt.key(i);
            if (Double.isNaN(expected[i])) {
                assertTrue(Double.isNaN(actual[i]), "expected NaN at " + at);
            } else {
                assertEquals(expected[i], actual[i], 1e-8, "at " + at);
            }
        }
    }
}
