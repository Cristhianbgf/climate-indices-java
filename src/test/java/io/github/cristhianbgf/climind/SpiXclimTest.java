package io.github.cristhianbgf.climind;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import io.github.cristhianbgf.climind.dist.FitMethod;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
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
        Table ppt = Table.read("/basin_ppt.csv");
        Table ref = Table.read("/xclim_spi.csv");
        MonthlySeries pr = new MonthlySeries(YearMonth.parse(ppt.months.get(0)), ppt.column("ppt_" + basin));
        double[] expected = ref.column("spi_" + basin + "_w" + window + "_" + tag);

        double[] actual = Spi.compute(pr, window, method, CAL_START, CAL_END).toArray();

        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            String at = ppt.months.get(i);
            if (Double.isNaN(expected[i])) {
                assertTrue(Double.isNaN(actual[i]), "expected NaN at " + at);
            } else {
                assertEquals(expected[i], actual[i], 1e-8, "at " + at);
            }
        }
    }

    /** Minimal reader for the CSV files in the test resources: first column is the month. */
    private record Table(List<String> header, List<String> months, List<double[]> rows) {

        static Table read(String resource) throws IOException {
            try (InputStream in = SpiXclimTest.class.getResourceAsStream(resource)) {
                List<String> lines = new String(in.readAllBytes(), StandardCharsets.UTF_8).lines().toList();
                List<String> header = List.of(lines.get(0).split(","));
                List<String> months = new ArrayList<>();
                List<double[]> rows = new ArrayList<>();
                for (String line : lines.subList(1, lines.size())) {
                    String[] f = line.split(",");
                    months.add(f[0]);
                    double[] r = new double[f.length - 1];
                    for (int k = 1; k < f.length; k++) {
                        r[k - 1] = Double.parseDouble(f[k]);
                    }
                    rows.add(r);
                }
                return new Table(header, months, rows);
            }
        }

        double[] column(String name) {
            int k = header.indexOf(name) - 1;
            if (k < 0) {
                throw new IllegalArgumentException("no column " + name);
            }
            return rows.stream().mapToDouble(r -> r[k]).toArray();
        }
    }
}
