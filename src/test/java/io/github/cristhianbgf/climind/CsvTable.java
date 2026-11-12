package io.github.cristhianbgf.climind;

import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

/**
 * Minimal reader for the CSV files in the test resources: a header, a text key in the first column
 * (a month or a basin code) and numbers in the others ({@code NaN} and {@code Infinity} allowed).
 */
public record CsvTable(List<String> header, List<String> keys, List<double[]> rows) {

    public static CsvTable read(String resource) throws IOException {
        try (InputStream in = CsvTable.class.getResourceAsStream(resource)) {
            if (in == null) {
                throw new IOException("missing test resource " + resource);
            }
            List<String> lines = new String(in.readAllBytes(), StandardCharsets.UTF_8).lines().toList();
            List<String> header = List.of(lines.get(0).split(","));
            List<String> keys = new ArrayList<>();
            List<double[]> rows = new ArrayList<>();
            for (String line : lines.subList(1, lines.size())) {
                String[] f = line.split(",");
                keys.add(f[0]);
                double[] r = new double[f.length - 1];
                for (int k = 1; k < f.length; k++) {
                    r[k - 1] = Double.parseDouble(f[k]);
                }
                rows.add(r);
            }
            return new CsvTable(header, keys, rows);
        }
    }

    public int size() {
        return rows.size();
    }

    public String key(int row) {
        return keys.get(row);
    }

    public double value(int row, String column) {
        return rows.get(row)[index(column)];
    }

    public double[] column(String name) {
        int k = index(name);
        return rows.stream().mapToDouble(r -> r[k]).toArray();
    }

    private int index(String column) {
        int k = header.indexOf(column) - 1;
        if (k < 0) {
            throw new IllegalArgumentException("no column " + column);
        }
        return k;
    }
}
