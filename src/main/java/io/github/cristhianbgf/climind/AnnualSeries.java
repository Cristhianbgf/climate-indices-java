package io.github.cristhianbgf.climind;

/** One value per calendar year, starting at {@code firstYear}. Missing years are {@code NaN}. */
public record AnnualSeries(int firstYear, double[] values) {

    public AnnualSeries {
        values = values.clone();
    }

    @Override
    public double[] values() {
        return values.clone();
    }

    public int size() {
        return values.length;
    }

    public int lastYear() {
        return firstYear + values.length - 1;
    }

    public double get(int year) {
        int i = year - firstYear;
        if (i < 0 || i >= values.length) {
            throw new IllegalArgumentException("no value for " + year);
        }
        return values[i];
    }
}
