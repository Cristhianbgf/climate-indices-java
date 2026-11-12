package io.github.cristhianbgf.climind.dist;

import java.util.Arrays;
import org.apache.commons.math3.special.Gamma;

/**
 * Two-parameter gamma distribution (location 0) fitted to positive values.
 *
 * <p>Zeros and missing values are ignored; the SPI handles zeros separately.
 */
public record GammaFit(double shape, double scale) {

    /** A fit with fewer than two usable values; its CDF is always {@code NaN}. */
    public static final GammaFit UNDEFINED = new GammaFit(Double.NaN, Double.NaN);

    public static GammaFit fit(double[] values, FitMethod method) {
        double[] x = Arrays.stream(values).filter(v -> v > 0 && Double.isFinite(v)).toArray();
        if (x.length < 2) {
            return UNDEFINED;
        }
        double mean = Arrays.stream(x).average().orElseThrow();
        double meanLog = Arrays.stream(x).map(Math::log).average().orElseThrow();
        double a = Math.log(mean) - meanLog;
        if (!(a > 0)) {
            return UNDEFINED; // all values equal: the shape is unbounded
        }
        double shape = (1 + Math.sqrt(1 + 4 * a / 3)) / (4 * a);
        if (method == FitMethod.MAXIMUM_LIKELIHOOD) {
            shape = solveShape(a, shape);
        }
        return new GammaFit(shape, mean / shape);
    }

    /**
     * Maximum likelihood shape: root of {@code ln(k) - digamma(k) = a}, by Newton's method
     * started from Thom's approximation. The function is convex and decreasing in {@code k},
     * so the iteration converges from that start.
     */
    private static double solveShape(double a, double start) {
        double k = start;
        for (int i = 0; i < 50; i++) {
            double f = Math.log(k) - Polygamma.digamma(k) - a;
            double df = 1 / k - Polygamma.trigamma(k);
            double next = k - f / df;
            if (next <= 0) {
                next = k / 2;
            }
            if (Math.abs(next - k) < 1e-12 * k) {
                return next;
            }
            k = next;
        }
        return k;
    }

    public double cdf(double x) {
        if (Double.isNaN(shape) || Double.isNaN(x)) {
            return Double.NaN;
        }
        return x <= 0 ? 0 : Gamma.regularizedGammaP(shape, x / scale);
    }
}
