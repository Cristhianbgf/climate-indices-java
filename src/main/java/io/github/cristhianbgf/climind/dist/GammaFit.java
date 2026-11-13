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

    /** Pi as written in Hosking's PELGAM. */
    private static final double HOSKING_PI = 3.1415927;

    public static GammaFit fit(double[] values, FitMethod method) {
        if (method == FitMethod.PROBABILITY_WEIGHTED_MOMENTS) {
            return fitPwm(values);
        }
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
     * Fit by unbiased probability weighted moments and L-moments, with Hosking's rational
     * approximation for the shape (routine PELGAM of his LMOMENTS package, also {@code pelgam} in
     * the R package lmom). As in the R package SPEI, fewer than four positive values give no fit.
     *
     * <p>The approximation is not the exact solution of the L-CV equation, but it is the one the R
     * package SPEI uses. It keeps Hosking's constant {@code PI = 3.1415927}, so that the shape matches
     * lmom to the last digits.
     */
    public static GammaFit fitPwm(double[] values) {
        double[] x = Arrays.stream(values).filter(v -> v > 0 && Double.isFinite(v)).sorted().toArray();
        int n = x.length;
        if (n < 4) {
            return UNDEFINED;
        }
        double b0 = 0;
        double b1 = 0;
        for (int j = 0; j < n; j++) {
            b0 += x[j];
            b1 += x[j] * j / (n - 1.0); // j is the 0-based rank in ascending order
        }
        b0 /= n;
        b1 /= n;
        double l1 = b0;
        double l2 = 2 * b1 - b0;
        if (!(l2 > 0)) {
            return UNDEFINED; // all values equal
        }
        double cv = l2 / l1;
        double shape;
        if (cv >= 0.5) {
            double t = 1 - cv;
            shape = t * (0.7213 + t * -0.5947) / (1 + t * (-2.1817 + t * 1.2113));
        } else {
            double t = HOSKING_PI * cv * cv;
            shape = (1 + -0.3080 * t) / (t * (1 + t * (-0.05812 + t * 0.01765)));
        }
        return new GammaFit(shape, l1 / shape);
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
