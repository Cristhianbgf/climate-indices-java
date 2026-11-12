package io.github.cristhianbgf.climind.dist;

import org.apache.commons.math3.special.Erf;

/**
 * Quantile function of the standard normal distribution, accurate in the tails.
 *
 * <p>Commons Math inverts the normal CDF through {@code erfInv(2p - 1)}, which loses most of its
 * precision when {@code p} is close to 0 (around 1e-12, an SPI near -7). Here Acklam's rational
 * approximation is refined with Halley's method on the complementary error function, which keeps
 * full relative precision in both tails.
 */
public final class StandardNormal {

    private static final double[] A = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
        1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    private static final double[] B = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
        6.680131188771972e+01, -1.328068155288572e+01};
    private static final double[] C = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
        -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    private static final double[] D = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
        3.754408661907416e+00};
    private static final double P_LOW = 0.02425;

    private StandardNormal() {}

    public static double quantile(double p) {
        if (Double.isNaN(p) || p < 0 || p > 1) {
            return Double.NaN;
        }
        if (p == 0) {
            return Double.NEGATIVE_INFINITY;
        }
        if (p == 1) {
            return Double.POSITIVE_INFINITY;
        }
        // The lower half is computed directly; the upper half by symmetry, so that p close to 1
        // is handled through 1 - p, which is exact for the values that reach this point.
        return p <= 0.5 ? lower(p) : -lower(1 - p);
    }

    private static double lower(double p) {
        double x;
        if (p < P_LOW) {
            double q = Math.sqrt(-2 * Math.log(p));
            x = (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
                    / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1);
        } else {
            double q = p - 0.5;
            double r = q * q;
            x = (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
                    / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1);
        }
        for (int i = 0; i < 2; i++) {
            double e = 0.5 * Erf.erfc(-x / Math.sqrt(2)) - p;
            double u = e * Math.sqrt(2 * Math.PI) * Math.exp(x * x / 2);
            x = x - u / (1 + x * u / 2);
        }
        return x;
    }
}
