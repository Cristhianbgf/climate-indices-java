package io.github.cristhianbgf.climind.dist;

import java.util.Arrays;

/**
 * Three-parameter log-logistic distribution, written as Hosking's generalized logistic (GLO) with
 * location {@code xi}, scale {@code alpha} and shape {@code kappa}. This is the form used by the R
 * package SPEI (Vicente-Serrano et al., 2010).
 *
 * <p>With {@code kappa = 0} it is the logistic distribution; with {@code kappa < 0} it has a lower
 * bound and with {@code kappa > 0} an upper bound at {@code xi + alpha / kappa}.
 */
public record LogLogisticFit(double xi, double alpha, double kappa) {

    /** A fit with too few or constant values; its CDF is always {@code NaN}. */
    public static final LogLogisticFit UNDEFINED = new LogLogisticFit(Double.NaN, Double.NaN, Double.NaN);

    /** Below this |kappa| the distribution is taken as logistic, as in Hosking's {@code pelglo}. */
    private static final double SMALL = 1e-6;

    /**
     * Fit by unbiased probability weighted moments and L-moments (Hosking, 1990), the default of
     * the R package SPEI. Missing values are ignored; fewer than four values give no fit.
     */
    public static LogLogisticFit fitPwm(double[] values) {
        double[] x = Arrays.stream(values).filter(Double::isFinite).sorted().toArray();
        int n = x.length;
        if (n < 4) {
            return UNDEFINED;
        }
        double b0 = 0;
        double b1 = 0;
        double b2 = 0;
        for (int j = 0; j < n; j++) {
            // weights of the unbiased estimators, with j the 0-based rank in ascending order
            b0 += x[j];
            b1 += x[j] * j / (n - 1.0);
            b2 += x[j] * j * (j - 1.0) / ((n - 1.0) * (n - 2.0));
        }
        b0 /= n;
        b1 /= n;
        b2 /= n;
        double l1 = b0;
        double l2 = 2 * b1 - b0;
        double l3 = 6 * b2 - 6 * b1 + b0;
        if (!(l2 > 0)) {
            return UNDEFINED; // constant sample
        }
        double t3 = l3 / l2;

        double k = -t3;
        if (Math.abs(k) <= SMALL) {
            return new LogLogisticFit(l1, l2, 0);
        }
        if (Math.abs(k) >= 1) {
            return UNDEFINED;
        }
        double gg = k * Math.PI / Math.sin(k * Math.PI);
        double a = l2 / gg;
        return new LogLogisticFit(l1 - a * (1 - gg) / k, a, k);
    }

    public double cdf(double x) {
        if (Double.isNaN(alpha) || Double.isNaN(x)) {
            return Double.NaN;
        }
        double y;
        if (kappa == 0) {
            y = (x - xi) / alpha;
        } else {
            // outside the support the log is -infinity and the CDF is exactly 0 or 1
            y = -Math.log(Math.max(0, 1 - kappa * (x - xi) / alpha)) / kappa;
        }
        return 1 / (1 + Math.exp(-y));
    }
}
