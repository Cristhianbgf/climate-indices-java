package io.github.cristhianbgf.climind.dist;

import java.util.Arrays;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;

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

    /**
     * Fit by maximum likelihood, with the shape kept in (-1, 1) as in {@link #fitPwm}.
     *
     * <p>In this form the logistic distribution is {@code kappa = 0}, an ordinary point of the
     * parameter space. In the usual log-logistic form (shape {@code c = -1/kappa}, as in scipy's
     * {@code fisk}) the same case sits at infinity, and samples close to it send the optimizer there.
     *
     * <p>Nelder-Mead on (xi, log alpha, atanh kappa), started from the PWM fit and restarted until the
     * log-likelihood stops improving.
     */
    public static LogLogisticFit fitMaximumLikelihood(double[] values) {
        double[] x = Arrays.stream(values).filter(Double::isFinite).toArray();
        LogLogisticFit start = fitPwm(x);
        if (start == UNDEFINED) {
            return UNDEFINED;
        }
        if (Double.isInfinite(start.logLikelihood(x))) {
            // the PWM fit leaves some values outside its support: start from the logistic instead
            double mean = Arrays.stream(x).average().orElseThrow();
            double sd = Math.sqrt(Arrays.stream(x).map(v -> (v - mean) * (v - mean)).sum() / x.length);
            start = new LogLogisticFit(mean, sd * Math.sqrt(3) / Math.PI, 0);
        }
        MultivariateFunction negLogLik = q -> {
            double ll = new LogLogisticFit(q[0], Math.exp(q[1]), Math.tanh(q[2])).logLikelihood(x);
            return Double.isFinite(ll) ? -ll : Double.MAX_VALUE;
        };
        double[] q = {start.xi, Math.log(start.alpha), atanh(start.kappa)};
        double best = negLogLik.value(q);
        SimplexOptimizer optimizer = new SimplexOptimizer(1e-15, 1e-15);
        for (int restart = 0; restart < 50; restart++) {
            PointValuePair r = optimizer.optimize(
                    new MaxEval(200_000),
                    new ObjectiveFunction(negLogLik),
                    GoalType.MINIMIZE,
                    new InitialGuess(q),
                    new NelderMeadSimplex(new double[] {0.1 * Math.exp(q[1]), 0.1, 0.1}));
            q = r.getPoint();
            if (best - r.getValue() < 1e-13) {
                break;
            }
            best = r.getValue();
        }
        return new LogLogisticFit(q[0], Math.exp(q[1]), Math.tanh(q[2]));
    }

    /** Log-likelihood of the finite values; {@code -Infinity} if any falls outside the support. */
    public double logLikelihood(double[] values) {
        if (!(alpha > 0) || !(Math.abs(kappa) < 1)) {
            return Double.NEGATIVE_INFINITY;
        }
        double sum = 0;
        for (double v : values) {
            if (!Double.isFinite(v)) {
                continue;
            }
            double u = (v - xi) / alpha;
            double y;
            if (Math.abs(kappa) < 1e-12) {
                y = u;
            } else {
                if (1 - kappa * u <= 0) {
                    return Double.NEGATIVE_INFINITY;
                }
                y = -Math.log1p(-kappa * u) / kappa;
            }
            // log(1 + exp(-y)) without overflow for large negative y
            double log1pExp = y > 0 ? Math.log1p(Math.exp(-y)) : -y + Math.log1p(Math.exp(y));
            sum += -Math.log(alpha) - (1 - kappa) * y - 2 * log1pExp;
        }
        return sum;
    }

    private static double atanh(double v) {
        return 0.5 * Math.log((1 + v) / (1 - v));
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
