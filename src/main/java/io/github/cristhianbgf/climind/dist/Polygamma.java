package io.github.cristhianbgf.climind.dist;

/**
 * Digamma and trigamma functions for positive arguments.
 *
 * <p>Commons Math's {@code Gamma.digamma} is only accurate to about 1e-8 for arguments near 1-2,
 * which moves the maximum likelihood gamma shape by a similar amount. Here the argument is shifted
 * above 10 with the recurrence and the asymptotic series is used, which is accurate to machine
 * precision.
 */
final class Polygamma {

    private Polygamma() {}

    static double digamma(double x) {
        double shift = 0;
        while (x < 10) {
            shift -= 1 / x;
            x += 1;
        }
        double r = 1 / (x * x);
        double series = r * (1.0 / 12 - r * (1.0 / 120 - r * (1.0 / 252 - r * (1.0 / 240
                - r * (1.0 / 132 - r * (691.0 / 32760 - r / 12))))));
        return shift + Math.log(x) - 0.5 / x - series;
    }

    static double trigamma(double x) {
        double shift = 0;
        while (x < 10) {
            shift += 1 / (x * x);
            x += 1;
        }
        double r = 1 / (x * x);
        double series = (1 + r * (1.0 / 6 - r * (1.0 / 30 - r * (1.0 / 42 - r * (1.0 / 30
                - r * (5.0 / 66 - r * (691.0 / 2730 - r * 7.0 / 6))))))) / x;
        return shift + series + 0.5 * r;
    }
}
