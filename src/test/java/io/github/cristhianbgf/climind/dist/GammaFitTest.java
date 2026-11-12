package io.github.cristhianbgf.climind.dist;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertSame;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

/** Expected values from scipy 1.18.0: {@code scipy.stats.gamma.fit(x, floc=0)} and Thom's formula. */
class GammaFitTest {

    private static final double[] X = {12.5, 40.1, 3.2, 77.0, 25.4, 18.9, 55.3, 9.7, 31.0, 64.2, 14.8, 22.6};

    @Test
    void maximumLikelihoodMatchesScipy() {
        GammaFit g = GammaFit.fit(X, FitMethod.MAXIMUM_LIKELIHOOD);
        assertEquals(1.7686037076279508, g.shape(), 1e-9);
        assertEquals(17.655170497114323, g.scale(), 1e-8);
        assertEquals(0.5796473361599042, g.cdf(30.0), 1e-10);
    }

    @Test
    void approximationMatchesThom() {
        GammaFit g = GammaFit.fit(X, FitMethod.APPROXIMATE);
        assertEquals(1.7725672629502307, g.shape(), 1e-12);
        assertEquals(17.61569259043499, g.scale(), 1e-10);
    }

    @Test
    void zerosAndMissingValuesAreIgnored() {
        double[] withGaps = {0, 12.5, 40.1, Double.NaN, 3.2, 77.0, 25.4, 18.9, 55.3, 9.7, 31.0, 0, 64.2, 14.8, 22.6};
        GammaFit a = GammaFit.fit(withGaps, FitMethod.MAXIMUM_LIKELIHOOD);
        GammaFit b = GammaFit.fit(X, FitMethod.MAXIMUM_LIKELIHOOD);
        assertEquals(b.shape(), a.shape(), 1e-14);
        assertEquals(b.scale(), a.scale(), 1e-12);
    }

    @Test
    void tooFewOrConstantValuesGiveNoFit() {
        assertSame(GammaFit.UNDEFINED, GammaFit.fit(new double[] {5.0}, FitMethod.APPROXIMATE));
        assertSame(GammaFit.UNDEFINED, GammaFit.fit(new double[] {0, 0, Double.NaN}, FitMethod.APPROXIMATE));
        assertSame(GammaFit.UNDEFINED, GammaFit.fit(new double[] {4, 4, 4}, FitMethod.MAXIMUM_LIKELIHOOD));
        assertTrue(Double.isNaN(GammaFit.UNDEFINED.cdf(1.0)));
    }

    @Test
    void cdfIsZeroAtOrBelowZero() {
        GammaFit g = GammaFit.fit(X, FitMethod.APPROXIMATE);
        assertEquals(0, g.cdf(0), 0);
        assertEquals(0, g.cdf(-1), 0);
    }
}
