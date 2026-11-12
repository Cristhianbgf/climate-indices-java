package io.github.cristhianbgf.climind.dist;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

class PolygammaTest {

    @Test
    void digammaKnownValues() {
        assertEquals(-0.5772156649015329, Polygamma.digamma(1), 1e-15); // minus Euler's constant
        assertEquals(-1.9635100260214235, Polygamma.digamma(0.5), 1e-15);
        assertEquals(2.2517525890667211, Polygamma.digamma(10), 1e-15);
    }

    @Test
    void trigammaKnownValues() {
        assertEquals(Math.PI * Math.PI / 6, Polygamma.trigamma(1), 1e-15);
        assertEquals(Math.PI * Math.PI / 2, Polygamma.trigamma(0.5), 1e-14);
        assertEquals(0.10516633568168575, Polygamma.trigamma(10), 1e-16);
    }
}
