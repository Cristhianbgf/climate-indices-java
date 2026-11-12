package io.github.cristhianbgf.climind.dist;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;

/** Expected values from {@code scipy.special.ndtri} (scipy 1.18.0). */
class StandardNormalTest {

    @ParameterizedTest
    @CsvSource({
        "1e-300, -37.0470962993612",
        "1e-12,  -7.034483825301131",
        "0.001,  -3.090232306167813",
        "0.02425, -1.972961051311885",
        "0.3,    -0.5244005127080409",
        "0.5,     0.0",
        "0.975,   1.959963984540054",
        "0.999999, 4.753424308817087",
    })
    void matchesScipy(double p, double expected) {
        assertEquals(expected, StandardNormal.quantile(p), 1e-13 * Math.max(1, Math.abs(expected)));
    }

    @Test
    void edges() {
        assertEquals(Double.NEGATIVE_INFINITY, StandardNormal.quantile(0));
        assertEquals(Double.POSITIVE_INFINITY, StandardNormal.quantile(1));
        assertTrue(Double.isNaN(StandardNormal.quantile(Double.NaN)));
        assertTrue(Double.isNaN(StandardNormal.quantile(-0.1)));
    }
}
