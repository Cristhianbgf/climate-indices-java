package io.github.cristhianbgf.climind.dist;

/** How the distribution parameters are estimated. */
public enum FitMethod {

    /**
     * Closed-form approximation of the maximum likelihood estimate (Thom 1958 for the gamma
     * distribution). Same as {@code method="APP"} in xclim. Gamma only.
     */
    APPROXIMATE,

    /**
     * Maximum likelihood. For the gamma distribution the location is fixed at zero, as with
     * {@code method="ML"} and {@code floc=0} in xclim.
     */
    MAXIMUM_LIKELIHOOD,

    /**
     * Unbiased probability weighted moments and L-moments (Hosking, 1990), the default of the R
     * package SPEI for both the log-logistic and the gamma distribution.
     */
    PROBABILITY_WEIGHTED_MOMENTS
}
