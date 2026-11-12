package io.github.cristhianbgf.climind.dist;

/** How the distribution parameters are estimated. */
public enum FitMethod {

    /**
     * Closed-form approximation of the maximum likelihood estimate (Thom 1958 for the gamma
     * distribution). Same as {@code method="APP"} in xclim.
     */
    APPROXIMATE,

    /** Maximum likelihood with the location fixed at zero. Same as {@code method="ML"} with {@code floc=0} in xclim. */
    MAXIMUM_LIKELIHOOD
}
