# Reference SPI values and gamma fits from the R package SPEI (Vicente-Serrano et al., 2010).
#
# Reads src/test/resources/basin_ppt.csv and writes
#   src/test/resources/r_spi.csv           SPI-1, SPI-3 and SPI-12 per basin
#   src/test/resources/r_gamma_params.csv  the fitted gamma parameters, with the share of zeros
# Gamma fitted by unbiased probability weighted moments (lmom::pelgam), the package default, on the
# non-zero values of the reference period 1991-2020.
#
# SPEI 1.8.1 computes the share of zeros but does not use it: every month is transformed with the
# gamma CDF alone, so a month without rain gets -Inf. The Java tests compare the SPI only where the
# reference period of that calendar month has no zeros.
#
# Run from the repository root:  Rscript validation/r_spi_reference.R
# Tested with R 4.6.0, SPEI 1.8.1 and lmom 3.3.

suppressMessages(library(SPEI))

res <- file.path("src", "test", "resources")
ppt <- read.csv(file.path(res, "basin_ppt.csv"))
basins <- c("1391", "132", "4994")
windows <- c(1, 3, 12)

out <- data.frame(month = ppt$month)
params <- data.frame()
for (b in basins) {
  x <- ts(ppt[[paste0("ppt_", b)]], start = c(1958, 1), frequency = 12)
  for (w in windows) {
    s <- spi(x, scale = w, ref.start = c(1991, 1), ref.end = c(2020, 12), verbose = FALSE)
    out[[paste0("spi_", b, "_w", w)]] <- as.numeric(s$fitted)
    co <- s$coefficients
    acc <- stats::filter(x, rep(1, w), sides = 1)
    ref <- window(ts(as.numeric(acc), start = c(1958, 1), frequency = 12), start = c(1991, 1), end = c(2020, 12))
    zeros <- sapply(1:12, function(m) { v <- ref[cycle(ref) == m]; mean(v[!is.na(v)] == 0) })
    params <- rbind(params, data.frame(basin = b, window = w, month = 1:12,
                                       shape = co["alpha", 1, ], scale = co["beta", 1, ], zeros = zeros))
  }
}
# 17 significant digits; NA, Inf and -Inf written as Java parses them
num <- function(x) {
  y <- sprintf("%.17g", x)
  y[is.na(x)] <- "NaN"
  y[x == Inf] <- "Infinity"
  y[x == -Inf] <- "-Infinity"
  y
}
for (k in names(out)[-1]) out[[k]] <- num(out[[k]])
for (k in c("shape", "scale", "zeros")) params[[k]] <- num(params[[k]])
write.csv(out, file.path(res, "r_spi.csv"), row.names = FALSE, quote = FALSE, eol = "\n")
write.csv(params, file.path(res, "r_gamma_params.csv"), row.names = FALSE, quote = FALSE, eol = "\n")

# PWM gamma fit of the sample in GammaFitTest
X <- c(12.5, 40.1, 3.2, 77.0, 25.4, 18.9, 55.3, 9.7, 31.0, 64.2, 14.8, 22.6)
print(lmom::pelgam(lmom::samlmu(X)), digits = 17)
