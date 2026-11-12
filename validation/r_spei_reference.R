# Reference SPEI values from the R package SPEI (Vicente-Serrano et al., 2010) for the Java tests.
#
# Reads src/test/resources/basin_ppt.csv and basin_pet.csv and writes
#   src/test/resources/r_spei.csv       SPEI-1, SPEI-3 and SPEI-12 per basin
#   src/test/resources/r_glo_params.csv the fitted generalized logistic parameters
# Log-logistic (GLO) fitted by unbiased probability weighted moments, the package default,
# with 1991-2020 as reference period.
#
# Run from the repository root:  Rscript validation/r_spei_reference.R
# Tested with R 4.6.0, SPEI 1.8.1 and lmom 3.3.

suppressMessages(library(SPEI))

res <- file.path("src", "test", "resources")
ppt <- read.csv(file.path(res, "basin_ppt.csv"))
pet <- read.csv(file.path(res, "basin_pet.csv"))
basins <- c("1391", "132", "4994")
windows <- c(1, 3, 12)

out <- data.frame(month = ppt$month)
params <- data.frame()
for (b in basins) {
  wb <- ts(ppt[[paste0("ppt_", b)]] - pet[[paste0("pet_", b)]], start = c(1958, 1), frequency = 12)
  for (w in windows) {
    s <- spei(wb, scale = w, ref.start = c(1991, 1), ref.end = c(2020, 12), verbose = FALSE)
    out[[paste0("spei_", b, "_w", w)]] <- as.numeric(s$fitted)
    co <- s$coefficients
    params <- rbind(params, data.frame(basin = b, window = w, month = 1:12,
                                       xi = co["xi", 1, ], alpha = co["alpha", 1, ], kappa = co["kappa", 1, ]))
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
for (k in c("xi", "alpha", "kappa")) params[[k]] <- num(params[[k]])
write.csv(out, file.path(res, "r_spei.csv"), row.names = FALSE, quote = FALSE, eol = "\n")
write.csv(params, file.path(res, "r_glo_params.csv"), row.names = FALSE, quote = FALSE, eol = "\n")
cat("SPEI", as.character(packageVersion("SPEI")), "lmom", as.character(packageVersion("lmom")), "\n")
