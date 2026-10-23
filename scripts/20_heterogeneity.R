# 20: Exposure x land/water contiguity, exposure x embedding (full panel), and exposure x log capital distance, each with explicit interaction CIs. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(cshapes) })
src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1, T = as.integer(tri))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
lin <- function(b, V, w) { est <- sum(b[names(w)] * w); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); c(est, se) }
fit <- function(x, y, dur, mod, label, extra = "") {
  x$l_dur <- log(x[[dur]]); f <- as.formula(paste(y, "~ z_mx *", mod, "+ z_mn + z_cap + l_dur + l_same + l_diff", if (mod != "contig") "+ contig" else "", extra))
  m <- glm(f, x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  inter <- paste0("z_mx:", mod); main <- lin(b, V, c(z_mx = 1)); it <- lin(b, V, setNames(1, inter)); both <- lin(b, V, setNames(c(1, 1), c("z_mx", inter)))
  cat(sprintf("%-52s base (%s=0) %s | interaction %s | implied at %s=1 %s | N=%d events=%d\n", label, mod, ci(main[1], main[2]), ci(it[1], it[2]), mod, ci(both[1], both[2]), nrow(x), sum(x[[y]])))
}
# capital distance for dyad-years of the embedded onset/ending samples (1886+)
yrs <- sort(unique(d$year[d$tri & d$year >= 1886]))
cap <- bind_rows(lapply(yrs, function(y) { mm <- distmatrix(as.Date(paste0(y, "-01-01")), type = "capdist", useGW = FALSE); ids <- as.numeric(rownames(mm)); ix <- which(upper.tri(mm), arr.ind = TRUE)
  tibble(year = y, lo = pmin(ids[ix[, 1]], ids[ix[, 2]]), hi = pmax(ids[ix[, 1]], ids[ix[, 2]]), capdist = mm[ix]) }))
d <- d %>% left_join(cap, by = c("year", "lo", "hi")) %>% mutate(z_lcap = as.numeric(scale(log(capdist + 1))))
sink("results/heterogeneity_output.txt", split = TRUE)
cat("(a) Embedded pairs: larger share x contiguity (land or water, any type)\n")
for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE")
  fit(d %>% filter(tri, m == 0), "onset_m", "dur_m", "contig", paste("onset,", s), sp); fit(d %>% filter(tri, m == 1), "end_m", "dur_m", "contig", paste("ending,", s), sp) }
cat("\n(b) Full panel: larger share x embedded pair (T)\n")
for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE")
  fit(d %>% filter(m == 0), "onset_m", "dur_m", "T", paste("onset,", s), sp); fit(d %>% filter(m == 1), "end_m", "dur_m", "T", paste("ending,", s), sp) }
cat("\n(c) Embedded pairs, 1886+: larger share x standardised log capital distance (controls include contiguity)\n")
x <- d %>% filter(tri, !is.na(capdist))
for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE"); fit(x %>% filter(m == 0), "onset_m", "dur_m", "z_lcap", paste("onset,", s), sp); fit(x %>% filter(m == 1), "end_m", "dur_m", "z_lcap", paste("ending,", s), sp) }
sink()
