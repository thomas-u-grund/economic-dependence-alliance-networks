# 84: Who are the distant exposed pairs? Interactions of the larger share with major-power involvement and same-region status among non-contiguous pairs. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
nmc <- read_csv("raw/nmc/NMCv7/abridged/NMC-70-abridged.csv", show_col_types = FALSE) %>% select(ccode, year, cinc) %>% group_by(year) %>% mutate(rk = rank(-cinc, ties.method = "min")) %>% ungroup()
MAJ <- c(2, 200, 220, 255, 300, 325, 365, 710, 740)
region <- function(cc) cut(cc, c(0, 199, 399, 626, 699, 999), labels = c("Americas", "Europe", "Africa", "MiddleEast", "Asia"))
d <- d %>% left_join(nmc %>% transmute(year, lo = ccode, rk_lo = rk), by = c("year", "lo")) %>% left_join(nmc %>% transmute(year, hi = ccode, rk_hi = rk), by = c("year", "hi")) %>%
  mutate(top5 = as.integer(coalesce(rk_lo <= 5, FALSE) | coalesce(rk_hi <= 5, FALSE)), ever_major = as.integer(lo %in% MAJ | hi %in% MAJ), same_region = as.integer(region(lo) == region(hi)),
         onset_m = as.integer(m1 == 1), form_a = as.integer(a1 == 1), tri = (n_same + n_diff) >= 1) %>% filter(contig == 0)
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
lin <- function(b, V, w) { est <- sum(b[names(w)] * w); se <- sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); c(est, se) }
fit <- function(x, y, dur, mod, label, extra = "") {
  x$l_dur <- log(x[[dur]]); m <- glm(as.formula(paste(y, "~ z_mx *", mod, "+ z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial)
  V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]; inter <- paste0("z_mx:", mod)
  a <- lin(b, V, c(z_mx = 1)); i <- lin(b, V, setNames(1, inter)); c1 <- lin(b, V, setNames(c(1, 1), c("z_mx", inter)))
  cat(sprintf("%-44s mod=0 %s | interaction %s | mod=1 %s | N=%d events=%d (mod=1: %d rows)\n", label, ci(a[1], a[2]), ci(i[1], i[2]), ci(c1[1], c1[2]), nrow(x), sum(x[[y]]), sum(x[[mod]])))
}
sub <- function(x, y, dur, cond, label, extra = "") {
  x <- x %>% filter(eval(parse(text = cond))); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%-44s larger share %s | N=%d events=%d\n", label, ci(ct["z_mx", 1], ct["z_mx", 2]), nrow(x), sum(x[[y]])))
}
sink("results/distant_pairs_output.txt", split = TRUE)
for (pop in c("embedded pairs (T)", "full panel")) {
  cat(sprintf("\n######## NON-CONTIGUOUS %s ########\n", toupper(pop)))
  dd <- if (pop == "full panel") d else d %>% filter(tri)
  for (oc in list(list("onset", "onset_m", "dur_m", dd %>% filter(m == 0)), list("alliance formation", "form_a", "dur_a", dd %>% filter(a == 0, m == 0)))) {
    cat(sprintf("\n== %s ==\n", oc[[1]]))
    for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE")
      for (mod in c("top5", "ever_major", "same_region")) fit(oc[[4]], oc[[2]], oc[[3]], mod, paste(mod, s), sp)
      sub(oc[[4]], oc[[2]], oc[[3]], "top5 == 0", paste("no top-5 state,", s), sp); sub(oc[[4]], oc[[2]], oc[[3]], "ever_major == 0", paste("no major power,", s), sp); sub(oc[[4]], oc[[2]], oc[[3]], "TRUE", paste("all non-contiguous,", s), sp) } } }
cat("\nShare of non-contiguous embedded onset risk-set rows with a top-5 state:", round(mean(d$top5[d$tri & d$m == 0]), 3), "; ever_major:", round(mean(d$ever_major[d$tri & d$m == 0]), 3), "; same region:", round(mean(d$same_region[d$tri & d$m == 0]), 3), "\n")
sink()
