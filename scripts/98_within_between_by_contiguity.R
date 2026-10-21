# 98: FULL PANEL: within-between decomposition of the larger share separately for contiguous and non-contiguous dyads. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(fixest) })
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), form_a = as.integer(a1 == 1))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
wb <- function(x, y, dur, label, extra = "") {
  x <- x %>% group_by(dyad) %>% mutate(mx_b = mean(z_mx), mx_w = z_mx - mx_b) %>% ungroup(); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ mx_b + mx_w + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  dif <- b["mx_b"] - b["mx_w"]; sed <- sqrt(V["mx_b", "mx_b"] + V["mx_w", "mx_w"] - 2 * V["mx_b", "mx_w"])
  cat(sprintf("%-46s between %s | within %s | difference %s | N=%d events=%d dyads=%d\n", label, ci(b["mx_b"], sqrt(V["mx_b", "mx_b"])), ci(b["mx_w"], sqrt(V["mx_w", "mx_w"])), ci(dif, sed), nrow(x), sum(x[[y]]), n_distinct(x$dyad))) }
wbi <- function(x, y, dur, label, extra = "") {
  x <- x %>% group_by(dyad) %>% mutate(mx_b = mean(z_mx), mx_w = z_mx - mx_b) %>% ungroup(); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ (mx_b + mx_w) * contig + z_mn + z_cap + l_dur + l_same + l_diff", extra)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  cat(sprintf("%-46s contiguity difference in between %s | in within %s\n", label, ci(b["mx_b:contig"], sqrt(V["mx_b:contig", "mx_b:contig"])), ci(b["mx_w:contig"], sqrt(V["mx_w:contig", "mx_w:contig"])))) }
sink("results/within_between_by_contiguity_output.txt", split = TRUE)
pops <- list(list("dispute onset", "onset_m", "dur_m", d %>% filter(m == 0)), list("alliance formation", "form_a", "dur_a", d %>% filter(a == 0, m == 0)))
for (p in pops) for (sp in c("", "+ z_vol + decade")) { s <- ifelse(sp == "", "baseline", "+ vol + decade FE")
  cat("\n")
  wb(p[[4]] %>% filter(contig == 0), p[[2]], p[[3]], paste(p[[1]], "non-contiguous,", s), sp)
  wb(p[[4]] %>% filter(contig == 1), p[[2]], p[[3]], paste(p[[1]], "contiguous,", s), sp)
  wbi(p[[4]], p[[2]], p[[3]], paste(p[[1]], "interaction model,", s), sp) }
cat("\nDyad fixed-effects logit among non-contiguous dyads, larger share (per SD)\n")
for (p in pops) for (sp in c("", "+ z_vol")) { x <- p[[4]] %>% filter(contig == 0); x$l_dur <- log(x[[p[[3]]]])
  m <- feglm(as.formula(paste(p[[2]], "~ z_mx + z_mn + z_cap + l_dur + l_same + l_diff", sp, "| dyad")), x, family = "binomial", cluster = ~dyad)
  cat(sprintf("%-20s %-10s %s (p=%.3f) | dyads used %d, obs %d\n", p[[1]], ifelse(sp == "", "baseline", sp), ci(coef(m)["z_mx"], se(m)["z_mx"]), pvalue(m)["z_mx"], length(unique(x$dyad[x[[p[[2]]]] == 1])), nobs(m))) }
sink()
