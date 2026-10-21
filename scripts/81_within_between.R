# 81: Within-between (Mundlak) decomposition of the larger share for dispute onset, dispute ending, and alliance formation in the embedded population T. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), form_a = as.integer(a1 == 1), tri = (n_same + n_diff) >= 1)
fit <- function(x, y, dur, label, extra = "") {
  x <- x %>% group_by(dyad) %>% mutate(mx_b = mean(z_mx), mx_w = z_mx - mx_b, n_dyad_rows = n()) %>% ungroup(); x$l_dur <- log(x[[dur]])
  m <- glm(as.formula(paste(y, "~ mx_b + mx_w + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)), x, family = binomial)
  V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
  dif <- b["mx_b"] - b["mx_w"]; sed <- sqrt(V["mx_b", "mx_b"] + V["mx_w", "mx_w"] - 2 * V["mx_b", "mx_w"])
  cat(sprintf("%-46s between %s | within %s | difference %s | N=%d events=%d dyads=%d (dyads with >1 row: %d, with variation in the outcome: %d)\n", label, ci(b["mx_b"], sqrt(V["mx_b", "mx_b"])), ci(b["mx_w"], sqrt(V["mx_w", "mx_w"])), ci(dif, sed),
      nrow(x), sum(x[[y]]), n_distinct(x$dyad), n_distinct(x$dyad[x$n_dyad_rows > 1]), x %>% group_by(dyad) %>% filter(n_distinct(.data[[y]]) > 1) %>% pull(dyad) %>% n_distinct()))
}
sink("results/within_between_output.txt", split = TRUE)
cat("Larger share, per SD (global scale): between-dyad mean (mx_b) and within-dyad deviation (mx_w); embedded pairs.\n\n")
pops <- list(list("dispute onset", "onset_m", d %>% filter(tri, m == 0), "dur_m"), list("dispute ending", "end_m", d %>% filter(tri, m == 1), "dur_m"), list("alliance formation", "form_a", d %>% filter(tri, a == 0, m == 0), "dur_a"))
for (p in pops) { fit(p[[3]], p[[2]], p[[4]], paste(p[[1]], "baseline")); fit(p[[3]], p[[2]], p[[4]], paste(p[[1]], "+ volume + decade FE"), "+ z_vol + decade") }
sink()
