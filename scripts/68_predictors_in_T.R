# 68: Predictor comparison and exposure/structure nested comparison within the embedded population T (dispute onset and ending). ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
set.seed(20261027)
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1)
auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
ll_obs <- function(y, p) { p <- pmin(pmax(p, 1e-9), 1 - 1e-9); y * log(p) + (1 - y) * log(1 - p) }
cvp <- function(x, y, specs, reps = 5, K = 10) { dy <- unique(x$dyad); P <- array(NA_real_, c(nrow(x), length(specs), reps), dimnames = list(NULL, names(specs), NULL))
  for (r in 1:reps) { f <- sample(rep(1:K, length.out = length(dy))); names(f) <- dy; fo <- f[x$dyad]
    for (nm in names(specs)) for (kk in 1:K) { m <- suppressWarnings(glm(as.formula(paste(y, "~", specs[[nm]])), x[fo != kk, ], family = binomial)); P[fo == kk, nm, r] <- predict(m, x[fo == kk, ], type = "response") } }
  P }
prep <- function(x, dur) x %>% mutate(z_dur = as.numeric(scale(log(.data[[dur]]))), z_same = as.numeric(scale(l_same)), z_diff = as.numeric(scale(l_diff)), z_mx = as.numeric(scale(mx)), z_mn = as.numeric(scale(mn)), z_cap = as.numeric(scale(cap)))
sink("results/predictors_in_T_output.txt", split = TRUE)
for (spec in list(list("onset_m", d %>% filter(tri, m == 0), "dur_m", "DISPUTE ONSET in T"), list("end_m", d %>% filter(tri, m == 1), "dur_m", "DISPUTE ENDING in T"))) {
  y <- spec[[1]]; x <- prep(spec[[2]], spec[[3]]); lab <- spec[[4]]
  cat(sprintf("\n######## %s: N=%d, events=%d, dyads=%d\n", lab, nrow(x), sum(x[[y]]), n_distinct(x$dyad)))
  preds <- c("z_dur", "z_mx", "z_mn", "z_cap", "z_same", "z_diff")
  m <- glm(as.formula(paste(y, "~", paste(c(preds, "contig"), collapse = " + "))), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  dyv <- unique(x$dyad)
  for (v in preds) { a <- numeric(3); for (r in 1:3) { f <- sample(rep(1:10, length.out = length(dyv))); names(f) <- dyv; fo <- f[x$dyad]; p <- rep(NA_real_, nrow(x))
      for (kk in 1:10) { mm <- glm(as.formula(paste(y, "~", v)), x[fo != kk, ], family = binomial); p[fo == kk] <- predict(mm, x[fo == kk, ], type = "response") }; a[r] <- auc(x[[y]], p) }
    cat(sprintf("%-8s joint coef (per SD) %6.3f (se %.3f, p=%.4f) | single-predictor held-out AUC %.3f\n", v, ct[v, 1], ct[v, 2], ct[v, 4], mean(a))) }
  cat(sprintf("%-8s joint coef %6.3f (se %.3f, p=%.4f)\n", "contig", ct["contig", 1], ct["contig", 2], ct["contig", 4]))
  specs <- list("Exposure + controls" = "z_mx + z_mn + z_cap + contig + z_dur", "+ third-state counts (structure)" = "z_mx + z_mn + z_cap + contig + z_dur + z_same + z_diff", "Structure + controls (no shares)" = "z_cap + contig + z_dur + z_same + z_diff", "Structure + controls, larger share only removed" = "z_mn + z_cap + contig + z_dur + z_same + z_diff")
  P <- cvp(x, y, specs); nullll <- mean(ll_obs(x[[y]], mean(x[[y]])))
  for (nm in names(specs)) { mm <- glm(as.formula(paste(y, "~", specs[[nm]])), x, family = binomial); cc <- coeftest(mm, vcov = vcovCL(mm, cluster = x$dyad))
    cat(sprintf("%-48s larger %s | same %s | diff %s | AIC %9.1f | held-out gain %.5f AUC %.4f\n", nm, if ("z_mx" %in% rownames(cc)) sprintf("%6.3f (%.3f)", cc["z_mx",1], cc["z_mx",2]) else "   --", if ("z_same" %in% rownames(cc)) sprintf("%6.3f (%.3f)", cc["z_same",1], cc["z_same",2]) else "   --", if ("z_diff" %in% rownames(cc)) sprintf("%6.3f (%.3f)", cc["z_diff",1], cc["z_diff",2]) else "   --", AIC(mm),
      mean(sapply(1:5, function(r) mean(ll_obs(x[[y]], P[, nm, r])))) - nullll, mean(sapply(1:5, function(r) auc(x[[y]], P[, nm, r]))))) }
  # incremental held-out value of the larger share (full minus larger-share-removed) with dyad bootstrap CI
  dll <- rowMeans(sapply(1:5, function(r) ll_obs(x[[y]], P[, "+ third-state counts (structure)", r]) - ll_obs(x[[y]], P[, "Structure + controls, larger share only removed", r])))
  idx <- split(seq_len(nrow(x)), x$dyad); bs <- replicate(2000, mean(dll[unlist(idx[sample(unique(x$dyad), replace = TRUE)])]))
  cat(sprintf("Incremental held-out LL per obs of the larger share (structure in): %.5f (95%% CI %.5f to %.5f)\n", mean(dll), quantile(bs, .025), quantile(bs, .975)))
  dls <- rowMeans(sapply(1:5, function(r) ll_obs(x[[y]], P[, "+ third-state counts (structure)", r]) - ll_obs(x[[y]], P[, "Exposure + controls", r]))); bs2 <- replicate(2000, mean(dls[unlist(idx[sample(unique(x$dyad), replace = TRUE)])]))
  cat(sprintf("Incremental held-out LL per obs of structure (exposure in): %.5f (95%% CI %.5f to %.5f)\n", mean(dls), quantile(bs2, .025), quantile(bs2, .975)))
}
sink()
