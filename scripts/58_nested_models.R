# 58: Nested comparison of exposure and network structure (sign change). Identical observations across models.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
set.seed(20261024)
d0 <- build_panel(base_signed(), "baseline")
auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
specs <- list("Exposure only" = "z_dep + z_mut + z_cap + l_dur",
              "+ unbalanced U" = "z_dep + z_mut + z_cap + l_dur + l_unbal",
              "+ balanced B (full)" = "z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal",
              "Network only" = "z_cap + l_dur + l_unbal + l_bal")
sink("results/nested_models_output.txt", split = TRUE)
for (sg in c(-1, 1)) {
  x <- d0 %>% filter(sign == sg, !is.na(dep), !is.na(mut), !is.na(cap)) %>%
    mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
  cat(sprintf("\n######## %s: N=%d, sign changes=%d, dyads=%d\n", ifelse(sg == -1, "DISPUTES", "ALLIANCES"), nrow(x), sum(x$flip_i), n_distinct(x$dyad)))
  null <- glm(flip_i ~ 1, x, family = binomial); ll0 <- as.numeric(logLik(null))
  dy <- unique(x$dyad); res <- list()
  for (nm in names(specs)) {
    m <- glm(as.formula(paste("flip_i ~", specs[[nm]])), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
    g <- function(v) if (v %in% rownames(ct)) sprintf("%6.3f (%.3f)", ct[v, 1], ct[v, 2]) else "   --"
    cat(sprintf("%-20s larger %s | smaller %s | U %s | B %s | AIC %9.1f BIC %9.1f McF-R2 %.3f\n", nm, g("z_dep"), g("z_mut"), g("l_unbal"), g("l_bal"), AIC(m), BIC(m), 1 - as.numeric(logLik(m)) / ll0))
    res[[nm]] <- list(m = m, ll = as.numeric(logLik(m)))
  }
  cat(sprintf("LR tests (naive, not cluster-robust): U added to exposure-only chi2=%.1f; B added after U chi2=%.1f; exposure added to network-only chi2=%.1f\n",
              2 * (res[[2]]$ll - res[[1]]$ll), 2 * (res[[3]]$ll - res[[2]]$ll), 2 * (res[[3]]$ll - res[[4]]$ll)))
  # CV: dyad-grouped 10-fold, 5 repeats
  cvres <- sapply(names(specs), function(nm) c(ll = 0, auc = 0)); cvll <- matrix(0, 5, 4, dimnames = list(NULL, names(specs))); cvauc <- cvll
  for (rep_ in 1:5) { f <- sample(rep(1:10, length.out = length(dy))); names(f) <- dy; x$fold <- f[x$dyad]
    pr <- sapply(names(specs), function(nm) { p <- rep(NA_real_, nrow(x))
      for (k in 1:10) { tr <- x[x$fold != k, ]; te <- x$fold == k
        m <- glm(as.formula(paste("flip_i ~", specs[[nm]])), tr, family = binomial); p[te] <- predict(m, x[te, ], type = "response") }; p })
    for (nm in names(specs)) { p <- pmin(pmax(pr[, nm], 1e-9), 1 - 1e-9); cvll[rep_, nm] <- mean(x$flip_i * log(p) + (1 - x$flip_i) * log(1 - p)); cvauc[rep_, nm] <- auc(x$flip_i, p) } }
  pn <- mean(x$flip_i); nullll <- mean(x$flip_i * log(pn) + (1 - x$flip_i) * log(1 - pn))
  out <- data.frame(model = names(specs), heldout_ll_per_obs = colMeans(cvll), gain_over_null = colMeans(cvll) - nullll, AUC = colMeans(cvauc), row.names = NULL)
  cat("Out-of-sample (dyad-grouped 10-fold, 5 repeats):\n"); print(out, digits = 6)
}
sink()
