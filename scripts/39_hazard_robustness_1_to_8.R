# 39: Robustness checks 1-8 for the dyad-year hazard (see results/hazard_robustness_1_8_output.txt).
suppressMessages({ library(sandwich); library(lmtest); library(survival); library(nnet) })
source("scripts/38_hazard_panel_lib.R"); set.seed(20261022)
d0 <- build_panel(base_signed(), "baseline")
tr <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), dep_l1 = pmax(trade_share_low, trade_share_high)) %>% mutate(year = year + 1)
d0 <- d0 %>% left_join(tr, by = c("year", "lo", "hi"))

zs <- function(x) as.numeric(scale(x))
prepx <- function(s, need = c("dep","mut","cap")) {
  s <- s[complete.cases(s[, need]), ]
  s %>% mutate(z_dep = zs(dep), z_mut = zs(mut), z_cap = zs(cap), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
}
BASE <- "z_mut + z_cap + l_dur + l_unbal + l_bal"
SG <- c(DISPUTES = -1, ALLIANCES = 1); YY <- c("flip_i", "gone_i")
est <- function(x, y, rhs, vars = "z_dep", link = "logit", label = "", vc = NULL) {
  m <- tryCatch(glm(as.formula(paste(y, "~", rhs)), data = x, family = binomial(link = link)), error = function(e) NULL)
  if (is.null(m)) { cat(sprintf("  %-46s FAILED\n", label)); return(invisible(NULL)) }
  V <- if (is.null(vc)) vcovCL(m, cluster = x$dyad) else vc(m, x)
  ct <- coeftest(m, vcov = V)
  s <- paste(sapply(vars, function(v) sprintf("%s=%6.3f (se %.3f, p=%.4f)", sub("z_", "", v), ct[v,1], ct[v,2], ct[v,4])), collapse = " | ")
  cat(sprintf("  %-46s %s | N=%d ev=%d\n", label, s, nobs(m), sum(x[[y]]))); invisible(ct)
}
hdr <- function(t) cat("\n\n==================== ", t, " ====================\n")
loop <- function(f) for (nm in names(SG)) for (y in YY) { cat(sprintf("\n-- %s / %s\n", nm, y)); tryCatch(f(d0 %>% filter(sign == SG[[nm]]), y, nm), error = function(e) cat("ERROR:", conditionMessage(e), "\n")) }
sink("results/hazard_robustness_1_8_output.txt", split = TRUE)

hdr("CHECK 1: panel structure / repeated spells")
loop(function(s, y, nm) {
  x <- prepx(s)
  est(x, y, paste("z_dep +", BASE), label = "baseline logit")
  est(x %>% mutate(l_prior = log1p(prior_spells)), y, paste("z_dep + l_prior +", BASE), label = "+ log(1+prior spells of dyad)")
  est(x, y, paste("z_dep +", BASE), link = "cloglog", label = "complementary log-log")
  x1 <- prepx(s %>% filter(run == 1)); est(x1, y, paste("z_dep +", BASE), label = "first spell of each dyad only")
  m <- coxph(as.formula(paste("Surv(dur - 1, dur,", y, ") ~ z_dep + z_mut + z_cap + l_unbal + l_bal + cluster(dyad)")), data = x, ties = "efron")
  ct <- summary(m)$coefficients["z_dep", ]
  cat(sprintf("  %-46s dep=%6.3f (robust se %.3f, p=%.4f) | N=%d\n", "Cox, counting process (duration = baseline)", ct["coef"], ct["robust se"], ct["Pr(>|z|)"], m$n))
})

hdr("CHECK 2: state-level dependence (dyadic-robust variance, Cameron-Miller / Aronow et al.)")
dyadic_vc <- function(m, x) {
  U <- estfun(m); B <- vcov(m)
  S1 <- rowsum(U, x$lo); S2 <- rowsum(U, x$hi); st <- union(rownames(S1), rownames(S2))
  S <- matrix(0, length(st), ncol(U), dimnames = list(st, colnames(U)))
  S[rownames(S1), ] <- S[rownames(S1), ] + S1; S[rownames(S2), ] <- S[rownames(S2), ] + S2
  Dm <- rowsum(U, x$dyad)
  B %*% (crossprod(S) - crossprod(Dm)) %*% B
}
loop(function(s, y, nm) {
  x <- prepx(s); r <- paste("z_dep +", BASE)
  est(x, y, r, label = "cluster: dyad")
  est(x, y, r, label = "dyadic-robust (shared state)", vc = dyadic_vc)
  est(x, y, r, label = "cluster: lower-ccode state (naive)", vc = function(m, x) vcovCL(m, cluster = x$lo))
})

hdr("CHECK 3: alliance / dispute definitions")
variants <- list("exclude alliance+dispute overlap (multiplex)" = NULL, "defense pacts only" = list(def = TRUE, h = 3),
                 "hihost >= 4" = list(def = FALSE, h = 4), "defense pacts + hihost >= 4" = list(def = TRUE, h = 4))
d_mult <- d0 %>% filter(!multiplex | is.na(multiplex))
for (vn in names(variants)) {
  cat("\n>>> VARIANT:", vn, "\n")
  dv <- if (is.null(variants[[vn]])) d_mult else build_panel(alt_signed(variants[[vn]]$def, variants[[vn]]$h), vn)
  for (nm in names(SG)) for (y in YY) {
    x <- prepx(dv %>% filter(sign == SG[[nm]])); if (sum(x[[y]]) < 15) { cat(sprintf("  %s/%s: too few events (%d)\n", nm, y, sum(x[[y]]))); next }
    est(x, y, paste("z_dep +", BASE), vars = c("z_dep", "l_unbal", "l_bal"), label = paste(nm, y))
  }
}

hdr("CHECK 4: level vs change in asymmetry")
loop(function(s, y, nm) {
  x <- prepx(s %>% filter(!is.na(dep_l1))) %>% mutate(z_dchg = zs(dep - dep_l1))
  est(x, y, paste("z_dep +", BASE), label = "level only (same sample)")
  est(x, y, paste("z_dchg +", BASE), vars = "z_dchg", label = "change only (dep_t - dep_t-1)")
  est(x, y, paste("z_dep + z_dchg +", BASE), vars = c("z_dep", "z_dchg"), label = "level + change")
})

hdr("CHECK 5: measurement of asymmetry")
loop(function(s, y, nm) {
  x <- prepx(s, c("dep","mut","cap","mx","asym")) %>% mutate(
    z_rank = zs(rank(dep)), z_logdep = zs(log(dep + 0.001)), z_lr = zs(abs(log((mx + 0.001) / (mut + 0.001)))),
    z_relasym = zs(asym / pmax(mx, 1e-9)), z_asym = zs(asym), z_max = zs(mx), z_min = zs(mut), z_wins = zs(pmin(dep, quantile(dep, .99))), z_logasym = zs(log(asym + 0.001)))
  est(x, y, paste("z_dep +", BASE), label = "larger share (baseline)")
  est(x, y, paste("z_asym +", BASE), vars = "z_asym", label = "asymmetry (abs difference), smaller share held fixed")
  est(x, y, paste("z_wins +", BASE), vars = "z_wins", label = "larger share, winsorized 99th")
  est(x, y, paste("z_logdep +", BASE), vars = "z_logdep", label = "log(larger share + .001)")
  est(x, y, paste("z_rank +", BASE), vars = "z_rank", label = "rank of larger share")
  est(x, y, paste("z_lr +", BASE), vars = "z_lr", label = "log-ratio |log(max/min)|")
  est(x, y, paste("z_relasym +", BASE), vars = "z_relasym", label = "relative asymmetry (max-min)/max")
  est(x, y, "z_max + z_min + z_cap + l_dur + l_unbal + l_bal", vars = c("z_max", "z_min"), label = "max share & min share separately")
})

hdr("CHECK 6: functional form of structural pressure")
loop(function(s, y, nm) {
  x <- prepx(s) %>% mutate(unb = n_unbal, bl = n_bal, sh_un = ifelse(any_tri, n_unbal / (n_unbal + n_bal), 0),
    l_tot = log1p(n_unbal + n_bal), q_un = cut(n_unbal, c(-1, 0, 1, 4, 9, Inf)), q_bl = cut(n_bal, c(-1, 0, 1, 4, 9, Inf)))
  B0 <- "z_mut + z_cap + l_dur"
  est(x, y, paste("z_dep + l_unbal + l_bal +", B0), label = "log(1+n) counts (baseline)")
  est(x, y, paste("z_dep + unb + bl +", B0), label = "raw counts")
  est(x, y, paste("z_dep + sqrt(unb) + sqrt(bl) +", B0), label = "sqrt counts")
  est(x, y, paste("z_dep + sh_un + l_tot + any_tri +", B0), label = "share unbalanced + log total + any-triad")
  est(x, y, paste("z_dep + q_un + q_bl +", B0), label = "categorical bins (0,1,2-4,5-9,10+)")
})

hdr("CHECK 7: competing risks")
for (nm in names(SG)) {
  x <- prepx(d0 %>% filter(sign == SG[[nm]]))
  cat(sprintf("\n-- %s: stacked cause-specific logit, contrast on asymmetry (sign change minus dissolution)\n", nm))
  L <- bind_rows(x %>% mutate(ev = flip_i, typ = 1), x %>% mutate(ev = gone_i, typ = 0)) %>% mutate(rid = rep(seq_len(nrow(x)), 2))
  m <- glm(ev ~ typ * (z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal), data = L, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = L$dyad))
  cat(sprintf("  dissolution: dep=%6.3f (p=%.4f) | sign change - dissolution contrast: %6.3f (se %.3f, p=%.4f)\n",
              ct["z_dep",1], ct["z_dep",4], ct["typ:z_dep",1], ct["typ:z_dep",2], ct["typ:z_dep",4]))
  x$out <- factor(ifelse(x$flip_i == 1, "flip", ifelse(x$gone_i == 1, "gone", "stay")), levels = c("stay", "flip", "gone"))
  mm <- multinom(out ~ z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal, data = x, trace = FALSE)
  cf <- summary(mm)$coefficients; sx <- summary(mm)$standard.errors
  cat(sprintf("  multinomial (ref = tie continues): dep -> flip %6.3f (naive se %.3f) | dep -> dissolve %6.3f (naive se %.3f)\n",
              cf["flip","z_dep"], sx["flip","z_dep"], cf["gone","z_dep"], sx["gone","z_dep"]))
}

hdr("CHECK 8: outcome windows (any event within k years)")
ev <- d0 %>% select(dyad, year, flip_i, gone_i)
win <- function(k) { out <- d0 %>% select(dyad, year)
  f <- rep(0L, nrow(out)); g <- rep(0L, nrow(out))
  for (j in 0:(k - 1)) { jj <- out %>% mutate(year = year + j) %>% left_join(ev, by = c("dyad", "year")); f <- pmax(f, replace_na(jj$flip_i, 0L)); g <- pmax(g, replace_na(jj$gone_i, 0L)) }
  tibble(!!paste0("flip_w", k) := f, !!paste0("gone_w", k) := g) }
dw <- bind_cols(d0, win(1), win(3), win(5)) %>% filter(year <= 2011 - 4)   # complete 5-year windows only
for (nm in names(SG)) for (k in c(1, 3, 5)) for (yv in c("flip", "gone")) {
  y <- paste0(yv, "_w", k); x <- prepx(dw %>% filter(sign == SG[[nm]]))
  est(x, y, paste("z_dep +", BASE), vars = c("z_dep", "l_unbal", "l_bal"), label = sprintf("%s %s within %d y", nm, yv, k))
}
sink()
