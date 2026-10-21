# 59: Does the network (U, B) information survive conditioning on concurrent alliance status? And incremental value of the larger share with CI.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
set.seed(20261025)
d0 <- build_panel(base_signed(), "baseline")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>% transmute(year, lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), allied = 1L) %>% distinct()
auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
prep <- function(sg) d0 %>% filter(sign == sg, !is.na(dep), !is.na(mut), !is.na(cap)) %>% left_join(al, by = c("year", "lo", "hi")) %>%
  mutate(allied = replace_na(allied, 0L), z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
cvpred <- function(x, specs, reps = 5, K = 10) {
  dy <- unique(x$dyad); P <- array(NA_real_, c(nrow(x), length(specs), reps), dimnames = list(NULL, names(specs), NULL))
  for (r in 1:reps) { f <- sample(rep(1:K, length.out = length(dy))); names(f) <- dy; fo <- f[x$dyad]
    for (nm in names(specs)) for (k in 1:K) { tr <- x[fo != k, ]; te <- fo == k
      m <- suppressWarnings(glm(as.formula(paste("flip_i ~", specs[[nm]])), tr, family = binomial)); P[te, nm, r] <- predict(m, x[te, ], type = "response") } }
  P }
ll_obs <- function(y, p) { p <- pmin(pmax(p, 1e-9), 1 - 1e-9); y * log(p) + (1 - y) * log(1 - p) }
sink("results/structure_vs_alliance_status_output.txt", split = TRUE)
# ---------------- D1: disputes ----------------
x <- prep(-1)
cat(sprintf("\n######## D1 DISPUTES: N=%d, sign changes=%d, allied dispute-years=%d (sign changes among them %d)\n", nrow(x), sum(x$flip_i), sum(x$allied), sum(x$flip_i[x$allied == 1])))
specs <- list("(e) allied + cap + dur" = "allied + z_cap + l_dur",
              "(a) allied + exposure" = "allied + z_dep + z_mut + z_cap + l_dur",
              "(b) + U" = "allied + z_dep + z_mut + z_cap + l_dur + l_unbal",
              "(c) + U + B (full)" = "allied + z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal",
              "(d) full minus larger share" = "allied + z_mut + z_cap + l_dur + l_unbal + l_bal",
              "(f) allied + network only" = "allied + z_cap + l_dur + l_unbal + l_bal")
null_ll <- mean(ll_obs(x$flip_i, mean(x$flip_i)))
P <- cvpred(x, specs)
for (nm in names(specs)) {
  m <- glm(as.formula(paste("flip_i ~", specs[[nm]])), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  g <- function(v) if (v %in% rownames(ct)) sprintf("%6.3f (%.3f)", ct[v, 1], ct[v, 2]) else "   --"
  cvll <- mean(sapply(1:5, function(r) mean(ll_obs(x$flip_i, P[, nm, r])))); cva <- mean(sapply(1:5, function(r) auc(x$flip_i, P[, nm, r])))
  cat(sprintf("%-28s allied %s | larger %s | U %s | B %s | AIC %7.1f | held-out gain %.4f AUC %.3f\n", nm, g("allied"), g("z_dep"), g("l_unbal"), g("l_bal"), AIC(m), cvll - null_ll, cva))
}
cat("\n-- separately by alliance status (U, B, exposure; same controls) --\n")
for (a in c(1, 0)) { s <- x[x$allied == a, ]
  m <- glm(flip_i ~ z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal, s, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = s$dyad))
  cat(sprintf("allied=%d: N=%d events=%d | larger %.3f (%.3f) | U %.3f (%.3f) | B %.3f (%.3f)\n", a, nrow(s), sum(s$flip_i), ct["z_dep",1], ct["z_dep",2], ct["l_unbal",1], ct["l_unbal",2], ct["l_bal",1], ct["l_bal",2])) }
# ---------------- D2: incremental value with bootstrap CI ----------------
boot_diff <- function(x, full, red, label) {
  sp <- list(full = full, red = red); P <- cvpred(x, sp)
  d <- rowMeans(sapply(1:5, function(r) ll_obs(x$flip_i, P[, "full", r]) - ll_obs(x$flip_i, P[, "red", r])))
  dy <- unique(x$dyad); idx <- split(seq_len(nrow(x)), x$dyad)
  bs <- replicate(2000, { s <- sample(dy, replace = TRUE); mean(d[unlist(idx[s])]) })
  da <- mean(sapply(1:5, function(r) auc(x$flip_i, P[, "full", r]) - auc(x$flip_i, P[, "red", r])))
  cat(sprintf("%-60s mean held-out LL gain per obs %.5f (95%% CI %.5f to %.5f) | AUC gain %.4f\n", label, mean(d), quantile(bs, .025), quantile(bs, .975), da)) }
cat("\n======== D2: incremental held-out value ========\n")
fullD <- "allied + z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal"
boot_diff(x, fullD, "allied + z_mut + z_cap + l_dur + l_unbal + l_bal", "Disputes: larger share added to (allied + smaller share + controls + U + B)")
boot_diff(x, fullD, "allied + z_dep + z_mut + z_cap + l_dur", "Disputes: U and B added to (allied + exposure + controls)")
boot_diff(x, "allied + z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal", "z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal", "Disputes: allied indicator added to (exposure + controls + U + B)")
y <- prep(1)
fullA <- "z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal"
boot_diff(y, fullA, "z_mut + z_cap + l_dur + l_unbal + l_bal", "Alliances: larger share added to (smaller share + controls + U + B)")
boot_diff(y, fullA, "z_dep + z_mut + z_cap + l_dur", "Alliances: U and B added to (exposure + controls)")
sink()
