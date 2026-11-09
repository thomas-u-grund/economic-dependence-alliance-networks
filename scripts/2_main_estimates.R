# 2_main_estimates: Section 5.2 and Table 3: transitions on the full panel, the controls path, and actor-overlap standard errors (Appendix A, Tables 1 and 2).
# Run from the repository root: Rscript scripts/2_main_estimates.R (or run everything with Rscript run_all.R).
# Each block below is self-contained and writes the results/ file named in its sink() call.
source("scripts/setup.R")

# ---- Full-panel versions of the main results, and lagged definition of embedding.
local({
cl <- get_cnt(d$year - 1, d$lo, d$hi); d$n_same_l <- cl$s; d$n_diff_l <- cl$d; d$l_same_l <- log1p(d$n_same_l); d$l_diff_l <- log1p(d$n_diff_l)
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), form_a = as.integer(a1 == 1), end_a = as.integer(a1 == 0), tri = (n_same + n_diff) >= 1, tri_l = (n_same_l + n_diff_l) >= 1)
mid <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb), fatlev = ifelse(fatlev < 0, NA, fatlev)) %>%
  group_by(year, lo, hi) %>% summarise(h = max(hihost), fat = suppressWarnings(max(fatlev, na.rm = TRUE)), .groups = "drop") %>% mutate(fat = ifelse(is.infinite(fat), NA, fat))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
est1 <- function(x, y, dur, label, sp = "", same = "l_same", diff = "l_diff", inter = NULL) {
  x$l_dur <- log(x[[dur]]); rhs <- if (is.null(inter)) "z_mx" else paste("z_mx *", inter)
  m <- glm(as.formula(paste(y, "~", rhs, "+ z_mn + z_cap + contig + l_dur +", same, "+", diff, sp)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
  out <- sprintf("%-46s larger share %s", label, ci(b["z_mx"], sqrt(V["z_mx", "z_mx"])))
  if (!is.null(inter)) { it <- paste0("z_mx:", inter); out <- paste0(out, sprintf(" | interaction %s", ci(b[it], sqrt(V[it, it])))) }
  cat(sprintf("%s | N=%d events=%d dyads=%d\n", out, nrow(x), sum(x[[y]]), n_distinct(x$dyad[x[[y]] == 1])))
}
pops <- function(dd) list(list("dispute onset", "onset_m", "dur_m", dd %>% filter(m == 0)), list("dispute ending", "end_m", "dur_m", dd %>% filter(m == 1)), list("alliance formation", "form_a", "dur_a", dd %>% filter(a == 0, m == 0)), list("alliance ending", "end_a", "dur_a", dd %>% filter(a == 1, m == 0)))
sink("results/population_selection_output.txt", split = TRUE)
cat("######## (a) FULL PANEL, contemporaneous partner counts as controls ########\n")
for (p in pops(d)) for (sp in c("", "+ z_vol + decade")) est1(p[[4]], p[[2]], p[[3]], paste(p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE")), sp)
cat("\n-- full panel, larger share x contiguity (onset, formation)\n")
for (p in pops(d)[c(1, 3)]) for (sp in c("", "+ z_vol + decade")) est1(p[[4]], p[[2]], p[[3]], paste(p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE")), sp, inter = "contig")
cat("\n-- full panel: severity of onsets (risk set: no MID at t)\n")
dm <- d %>% left_join(mid, by = c("year", "lo", "hi")) %>% mutate(h = replace_na(h, 0)) %>% filter(h == 0) %>% left_join(mid %>% transmute(year = year - 1, lo, hi, h1 = h, fat1 = fat), by = c("year", "lo", "hi")) %>%
  mutate(h1 = replace_na(h1, 0), fat1 = replace_na(fat1, 0), y_dis = as.integer(h1 == 3), y_force = as.integer(h1 >= 4), y_fat = as.integer(fat1 >= 1), y_any = as.integer(h1 >= 2))
dm <- dm %>% arrange(dyad, year) %>% group_by(dyad) %>% mutate(run = cumsum(c(TRUE, diff(year) != 1))) %>% group_by(dyad, run) %>% mutate(dur_q = row_number()) %>% ungroup()
for (y in c("y_dis", "y_force", "y_fat", "y_any")) for (sp in c("", "+ z_vol + decade")) est1(dm, y, "dur_q", paste("severity", y, ifelse(sp == "", "baseline", "+ vol + decade FE")), sp)
cat("\n######## (b) LAGGED EMBEDDING: embedded pair and partner counts from the network of year t-1 ########\n")
for (p in pops(d %>% filter(tri_l))) for (sp in c("", "+ z_vol + decade")) est1(p[[4]], p[[2]], p[[3]], paste("lagged-embedded", p[[1]], ifelse(sp == "", "baseline", "+ vol + decade FE")), sp, same = "l_same_l", diff = "l_diff_l")
cat(sprintf("\nShare of contemporaneously embedded dyad-years that are also embedded at t-1: %.3f; of t-1 embedded that are embedded at t: %.3f\n", mean(d$tri_l[d$tri]), mean(d$tri[d$tri_l])))
sink()
})

# ---- Proximity / interaction / pact diagnostic for dispute onset and ending, FULL PANEL (all dyad-years, 1886+).
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0)) %>% filter(year >= 1886, year <= 2011)
yrs <- sort(unique(d$year))
cap <- bind_rows(lapply(yrs, function(y) { m <- distmatrix(as.Date(paste0(y, "-01-01")), type = "capdist", useGW = FALSE); ids <- as.numeric(rownames(m)); ix <- which(upper.tri(m), arr.ind = TRUE)
  tibble(year = y, lo = pmin(ids[ix[,1]], ids[ix[,2]]), hi = pmax(ids[ix[,1]], ids[ix[,2]]), capdist = m[ix]) }))
ct <- read_csv("raw/contiguity/DirectContiguity320/contdird.csv", show_col_types = FALSE) %>% transmute(year, lo = pmin(state1no, state2no), hi = pmax(state1no, state2no), conttype) %>%
  group_by(year, lo, hi) %>% summarise(land = as.integer(any(conttype == 1)), water = as.integer(any(conttype %in% 2:5)), .groups = "drop")
mi2 <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>% filter(hihost >= 2) %>% distinct(year, lo, hi) %>% mutate(act = 1L)
first <- mi2 %>% group_by(lo, hi) %>% summarise(first_mid = min(year), .groups = "drop")
prior10 <- bind_rows(lapply(1:10, function(k) mi2 %>% transmute(year = year + k, lo, hi, act))) %>% group_by(year, lo, hi) %>% summarise(prior10 = sum(act), .groups = "drop")
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE); size <- mem %>% group_by(version4id) %>% summarise(n_mem = n_distinct(ccode), .groups = "drop")
pact <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>% filter(defense == 1) %>% left_join(size, by = "version4id") %>% mutate(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>%
  group_by(year, lo, hi) %>% arrange(desc(n_mem)) %>% slice(1) %>% ungroup() %>% transmute(year, lo, hi, pact = ifelse(n_mem >= 5, paste0("p", version4id), "other"))
d <- d %>% left_join(cap, by = c("year", "lo", "hi")) %>% left_join(ct, by = c("year", "lo", "hi")) %>% left_join(prior10, by = c("year", "lo", "hi")) %>% left_join(first, by = c("lo", "hi")) %>% left_join(pact, by = c("year", "lo", "hi")) %>%
  mutate(land = replace_na(land, 0L), water = replace_na(water, 0L), prior10 = replace_na(prior10, 0L), any_prior = as.integer(!is.na(first_mid) & first_mid < year), pact = replace_na(pact, "none"), lcap = log(capdist + 1), z_lcap = as.numeric(scale(lcap)), z_pri = as.numeric(scale(log1p(prior10))))
run <- function(x, y, label, dur = "dur_m") {
  x <- x %>% filter(!is.na(capdist)) %>% mutate(l_dur = log(.data[[dur]]), zz_vol = z_vol)
  cat(sprintf("\n######## %s: N=%d, events=%d, dyads=%d, corr(larger share, log distance)=%.2f\n", label, nrow(x), sum(x[[y]]), n_distinct(x$dyad), cor(x$z_mx, x$z_lcap)))
  base <- "z_mx + z_mn + z_cap + l_dur + l_same + l_diff"
  specs <- list("A main (contiguity)" = paste(base, "+ contig"), "B1 + land, water contiguity, log capital distance" = paste(base, "+ land + water + z_lcap"), "B2 + log trade volume" = paste(base, "+ land + water + z_lcap + zz_vol"),
                "C + prior militarised activity" = paste(base, "+ land + water + z_lcap + zz_vol + z_pri + any_prior"), "D + pact fixed effects" = paste(base, "+ land + water + z_lcap + zz_vol + z_pri + any_prior | pact"))
  for (nm in names(specs)) { m <- tryCatch(feglm(as.formula(paste(y, "~", specs[[nm]])), data = x, family = "binomial", cluster = ~dyad), error = function(e) NULL)
    if (is.null(m)) { cat(sprintf("  %-55s (failed)\n", nm)); next }; co <- coeftable(m); cat(sprintf("  %-55s larger share %6.3f (se %.3f, p=%.4f) | N used=%d\n", nm, co["z_mx",1], co["z_mx",2], co["z_mx",4], m$nobs)) } }
sink("results/fullpanel_controls_output.txt", split = TRUE)
x1 <- d %>% filter(m == 0); run(x1, "onset_m", "(i) DISPUTE ONSET, full panel")
run(x1 %>% filter(a == 1), "onset_m", "(i-a) DISPUTE ONSET, allied pairs")
run(x1 %>% filter(a == 0), "onset_m", "(i-n) DISPUTE ONSET, non-allied pairs")
run(x1 %>% filter(land == 1), "onset_m", "(i-b) DISPUTE ONSET, land-border pairs")
run(x1 %>% filter(contig == 0), "onset_m", "(i-c) DISPUTE ONSET, non-contiguous pairs")
x2 <- d %>% filter(m == 1); run(x2, "end_m", "(ii) DISPUTE ENDING, full panel")
run(x2 %>% filter(a == 1), "end_m", "(ii-a) DISPUTE ENDING, allied pairs")
sink()
})

# ---- FULL PANEL: dyad-clustered vs actor-overlap (dyadic-robust) SEs for the larger share.
local({
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), form_a = as.integer(a1 == 1), end_a = as.integer(a1 == 0))
robust_se <- function(x, y, extra, dur) {
  x$l_dur <- log(x[[dur]]); f <- as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)); m <- glm(f, x, family = binomial); X <- model.matrix(m); keep <- !is.na(coef(m)); X <- X[, keep, drop = FALSE]
  fit <- glm.fit(X, x[[y]], family = binomial()); mu <- fit$fitted.values; U <- X * (x[[y]] - mu); n <- nrow(X); A <- solve(crossprod(X, X * (mu * (1 - mu))))
  Sd <- rowsum(U, x$dyad); Gd <- nrow(Sd); Vd <- A %*% crossprod(Sd) %*% A * (Gd / (Gd - 1)); chk <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad)); j <- which(colnames(X) == "z_mx")
  stopifnot(abs(sqrt(Vd[j, j]) - chk["z_mx", 2]) / chk["z_mx", 2] < 1e-3)
  st <- rbind(data.frame(row = seq_len(n), s = x$lo), data.frame(row = seq_len(n), s = x$hi)); Ss <- rowsum(U[st$row, , drop = FALSE], group = st$s); Gs <- length(unique(st$s)); Vs <- A %*% (crossprod(Ss) - crossprod(Sd)) %*% A * (Gs / (Gs - 1))
  list(est = unname(coef(fit)["z_mx"]), se_dyad = sqrt(Vd[j, j]), se_dyadic = sqrt(max(Vs[j, j], 0)), N = n, ev = sum(x[[y]])) }
pops <- list(list("dispute onset", "onset_m", d %>% filter(m == 0), "dur_m"), list("dispute ending", "end_m", d %>% filter(m == 1), "dur_m"), list("alliance formation", "form_a", d %>% filter(a == 0, m == 0), "dur_a"), list("alliance ending", "end_a", d %>% filter(a == 1, m == 0), "dur_a"))
res <- list(); for (p in pops) for (sp in c("baseline", "+ volume + decade FE")) { r <- robust_se(p[[3]], p[[2]], if (sp == "baseline") "" else "+ z_vol + decade", p[[4]])
  res[[length(res) + 1]] <- data.frame(outcome = p[[1]], spec = sp, N = r$N, events = r$ev, est = r$est, se_dyad = r$se_dyad, se_dyadic = r$se_dyadic, ratio = r$se_dyadic / r$se_dyad, p_dyadic = 2 * pnorm(-abs(r$est / r$se_dyadic)), ci_lo_dyadic = r$est - 1.96 * r$se_dyadic, ci_hi_dyadic = r$est + 1.96 * r$se_dyadic) }
tab <- bind_rows(res) %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))
sink("results/fullpanel_shared_se_output.txt", split = TRUE); cat("Larger share, FULL PANEL: dyad-clustered vs actor-overlap SE\n"); print(as.data.frame(tab), row.names = FALSE); sink(); write_csv(tab, "results/fullpanel_shared_se.csv")
})
