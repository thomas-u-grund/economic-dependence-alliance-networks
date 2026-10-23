# 10: Resolution of a negative tie (dispute) inside an unbalanced (+,+,-) triad: hazard of the dispute ending. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(fixest) })
suppressMessages(library(cshapes))

src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]
eval(parse(text = src[1:(k - 1)]))
pn <- lapply(split(raw_edges, raw_edges$year), function(e) {
  ids <- sort(unique(c(e$lo, e$hi))); n <- length(ids); ix <- setNames(seq_len(n), ids); p <- e$s == 1L
  P <- matrix(0, n, n); N <- matrix(0, n, n)
  P[cbind(ix[as.character(e$lo[p])], ix[as.character(e$hi[p])])] <- 1; P <- pmax(P, t(P))
  N[cbind(ix[as.character(e$lo[!p])], ix[as.character(e$hi[!p])])] <- 1; N <- pmax(N, t(N))
  list(ix = ix, pp = P %*% P, nn = N %*% N) })
g <- function(yr, lo, hi, w) { o <- rep(0, length(yr)); for (y in unique(yr)) { kk <- which(yr == y); c0 <- pn[[as.character(y)]]; if (is.null(c0)) next
  i <- c0$ix[as.character(lo[kk])]; j <- c0$ix[as.character(hi[kk])]; ok <- !is.na(i) & !is.na(j); o[kk[ok]] <- c0[[w]][cbind(i[ok], j[ok])] }; o }
dm <- d %>% filter(m == 1) %>% mutate(end_m = as.integer(m1 == 0), n_pp = g(year, lo, hi, "pp"), n_nn = g(year, lo, hi, "nn"))
# capital distance (CShapes, 1886+) for dispute dyad-years only
yrs <- sort(unique(dm$year[dm$year >= 1886]))
cap <- bind_rows(lapply(yrs, function(y) { mm <- distmatrix(as.Date(paste0(y, "-01-01")), type = "capdist", useGW = FALSE)
  ids <- as.numeric(rownames(mm)); ix <- which(upper.tri(mm), arr.ind = TRUE)
  tibble(year = y, lo = pmin(ids[ix[, 1]], ids[ix[, 2]]), hi = pmax(ids[ix[, 1]], ids[ix[, 2]]), capdist = mm[ix]) }))
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE); size <- mem %>% group_by(version4id) %>% summarise(n_mem = n_distinct(ccode), .groups = "drop")
alp <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>% filter(defense == 1) %>% left_join(size, by = "version4id") %>%
  mutate(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>% group_by(year, lo, hi) %>% arrange(desc(n_mem)) %>% slice(1) %>% ungroup() %>% transmute(year, lo, hi, pact = ifelse(n_mem >= 5, paste0("p", version4id), "other"))
dm <- dm %>% left_join(cap, by = c("year", "lo", "hi")) %>% left_join(alp, by = c("year", "lo", "hi")) %>%
  mutate(pact = replace_na(pact, "none"), lcap = log(capdist + 1), l_pp = log1p(n_pp), l_nn = log1p(n_nn), l_dur = log(dur_m), asym = mx - mn, z_asym = as.numeric(scale(asym)),
         z_vol = as.numeric(scale(lvol)), z_lcap = as.numeric(scale(lcap)), allied = as.integer(a == 1), era = cut(year, c(1869, 1945, 1989, 2011), labels = c("pre1946", "1946-89", "1990+")))
S <- dm %>% filter(n_pp >= 1, !is.na(capdist))
sink("results/negative_tie_resolution_output.txt", split = TRUE)
cat(sprintf("Sample (n_pp>=1, capdist available): %d dyad-years, %d endings, %d dyads, %d spells; pact groups %d\n", nrow(S), sum(S$end_m), n_distinct(S$dyad), n_distinct(S$dyad, S$rm_), n_distinct(S$pact)))
fit <- function(x, rhs, label, fe = "") {
  f <- as.formula(paste("end_m ~", rhs, fe)); m <- tryCatch(feglm(f, data = x, family = "binomial", cluster = ~dyad), error = function(e) NULL)
  if (is.null(m)) { cat("\n--", label, ": FAILED\n"); return(invisible()) }
  ct <- coeftest(m); cat(sprintf("\n-- %s | N=%d events=%d AIC=%.1f\n", label, nobs(m), sum(x$end_m[x$dyad %in% x$dyad]), AIC(m)))
  keep <- grep("l_pp|l_diff|l_nn|z_mx|z_mn|z_asym|z_cap|z_vol|contig|z_lcap|allied|l_dur", rownames(ct), value = TRUE)
  for (r in keep) cat(sprintf("   %-14s %7.3f (se %.3f, p=%.4f)\n", r, ct[r, 1], ct[r, 2], ct[r, 4])); invisible(m)
}
core <- "l_pp + l_diff + l_nn + l_dur"
fit(S, core, "A pressure only")
fit(S, paste(core, "+ z_mx + z_mn + z_cap"), "B + exposure/dependence + capability")
fit(S, paste(core, "+ z_mx + z_mn + z_cap + z_vol + contig + z_lcap"), "C + volume, contiguity, capital distance")
fit(S, paste(core, "+ z_mx + z_mn + z_cap + z_vol + contig + z_lcap + allied"), "D + allied disputants")
fit(S, paste(core, "+ z_mx + z_mn + z_cap + z_vol + contig + z_lcap + allied"), "E + pact FE", "| pact")
fit(S, paste(core, "+ z_asym + z_cap + z_vol + contig + z_lcap + allied"), "F asymmetry instead of shares")
fit(S, paste(core, "+ z_mx * l_pp + z_mn + z_cap + z_vol + contig + z_lcap + allied"), "G larger share x log1p(common allies)")
fit(S, paste(core, "+ z_mx + z_mn + z_cap + z_vol + contig + z_lcap + allied"), "E2 + era FE (3 eras)", "| era")
cat("\n===== H sensitivity: ALL dispute spells with capdist, n_pp as covariate =====\n")
A <- dm %>% filter(!is.na(capdist)); fit(A, paste(core, "+ z_mx + z_mn + z_cap + z_vol + contig + z_lcap + allied"), "H all disputes (N includes n_pp=0)")
fit(A, paste(core, "+ z_mx + z_mn + z_cap + z_vol + contig + z_lcap + allied"), "H2 all disputes + pact FE", "| pact")
cat("\n===== I stratified (spec D) =====\n")
for (st in list(list("allied disputants", S %>% filter(allied == 1)), list("non-allied disputants", S %>% filter(allied == 0)), list("contiguous", S %>% filter(contig == 1)), list("non-contiguous", S %>% filter(contig == 0)), list("1946+", S %>% filter(year >= 1946))))
  { x <- st[[2]]; if (sum(x$end_m) < 25 || sum(x$end_m == 0) < 15) cat(sprintf("\n-- %s: too few events/non-events (%d/%d)\n", st[[1]], sum(x$end_m), sum(x$end_m == 0))) else
    fit(x, paste(core, "+ z_mx + z_mn + z_cap + z_vol + z_lcap", if (length(unique(x$allied)) > 1) "+ allied" else ""), paste("I", st[[1]])) }
cat("\n===== descriptive: mean predictors by ending =====\n"); print(as.data.frame(S %>% group_by(end_m) %>% summarise(n = n(), n_pp = mean(n_pp), n_diff = mean(n_diff), n_nn = mean(n_nn), larger = mean(mx), smaller = mean(mn), dist_km = mean(capdist), contig = mean(contig), allied = mean(allied), dur = mean(dur_m))))
sink()
