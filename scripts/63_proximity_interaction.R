# 63: Does the exposure association survive proximity, prior interaction and pact fixed effects? (dispute onset among allies)
suppressMessages({ library(sandwich); library(lmtest); library(fixest) })
suppressMessages(library(cshapes))
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
# --- capital distance by year ---
yrs <- sort(unique(d0$year[d0$year >= 1886]))
cap <- bind_rows(lapply(yrs, function(y) { m <- distmatrix(as.Date(paste0(y, "-01-01")), type = "capdist", useGW = FALSE)
  ids <- as.numeric(rownames(m)); ix <- which(upper.tri(m), arr.ind = TRUE)
  tibble(year = y, lo = pmin(ids[ix[,1]], ids[ix[,2]]), hi = pmax(ids[ix[,1]], ids[ix[,2]]), capdist = m[ix]) }))
# --- contiguity by type ---
ct <- read_csv("raw/contiguity/DirectContiguity320/contdird.csv", show_col_types = FALSE) %>% transmute(year, lo = pmin(state1no, state2no), hi = pmax(state1no, state2no), conttype) %>%
  group_by(year, lo, hi) %>% summarise(land = as.integer(any(conttype == 1)), water = as.integer(any(conttype %in% 2:5)), .groups = "drop")
# --- prior militarised activity (any MID with hihost>=2) ---
mi <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>% filter(hihost >= 2) %>% distinct(year, lo, hi) %>% mutate(act = 1L)
pri <- mi %>% group_by(lo, hi) %>% arrange(year) %>% mutate(first_mid = min(year)) %>% ungroup()
first <- pri %>% distinct(lo, hi, first_mid)
prior10 <- bind_rows(lapply(1:10, function(k) mi %>% transmute(year = year + k, lo, hi, act))) %>% group_by(year, lo, hi) %>% summarise(prior10 = sum(act), .groups = "drop")
# --- pact id: largest covering alliance with >=5 members ---
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE); size <- mem %>% group_by(version4id) %>% summarise(n_mem = n_distinct(ccode), .groups = "drop")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>% filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>%
  left_join(size, by = "version4id") %>% mutate(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>% group_by(year, lo, hi) %>% arrange(desc(n_mem)) %>% slice(1) %>% ungroup() %>%
  transmute(year, lo, hi, allied = 1L, pact = ifelse(n_mem >= 5, paste0("p", version4id), "other"))
tt <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>% transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), lvol = log1p(total_dyadic_trade))
d <- d0 %>% filter(year >= 1886) %>% left_join(cap, by = c("year", "lo", "hi")) %>% left_join(ct, by = c("year", "lo", "hi")) %>% left_join(prior10, by = c("year", "lo", "hi")) %>%
  left_join(first, by = c("lo", "hi")) %>% left_join(al, by = c("year", "lo", "hi")) %>% left_join(tt, by = c("year", "lo", "hi")) %>%
  mutate(land = replace_na(land, 0L), water = replace_na(water, 0L), prior10 = replace_na(prior10, 0L), any_prior = as.integer(!is.na(first_mid) & first_mid < year), allied = replace_na(allied, 0L), pact = replace_na(pact, "none"),
         lcap = log(capdist + 1), l_prior10 = log1p(prior10))
z <- function(x) as.numeric(scale(x))
run <- function(x, label) {
  x <- x %>% filter(!is.na(dep), !is.na(mut), !is.na(cap), !is.na(capdist), !is.na(lvol)) %>%
    mutate(z_dep = z(dep), z_mut = z(mut), z_cap = z(cap), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal), z_lcap = z(lcap), z_vol = z(lvol), z_pri = z(l_prior10))
  cat(sprintf("\n######## %s: N=%d, events=%d, dyads=%d, pact groups=%d, corr(larger share, log distance)=%.2f\n", label, nrow(x), sum(x$flip_i), n_distinct(x$dyad), n_distinct(x$pact), cor(x$z_dep, x$z_lcap)))
  base <- "z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal"
  specs <- list("A main" = base, "B1 + land border, water contiguity, log capital distance" = paste(base, "+ land + water + z_lcap"),
                "B2 + log trade volume" = paste(base, "+ land + water + z_lcap + z_vol"), "C + prior militarised activity" = paste(base, "+ land + water + z_lcap + z_vol + z_pri + any_prior"),
                "D + pact fixed effects" = paste(base, "+ land + water + z_lcap + z_vol + z_pri + any_prior | pact"))
  for (nm in names(specs)) {
    f <- as.formula(paste("flip_i ~", specs[[nm]]))
    m <- tryCatch(feglm(f, data = x, family = "binomial", cluster = ~dyad), error = function(e) NULL)
    if (is.null(m)) { cat(sprintf("  %-58s (failed)\n", nm)); next }
    co <- coeftable(m); cat(sprintf("  %-58s larger share %6.3f (se %.3f, p=%.4f) | N used=%d\n", nm, co["z_dep",1], co["z_dep",2], co["z_dep",4], m$nobs)) }
  invisible(x)
}
sink("results/proximity_interaction_output.txt", split = TRUE)
x1 <- d %>% filter(sign == 1)
run(x1, "(i) DISPUTE BEGINS INSIDE ALLIANCE (alliance-years)")
run(x1 %>% filter(land == 1), "(i-b) land-border dyads only")
run(x1 %>% filter(pact != "other"), "(i-c) within large pacts only (largest covering alliance >=5 members)")
x2 <- d %>% filter(sign == -1, allied == 1)
run(x2, "(ii) DISPUTE ENDS INSIDE ALLIANCE (allied dispute-years)")
sink()
