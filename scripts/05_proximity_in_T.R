# 05: Proximity / interaction / pact diagnostic for dispute onset and ending within the embedded population T. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(fixest) })
suppressMessages(library(cshapes))

src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1) %>% filter(year >= 1886, year <= 2011)
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
sink("results/proximity_in_T_output.txt", split = TRUE)
x1 <- d %>% filter(tri, m == 0); run(x1, "onset_m", "(i) DISPUTE ONSET in T")
run(x1 %>% filter(a == 1), "onset_m", "(i-a) DISPUTE ONSET in T, allied pairs")
run(x1 %>% filter(a == 0), "onset_m", "(i-n) DISPUTE ONSET in T, non-allied pairs")
run(x1 %>% filter(land == 1), "onset_m", "(i-b) DISPUTE ONSET in T, land-border pairs")
x2 <- d %>% filter(tri, m == 1); run(x2, "end_m", "(ii) DISPUTE ENDING in T")
run(x2 %>% filter(a == 1), "end_m", "(ii-a) DISPUTE ENDING in T, allied pairs")
sink()
