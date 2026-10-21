# 61: ATOP 5.1 swap and alliance-type splits.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
base <- "raw/atop_download/ex/ATOP 5.1 (.csv)/"
dy <- read_csv(paste0(base, "atop5_1dy.csv"), show_col_types = FALSE) %>% filter(atopally == 1, year <= 2012) %>%
  mutate(lo = pmin(mem1, mem2), hi = pmax(mem1, mem2))
# alliance-level flags; dyad-year table of ATOP alliances restricted to DEFENCE/OFFENCE pledges (primary definition) and, separately, any ATOP alliance
al <- read_csv(paste0(base, "atop5_1a.csv"), show_col_types = FALSE) %>% transmute(atopid, defoff_al = as.integer(defense == 1 | offense == 1), bilat_al = as.integer(bilat == 1), conf = as.integer(coalesce(conwtin, 0) >= 1 | coalesce(medarb, 0) == 1))
ids_all <- dy %>% select(year, lo, hi, atopid1:atopid9) %>% pivot_longer(starts_with("atopid"), values_to = "atopid") %>% filter(!is.na(atopid)) %>% left_join(al, by = "atopid")
dyx_all <- ids_all %>% group_by(year, lo, hi) %>% summarise(defoff = max(defoff_al, na.rm = TRUE), bilat = max(bilat_al, na.rm = TRUE), conf_prov = max(conf, na.rm = TRUE), .groups = "drop")
dyx_def <- ids_all %>% filter(defoff_al == 1) %>% group_by(year, lo, hi) %>% summarise(defoff = 1L, bilat = max(bilat_al, na.rm = TRUE), conf_prov = max(conf, na.rm = TRUE), .groups = "drop")
dyx <- dyx_all
mi <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>% mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>%
  filter(hihost >= 3) %>% transmute(year, lo, hi, sign = -1) %>% distinct()
mk <- function(alli) bind_rows(alli %>% transmute(year, lo, hi, sign = 1), mi) %>% group_by(year, lo, hi) %>%
  summarise(multiplex = n_distinct(sign) > 1, sign = ifelse(any(sign == -1), -1, 1), .groups = "drop") %>% filter(year >= 1816, year <= 2012)
prep <- function(x) x %>% filter(!is.na(dep), !is.na(mut), !is.na(cap)) %>% mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
rhs <- "z_dep + z_mut + z_cap + l_dur + l_unbal + l_bal"
sink("results/atop_checks_output.txt", split = TRUE)
cat("ATOP dyad-years (<=2012), any alliance:", nrow(dyx_all), "| defence/offence alliances:", nrow(dyx_def), "| of which covered by a bilateral alliance:", sum(dyx_def$bilat), "| with conflict-management provision:", sum(dyx_def$conf_prov == 1, na.rm = TRUE), "\n")
cat("\n== (S) Swap COW alliances for ATOP ==\n")
panels <- list()
for (v in c("S1 any ATOP alliance", "S2 ATOP defense/offense pledge")) {
  dd <- if (grepl("S1", v)) dyx_all else dyx_def
  d <- build_panel(mk(dd %>% select(year, lo, hi)), v) %>% left_join(dd, by = c("year", "lo", "hi")) %>% mutate(allied = as.integer(!is.na(defoff)))
  panels[[v]] <- d
  cat(sprintf("%s: allied dispute-years %d of %d; dispute sign changes %d (inside alliance %d)\n", v, sum(d$sign == -1 & d$allied == 1), sum(d$sign == -1), sum(d$sign == -1 & d$flip_i == 1), sum(d$sign == -1 & d$flip_i == 1 & d$allied == 1)))
  for (sg in c(-1, 1)) { x <- prep(d %>% filter(sign == sg)); m <- glm(as.formula(paste("flip_i ~", rhs)), x, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
    cat(sprintf("   %s sign change: larger share %.3f (se %.3f, p=%.4f) | N=%d events=%d\n", ifelse(sg == -1, "disputes ", "alliances"), ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], nobs(m), sum(x$flip_i))) } }
d <- panels[[2]]
cat("\n== (B) bilateral vs multilateral, and (C) conflict-management provisions (S2 panel: ATOP defence/offence alliances) ==\n")
for (spec in list(list(-1, "dispute ends inside alliance", TRUE), list(1, "alliance turns hostile (dispute onset)", FALSE))) {
  x <- prep(d %>% filter(sign == spec[[1]], if (spec[[3]]) allied == 1 else TRUE)) %>% mutate(bil = as.integer(bilat == 1), prov = as.integer(conf_prov == 1))
  for (f in c("bil", "prov")) { xx <- x %>% filter(!is.na(.data[[f]]))
    m <- glm(as.formula(paste("flip_i ~ z_dep *", f, "+ z_mut + z_cap + l_dur + l_unbal + l_bal")), xx, family = binomial); ct <- coeftest(m, vcov = vcovCL(m, cluster = xx$dyad)); nm <- paste0("z_dep:", f)
    cat(sprintf("%-40s %-5s: larger share (flag=0) %.3f (p=%.4f) | x flag %.3f (se %.3f, p=%.4f) | N=%d events=%d, share with flag %.2f\n", spec[[2]], f, ct["z_dep",1], ct["z_dep",4], ct[nm,1], ct[nm,2], ct[nm,4], nobs(m), sum(xx$flip_i), mean(xx[[f]]))) } }
sink()
