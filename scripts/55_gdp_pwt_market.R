# 55: Trade/GDP with MARKET-rate GDP (Penn World Table 10.01: cgdpo * pl_gdpo), 1950+. Pre-specified before running.
# log(bilateral/GDP) is standardised WITHIN YEAR to remove the common price-level factor.
# Note: bilateral trade is common to both states, so the larger bilateral trade/GDP always belongs to the lower-GDP state.
suppressMessages({ library(sandwich); library(lmtest); library(countrycode) })
library(pwt10)
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
p <- as.data.frame(pwt10.01) %>% transmute(iso = as.character(isocode), year, gdp = cgdpo * pl_gdpo) %>% filter(!is.na(gdp), year >= 1950) %>%
  mutate(ccode = suppressWarnings(countrycode(iso, "iso3c", "cown")), ccode = case_when(iso == "DEU" & year <= 1990 ~ 260, iso == "DEU" ~ 255, iso == "SRB" ~ 345, iso == "YUG" ~ 345, TRUE ~ ccode))
g <- p %>% filter(!is.na(ccode)) %>% distinct(year, ccode, .keep_all = TRUE) %>% select(year, ccode, gdp)
cat("PWT country-years mapped:", nrow(g), " unmapped isocodes:", paste(unique(p$iso[is.na(p$ccode)]), collapse = " "), "\n")
tt <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>% transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), vol = total_dyadic_trade, sl = trade_share_low, sh_ = trade_share_high)
d <- d0 %>% left_join(tt, by = c("year", "lo", "hi")) %>%
  left_join(g %>% rename(lo = ccode, g_lo = gdp), by = c("year", "lo")) %>% left_join(g %>% rename(hi = ccode, g_hi = gdp), by = c("year", "hi")) %>%
  filter(!is.na(g_lo), !is.na(g_hi), !is.na(vol), vol > 0, !is.na(dep), !is.na(mut), !is.na(cap), sl != sh_) %>%
  mutate(lexp_max = log(vol / pmin(g_lo, g_hi)), lexp_min = log(vol / pmax(g_lo, g_hi)),
         weak_gdp = as.integer((sl > sh_ & g_lo < g_hi) | (sh_ > sl & g_hi < g_lo)), lgr = abs(log(g_lo / g_hi)))
zw <- function(x, yr) ave(x, yr, FUN = function(v) if (length(v) > 1 && isTRUE(sd(v) > 0)) (v - mean(v)) / sd(v) else 0)
d <- d %>% mutate(z_emax = zw(lexp_max, year), z_emin = zw(lexp_min, year), z_lvol = zw(log(vol), year), z_lgmin = zw(log(pmin(g_lo, g_hi)), year))
ctrl <- "z_mut + z_cap + l_dur + l_unbal + l_bal"
rep_ <- function(m, x, terms, label) { ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad)); cat(sprintf("%-40s", label))
  for (t in terms) cat(sprintf(" | %s %6.3f (se %.3f, p=%.4f)", t, ct[t,1], ct[t,2], ct[t,4])); cat(sprintf(" | N=%d events=%d\n", nobs(m), sum(m$y))) }
sink("results/gdp_pwt_market_output.txt", split = TRUE)
for (sg in c(-1, 1)) {
  x <- d %>% filter(sign == sg) %>% mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur),
        l_unbal = log1p(n_unbal), l_bal = log1p(n_bal), z_gr = as.numeric(scale(lgr)))
  cat(sprintf("\n######## %s: PWT market-GDP sample 1950+ N=%d, sign changes=%d, share where the more share-exposed state has the lower GDP: %.2f\n", ifelse(sg == -1, "DISPUTES", "ALLIANCES"), nrow(x), sum(x$flip_i), mean(x$weak_gdp)))
  m <- glm(as.formula(paste("flip_i ~ z_dep +", ctrl)), x, family = binomial); rep_(m, x, "z_dep", "G0 larger trade share (baseline, same sample)")
  m <- glm(as.formula(paste("flip_i ~ z_emax + z_emin +", ctrl)), x, family = binomial); rep_(m, x, c("z_emax", "z_emin"), "G1 bilateral trade / GDP (max, min)")
  m <- glm(as.formula(paste("flip_i ~ z_dep + z_emax + z_emin +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "z_emax", "z_emin"), "G1 both measures")
  m <- glm(as.formula(paste("flip_i ~ z_dep + z_lvol + z_lgmin +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "z_lvol", "z_lgmin"), "G1b share + log volume + log min GDP")
  m <- glm(as.formula(paste("flip_i ~ z_dep * weak_gdp +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "weak_gdp", "z_dep:weak_gdp"), "G2 x (exposed state has lower GDP)")
  m <- glm(as.formula(paste("flip_i ~ z_dep * z_gr +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "z_gr", "z_dep:z_gr"), "G3 x |log GDP ratio|")
}
sink()
