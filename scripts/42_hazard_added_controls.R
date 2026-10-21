# 42: Add contiguity, economic-size asymmetry, and regime type (Polity2) to the dyad-year hazard.
suppressMessages({ library(sandwich); library(lmtest); library(readxl) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
contig <- read_csv("raw/contiguity/DirectContiguity320/contdird.csv", show_col_types = FALSE) %>%
  transmute(year, lo = pmin(state1no, state2no), hi = pmax(state1no, state2no), contig = 1L) %>% distinct()
nat <- read_csv("raw/trade/COW_Trade_4.0/National_COW_4.0.csv", show_col_types = FALSE) %>%
  mutate(tt = ifelse(imports < 0 | exports < 0, NA, imports + exports)) %>% filter(!is.na(tt), tt > 0) %>% select(ccode, year, tt)
pol <- read_excel("raw/polity/p5v2018.xls") %>% transmute(ccode, year, polity2 = ifelse(polity2 < -10, NA, polity2)) %>% filter(!is.na(polity2))
d <- d0 %>% left_join(contig, by = c("year", "lo", "hi")) %>% mutate(contig = replace_na(contig, 0L)) %>%
  left_join(nat %>% rename(lo = ccode, tt_lo = tt), by = c("year", "lo")) %>% left_join(nat %>% rename(hi = ccode, tt_hi = tt), by = c("year", "hi")) %>%
  left_join(pol %>% rename(lo = ccode, p_lo = polity2), by = c("year", "lo")) %>% left_join(pol %>% rename(hi = ccode, p_hi = polity2), by = c("year", "hi")) %>%
  mutate(size_asym = abs(log(tt_lo) - log(tt_hi)), min_polity = pmin(p_lo, p_hi), polity_gap = abs(p_lo - p_hi),
         joint_dem = as.integer(p_lo >= 6 & p_hi >= 6))
cat("contiguous share:", mean(d$contig), " size_asym avail:", mean(!is.na(d$size_asym)), " polity avail:", mean(!is.na(d$min_polity)), "\n")
zs <- function(x) as.numeric(scale(x))
BASE <- "z_mut + z_cap + l_dur + l_unbal + l_bal"
sink("results/hazard_added_controls_output.txt", split = TRUE)
specs <- list("baseline" = "", "+ contiguity" = "+ contig", "+ economic-size asymmetry" = "+ z_size",
              "+ regime type (min Polity2, gap)" = "+ z_minpol + z_polgap", "+ joint democracy" = "+ joint_dem",
              "+ all" = "+ contig + z_size + z_minpol + z_polgap")
for (sg in c(-1, 1)) for (y in c("flip_i", "gone_i")) {
  cat(sprintf("\n-- %s / %s\n", ifelse(sg == -1, "DISPUTES", "ALLIANCES"), y))
  # common sample: all controls available, so that specs are comparable
  x <- d %>% filter(sign == sg, !is.na(dep), !is.na(mut), !is.na(cap), !is.na(size_asym), !is.na(min_polity)) %>%
    mutate(z_dep = zs(dep), z_mut = zs(mut), z_cap = zs(cap), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal),
           z_size = zs(size_asym), z_minpol = zs(min_polity), z_polgap = zs(polity_gap))
  for (nm in names(specs)) {
    m <- glm(as.formula(paste(y, "~ z_dep +", BASE, specs[[nm]])), data = x, family = binomial)
    ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
    extra <- if (nm == "+ all") sprintf(" | contig=%.3f (p=%.3f) size=%.3f (p=%.3f) minpol=%.3f (p=%.3f)", ct["contig",1], ct["contig",4], ct["z_size",1], ct["z_size",4], ct["z_minpol",1], ct["z_minpol",4]) else ""
    cat(sprintf("  %-34s dep=%6.3f (se %.3f, p=%.4f) | N=%d ev=%d%s\n", nm, ct["z_dep",1], ct["z_dep",2], ct["z_dep",4], nobs(m), sum(x[[y]]), extra))
  }
}
sink()
