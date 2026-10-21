# 46: Larger share (maximum exposure) as primary predictor; compare with conventional trade-intensity measures.
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
  mutate(vol = shl * tt_lo, log_vol = log(vol), comb = vol * 2 / (tt_lo + tt_hi), log_comb = log(comb),
         size_asym = abs(log(tt_lo) - log(tt_hi)), min_polity = pmin(p_lo, p_hi), polity_gap = abs(p_lo - p_hi),
         decade = factor(floor(year / 10)))
zs <- function(x) as.numeric(scale(x))
CTRL <- "z_cap + l_dur + l_unbal + l_bal"
sink("results/max_exposure_rivals_output.txt", split = TRUE)
for (sg in c(-1, 1)) for (y in c("flip_i", "gone_i")) {
  cat(sprintf("\n==== %s / %s\n", ifelse(sg == -1, "DISPUTES", "ALLIANCES"), y))
  x <- d %>% filter(sign == sg, !is.na(mx), !is.na(mut), !is.na(cap), !is.na(log_vol), !is.na(size_asym), !is.na(min_polity), vol > 0) %>%
    mutate(z_mx = zs(mx), z_mut = zs(mut), z_dep = zs(mx - mut), z_cap = zs(cap), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal),
           z_vol = zs(log_vol), z_comb = zs(log_comb), z_size = zs(size_asym), z_minpol = zs(min_polity), z_polgap = zs(polity_gap),
           z_mxmut = z_mx * z_mut)
  cat("N =", nrow(x), " events =", sum(x[[y]]), "\n")
  specs <- list(
    "A asymmetry + min (published)"          = "z_dep + z_mut",
    "B larger + smaller share"               = "z_mx + z_mut",
    "C larger + smaller + larger x smaller"  = "z_mx + z_mut + z_mxmut",
    "D B + log dyadic trade volume"          = "z_mx + z_mut + z_vol",
    "E B + log trade intensity (comb. share)"= "z_mx + z_mut + z_comb",
    "F B + volume + all controls"            = "z_mx + z_mut + z_vol + contig + z_size + z_minpol + z_polgap",
    "G B + decade FE"                        = "z_mx + z_mut + decade",
    "H volume only (no shares)"              = "z_vol",
    "I intensity only (no shares)"           = "z_comb")
  for (nm in names(specs)) {
    m <- glm(as.formula(paste(y, "~", specs[[nm]], "+", CTRL)), data = x, family = binomial)
    ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
    keep <- intersect(c("z_dep", "z_mx", "z_mut", "z_mxmut", "z_vol", "z_comb"), rownames(ct))
    cat(sprintf("%-42s AIC=%.1f | %s\n", nm, AIC(m), paste(sprintf("%s=%.3f (p=%.3f)", sub("^z_", "", keep), ct[keep, 1], ct[keep, 4]), collapse = "  ")))
  }
}
sink()
