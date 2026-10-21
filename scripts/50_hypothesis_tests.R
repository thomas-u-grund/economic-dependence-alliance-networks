# 50: Pre-specified tests: T1 weaker vs stronger exposed state, T2 alliance depth, T3 timing placebo.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
tt <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), shl = trade_share_low, shh = trade_share_high, cl = cinc_low, ch = cinc_high)
# lag 1 and lead 2 of the larger share
sh <- tt %>% transmute(year, lo, hi, mx_ = pmax(shl, shh))
d0 <- d0 %>% left_join(sh %>% transmute(year = year + 1, lo, hi, mx_lag1 = mx_), by = c("year", "lo", "hi")) %>%
  left_join(sh %>% transmute(year = year - 2, lo, hi, mx_lead2 = mx_), by = c("year", "lo", "hi")) %>%
  left_join(tt %>% select(year, lo, hi, shl, shh, cl, ch) %>% rename(sl = shl, sh_ = shh), by = c("year", "lo", "hi"))
d0 <- d0 %>% mutate(weak_exp = case_when(is.na(cl) | is.na(ch) | sl == sh_ ~ NA_real_,
                                          (sl > sh_ & cl < ch) | (sh_ > sl & ch < cl) ~ 1, TRUE ~ 0),
                    end_i = as.integer(flip | gone))
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>%
  mutate(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2)) %>%
  group_by(year, lo, hi) %>% summarise(def = as.integer(max(defense == 1)), .groups = "drop")
mem <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_member.csv", show_col_types = FALSE)
big <- mem %>% group_by(version4id) %>% summarise(n = n_distinct(ccode)) %>% filter(n >= 8) %>% pull(version4id)
mb <- mem %>% filter(version4id %in% big) %>% transmute(id = version4id, ccode, y0 = mem_st_year, y1 = ifelse(is.na(mem_end_year) | mem_end_year == 0, 9999, mem_end_year))
inbig <- function(x) { r <- rep(FALSE, nrow(x)); for (id in big) { m <- mb[mb$id == id, ]
  f <- function(cc, yr) sapply(seq_along(cc), function(i) any(m$ccode == cc[i] & m$y0 <= yr[i] & m$y1 >= yr[i])); r <- r | (f(x$lo, x$year) & f(x$hi, x$year)) }; as.integer(r) }
prep <- function(x) x %>% filter(!is.na(dep), !is.na(mut), !is.na(cap)) %>%
  mutate(z_dep = as.numeric(scale(dep)), z_mut = as.numeric(scale(mut)), z_cap = as.numeric(scale(cap)), l_dur = log(dur), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))
ctrl <- "z_mut + z_cap + l_dur + l_unbal + l_bal"
rep_ <- function(m, x, terms, label) { ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%-46s", label)); for (t in terms) cat(sprintf(" | %s %6.3f (se %.3f, p=%.4f)", t, ct[t,1], ct[t,2], ct[t,4])); cat(sprintf(" | N=%d events=%d\n", nobs(m), sum(m$y))); }
sink("results/hypothesis_tests_output.txt", split = TRUE)
cat("== T1: larger share x (more exposed state is the weaker state) ==\n")
for (sg in c(-1, 1)) { x <- prep(d0 %>% filter(sign == sg, !is.na(weak_exp)))
  cat(sprintf("%s: share of dyad-years where the more exposed state is weaker: %.2f\n", ifelse(sg == -1, "Disputes", "Alliances"), mean(x$weak_exp)))
  m <- glm(as.formula(paste("flip_i ~ z_dep * weak_exp +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "weak_exp", "z_dep:weak_exp"), ifelse(sg == -1, "disputes sign change", "alliances sign change")) }
x <- prep(d0 %>% filter(sign == -1, !is.na(weak_exp)) %>% inner_join(al, by = c("year", "lo", "hi")))
m <- glm(as.formula(paste("end_i ~ z_dep * weak_exp +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "z_dep:weak_exp"), "disputes ending inside alliance")
cat("\n== T2: alliance depth ==\n")
x <- prep(d0 %>% filter(sign == -1) %>% inner_join(al, by = c("year", "lo", "hi"))); x$big <- inbig(x)
if (length(unique(x$def)) > 1) { m <- glm(as.formula(paste("end_i ~ z_dep * def +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "def", "z_dep:def"), "disputes ending x defense pact") } else cat("disputes ending x defense pact: skipped (all allied dyads have a defence pact under this definition)\n")
m <- glm(as.formula(paste("end_i ~ z_dep * big +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "big", "z_dep:big"), "disputes ending x both in 8+ alliance")
x <- prep(d0 %>% filter(sign == 1) %>% inner_join(al, by = c("year", "lo", "hi"))); x$big <- inbig(x)
if (length(unique(x$def)) > 1) { m <- glm(as.formula(paste("flip_i ~ z_dep * def +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "def", "z_dep:def"), "alliances turn hostile x defense pact") } else cat("alliances turn hostile x defense pact: skipped (all allied dyads have a defence pact under this definition)\n")
m <- glm(as.formula(paste("flip_i ~ z_dep * big +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "big", "z_dep:big"), "alliances turn hostile x 8+ alliance")
cat("\n== T3: timing placebo (current, lag 1, lead 2 of the larger share; common sample) ==\n")
for (sg in c(-1, 1)) { x <- d0 %>% filter(sign == sg, !is.na(mx_lag1), !is.na(mx_lead2)); x <- prep(x)
  m0 <- 1; mu <- mean(x$dep); sdv <- sd(x$dep); x$z_lag <- (x$mx_lag1 - mu) / sdv; x$z_lead <- (x$mx_lead2 - mu) / sdv
  cat(sprintf("%s: corr(current, lead2) = %.2f\n", ifelse(sg == -1, "Disputes", "Alliances"), cor(x$z_dep, x$z_lead)))
  m <- glm(as.formula(paste("flip_i ~ z_dep + z_lag + z_lead +", ctrl)), x, family = binomial); rep_(m, x, c("z_dep", "z_lag", "z_lead"), ifelse(sg == -1, "disputes sign change", "alliances sign change"))
  m <- glm(as.formula(paste("flip_i ~ z_dep +", ctrl)), x, family = binomial); rep_(m, x, "z_dep", "  (current only, same sample)") }
sink()
