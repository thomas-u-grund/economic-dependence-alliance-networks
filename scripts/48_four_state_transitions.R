# 48: Four relational states (A alliance only, B alliance + dispute, C dispute only) and transitions out of each,
#     with the larger trade share as predictor. D (neither) is the absorbing 'no signed tie' state; onset from D is not modelled.
suppressMessages({ library(sandwich); library(lmtest) })
source("scripts/38_hazard_panel_lib.R")
d0 <- build_panel(base_signed(), "baseline")
al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>%
  transmute(year, lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), a = 1L) %>% distinct()
mi <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>%
  mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>% filter(hihost >= 3) %>%
  transmute(year, lo, hi, m = 1L) %>% distinct()
st <- function(a, m) ifelse(a == 1 & m == 1, "B", ifelse(a == 1, "A", ifelse(m == 1, "C", "D")))
d <- d0 %>% left_join(al, by = c("year", "lo", "hi")) %>% left_join(mi, by = c("year", "lo", "hi")) %>%
  mutate(a = replace_na(a, 0L), m = replace_na(m, 0L), s0 = st(a, m)) %>%
  left_join(al %>% transmute(year = year - 1, lo, hi, a1 = a), by = c("year", "lo", "hi")) %>%
  left_join(mi %>% transmute(year = year - 1, lo, hi, m1 = m), by = c("year", "lo", "hi")) %>%
  mutate(a1 = replace_na(a1, 0L), m1 = replace_na(m1, 0L), s1 = st(a1, m1)) %>%
  filter(!is.na(dep), !is.na(mut), !is.na(cap))
sink("results/four_state_transitions_output.txt", split = TRUE)
cat("Transition counts (rows = state at t, columns = state at t+1), sample with trade data\n")
print(table(from = d$s0, to = d$s1))
cat("\nRow shares\n"); print(round(prop.table(table(from = d$s0, to = d$s1), 1), 4))
zs <- function(x) as.numeric(scale(x))
d <- d %>% mutate(z_mx = zs(mx), z_mut = zs(mut), z_cap = zs(cap), l_unbal = log1p(n_unbal), l_bal = log1p(n_bal))  # one scale for all rows
cat("\nLarger-share coefficient, transition vs all other outcomes from the same starting state (logit, dyad-clustered SE)\n")
for (from in c("A", "B", "C")) for (to in c("A", "B", "C", "D")) {
  if (from == to) next
  x <- d %>% filter(s0 == from) %>% mutate(ev = as.integer(s1 == to))
  if (sum(x$ev) < 25) { cat(sprintf("%s -> %s: %d events (too few)\n", from, to, sum(x$ev))); next }
  x <- x %>% mutate(l_dur = log(dur))
  m <- glm(ev ~ z_mx + z_mut + z_cap + l_dur + l_unbal + l_bal, x, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%s -> %s: larger share %6.3f (se %.3f, p=%.4f) | smaller %6.3f | events %d of %d, dyads with event %d\n", from, to, ct["z_mx",1], ct["z_mx",2], ct["z_mx",4], ct["z_mut",1], sum(x$ev), nrow(x), n_distinct(x$dyad[x$ev == 1])))
}
sink()
