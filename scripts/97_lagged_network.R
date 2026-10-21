# 97: FULL PANEL: partner counts measured in the network of year t-1 (lagged) vs year t, onset and alliance formation. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
cl <- get_cnt(d$year - 1, d$lo, d$hi); d$l_same_l <- log1p(cl$s); d$l_diff_l <- log1p(cl$d)
d <- d %>% mutate(onset_m = as.integer(m1 == 1), form_a = as.integer(a1 == 1))
ci <- function(est, se) sprintf("%6.3f [%6.3f, %6.3f]", est, est - 1.96 * se, est + 1.96 * se)
sink("results/lagged_network_output.txt", split = TRUE)
cat("FULL PANEL: contemporaneous (t) vs lagged (t-1) partner counts; dyad-clustered SEs\n")
for (p in list(list("dispute onset", "onset_m", "dur_m", d %>% filter(m == 0)), list("alliance formation", "form_a", "dur_a", d %>% filter(a == 0, m == 0)))) for (sp in c("", "+ z_vol + decade")) for (tm in c("t", "t-1")) {
  x <- p[[4]]; x$l_dur <- log(x[[p[[3]]]]); s <- if (tm == "t") "l_same" else "l_same_l"; df <- if (tm == "t") "l_diff" else "l_diff_l"
  m <- glm(as.formula(paste(p[[2]], "~ z_mx + z_mn + z_cap + contig + l_dur +", s, "+", df, sp)), x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m)
  g <- function(v) ci(b[v], sqrt(V[v, v]))
  cat(sprintf("%-19s %-18s partners at %-3s | larger share %s | split %s | common %s | N=%d events=%d\n", p[[1]], ifelse(sp == "", "baseline", "+ vol + decade"), tm, g("z_mx"), g(df), g(s), nrow(x), sum(x[[p[[2]]]]))) }
sink()
