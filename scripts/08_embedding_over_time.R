# 08: Does structural embedding (T) and balance (common/split partners) matter differently across eras? Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest); library(car) })

src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]
eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1), end_m = as.integer(m1 == 0), tri = (n_same + n_diff) >= 1, T = as.integer(tri),
  era = cut(year, c(1869, 1913, 1945, 1989, 2011), labels = c("1870-1913", "1914-1945", "1946-1989", "1990-2011")), era = relevel(era, "1946-1989"), yc = (year - 1950) / 50)
sink("results/embedding_over_time_output.txt", split = TRUE)
cat("== (a) Descriptives by era ==\n")
a <- d %>% group_by(era) %>% summarise(dyad_years = n(), pct_in_T = round(100 * mean(T), 1), onsets = sum(m == 0 & onset_m == 1),
  pct_onsets_in_T = round(100 * mean(T[m == 0 & onset_m == 1]), 1), rate_T = round(100 * mean(onset_m[m == 0 & T == 1]), 2), rate_notT = round(100 * mean(onset_m[m == 0 & T == 0]), 2), .groups = "drop") %>% mutate(rate_ratio = round(rate_T / rate_notT, 1))
print(as.data.frame(a))
cat("\n== by decade (T share, onsets, rate ratio) ==\n")
print(as.data.frame(d %>% group_by(decade) %>% summarise(pct_T = round(100 * mean(T), 1), onsets = sum(m == 0 & onset_m == 1), rate_T = round(100 * mean(onset_m[m == 0 & T == 1]), 2), rate_notT = round(100 * mean(onset_m[m == 0 & T == 0]), 2), .groups = "drop")))
on <- d %>% filter(m == 0) %>% mutate(l_dur = log(dur_m))
wald <- function(b, V, idx) { i <- match(idx, names(b)); w <- drop(t(b[i]) %*% solve(V[i, i]) %*% b[i]); cat(sprintf("  joint Wald chi2(%d) = %.2f, p = %.4f\n", length(i), w, pchisq(w, length(i), lower.tail = FALSE))) }
cl <- function(m, x) coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
cat("\n== (b) onset ~ T*era + controls + decade FE ==\n")
m1 <- glm(onset_m ~ T * era + z_mx + z_mn + z_cap + contig + l_dur + decade, on, family = binomial); b <- coef(m1); b <- b[!is.na(b)]; V <- vcovCL(m1, cluster = on$dyad)[names(b), names(b)]
for (e in levels(on$era)) { nm <- if (e == "1946-1989") "T" else paste0("T:era", e); idx <- c("T", if (e != "1946-1989") nm)
  est <- sum(b[idx]); L <- setNames(rep(0, length(b)), names(b)); L[idx] <- 1; se <- sqrt(drop(t(L) %*% V %*% L)); cat(sprintf("  %-10s T log-odds %6.3f (se %.3f)  OR %.2f\n", e, est, se, exp(est))) }
wald(b, V, grep("^T:era", names(b), value = TRUE))
cat("\n== (c) balance by era within T: onset ==\n")
t1 <- on %>% filter(T == 1)
m2 <- glm(onset_m ~ (l_diff + l_same) * era + z_mx + z_mn + z_cap + contig + l_dur + decade, t1, family = binomial); b2 <- coef(m2); b2 <- b2[!is.na(b2)]; V2 <- vcovCL(m2, cluster = t1$dyad)[names(b2), names(b2)]
for (v in c("l_diff", "l_same")) for (e in levels(t1$era)) { idx <- c(v, if (e != "1946-1989") paste0(v, ":era", e)); L <- setNames(rep(0, length(b2)), names(b2)); L[idx] <- 1
  cat(sprintf("  %-7s %-10s %6.3f (se %.3f)\n", v, e, sum(b2[idx]), sqrt(drop(t(L) %*% V2 %*% L)))) }
wald(b2, V2, grep("^l_(diff|same):era", names(b2), value = TRUE))
m3 <- glm(onset_m ~ (l_diff + l_same) * yc + z_mx + z_mn + z_cap + contig + l_dur, t1, family = binomial); print(cl(m3, t1)[c("l_diff", "l_same", "yc", "l_diff:yc", "l_same:yc"), ])
cat("\n== (c2) balance by era within T: ending ==\n")
en <- d %>% filter(m == 1, T == 1) %>% mutate(l_dur = log(dur_m))
m4 <- glm(end_m ~ (l_diff + l_same) * era + z_mx + z_mn + z_cap + contig + l_dur, en, family = binomial); b4 <- coef(m4); b4 <- b4[!is.na(b4)]; V4 <- vcovCL(m4, cluster = en$dyad)[names(b4), names(b4)]
for (v in c("l_diff", "l_same")) for (e in levels(en$era)) { idx <- c(v, if (e != "1946-1989") paste0(v, ":era", e)); L <- setNames(rep(0, length(b4)), names(b4)); L[idx] <- 1
  cat(sprintf("  %-7s %-10s %6.3f (se %.3f)\n", v, e, sum(b4[idx]), sqrt(drop(t(L) %*% V4 %*% L)))) }
wald(b4, V4, grep("^l_(diff|same):era", names(b4), value = TRUE)); cat("ending events by era:\n"); print(table(en$era, en$end_m))
cat("\n== (d) exposure (larger share) by era within T: onset ==\n")
m5 <- glm(onset_m ~ z_mx * era + z_mn + z_cap + contig + l_dur + l_same + l_diff + decade, t1, family = binomial); b5 <- coef(m5); b5 <- b5[!is.na(b5)]; V5 <- vcovCL(m5, cluster = t1$dyad)[names(b5), names(b5)]
for (e in levels(t1$era)) { idx <- c("z_mx", if (e != "1946-1989") paste0("z_mx:era", e)); L <- setNames(rep(0, length(b5)), names(b5)); L[idx] <- 1; cat(sprintf("  %-10s %6.3f (se %.3f)\n", e, sum(b5[idx]), sqrt(drop(t(L) %*% V5 %*% L)))) }
wald(b5, V5, grep("^z_mx:era", names(b5), value = TRUE))
cat("\nonsets in T by era:\n"); print(table(t1$era, t1$onset_m))
sink()
