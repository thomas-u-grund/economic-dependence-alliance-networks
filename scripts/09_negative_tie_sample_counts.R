# 09: Feasibility counts for the negative-tie resolution design (descriptive only; no outcome models). Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr) })

src <- readLines("scripts/03_layered_hazards.R"); k <- grep("^sink\\(", src)[1]
eval(parse(text = src[1:(k - 1)]))
# split common partners: PP (common allies, +,+ with the dispute = unbalanced), NN (common enemies, -,-), PN (split = balanced)
pn <- lapply(split(raw_edges, raw_edges$year), function(e) {
  ids <- sort(unique(c(e$lo, e$hi))); n <- length(ids); ix <- setNames(seq_len(n), ids); p <- e$s == 1L
  P <- matrix(0, n, n); N <- matrix(0, n, n)
  P[cbind(ix[as.character(e$lo[p])], ix[as.character(e$hi[p])])] <- 1; P <- pmax(P, t(P))
  N[cbind(ix[as.character(e$lo[!p])], ix[as.character(e$hi[!p])])] <- 1; N <- pmax(N, t(N))
  list(ix = ix, pp = P %*% P, nn = N %*% N) })
g <- function(yr, lo, hi, w) { o <- rep(0, length(yr)); for (y in unique(yr)) { k <- which(yr == y); c0 <- pn[[as.character(y)]]; if (is.null(c0)) next
  i <- c0$ix[as.character(lo[k])]; j <- c0$ix[as.character(hi[k])]; ok <- !is.na(i) & !is.na(j); o[k[ok]] <- c0[[w]][cbind(i[ok], j[ok])] }; o }
d$n_pp <- g(d$year, d$lo, d$hi, "pp"); d$n_nn <- g(d$year, d$lo, d$hi, "nn")
dm <- d %>% filter(m == 1) %>% mutate(end_m = as.integer(m1 == 0))
cat(sprintf("Dispute dyad-years (next year observed): %d; spells (dyad x run): %d; ending events: %d\n", nrow(dm), n_distinct(dm$dyad, dm$rm_), sum(dm$end_m)))
cat("\nBy common allies (n_pp) category:\n"); print(as.data.frame(dm %>% mutate(pp = pmin(n_pp, 3)) %>% group_by(pp) %>% summarise(dyad_years = n(), spells = n_distinct(dyad, rm_), ends = sum(end_m), end_rate_pct = round(100 * mean(end_m), 1), .groups = "drop")))
cat("\nBy common enemies (n_nn) category:\n"); print(as.data.frame(dm %>% mutate(nn = pmin(n_nn, 3)) %>% group_by(nn) %>% summarise(dyad_years = n(), spells = n_distinct(dyad, rm_), ends = sum(end_m), end_rate_pct = round(100 * mean(end_m), 1), .groups = "drop")))
cat("\nSample for the design (n_pp >= 1): dyad-years, spells, ends, and how many are allied dyads (a==1):\n")
s <- dm %>% filter(n_pp >= 1); cat(sprintf("  %d dyad-years, %d spells, %d ends, %d dyads; allied dyad-years %d (%.0f%%); contiguous %d (%.0f%%)\n", nrow(s), n_distinct(s$dyad, s$rm_), sum(s$end_m), n_distinct(s$dyad), sum(s$a == 1), 100 * mean(s$a == 1), sum(s$contig == 1), 100 * mean(s$contig == 1)))
cat("  n_diff (balanced triads an ending would break) among them: "); print(summary(s$n_diff))
cat("  correlation n_pp vs n_diff:", round(cor(s$n_pp, s$n_diff), 2), "\n")
cat("  by era (year):\n"); print(table(cut(s$year, c(1869, 1913, 1945, 1989, 2011)), s$end_m))
cat("\nSpells whose n_pp drops to 0 while the dispute continues (informative exit from the risk set):\n")
dm2 <- dm %>% arrange(dyad, rm_, year) %>% group_by(dyad, rm_) %>% mutate(had = cummax(n_pp >= 1)) %>% ungroup()
cat(sprintf("  dyad-years with n_pp==0 after an earlier n_pp>=1 in the same spell: %d\n", sum(dm2$had == 1 & dm2$n_pp == 0)))
