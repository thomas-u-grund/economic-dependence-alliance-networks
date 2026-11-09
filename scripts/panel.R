# panel.R: builds the dyad-year analysis panel `d` from raw/ and data/tie_trade.csv.gz. Loaded once by setup.R; writes no output.
#     Panel: alliance (a) and dispute (m) status per dyad-year, larger/smaller trade share, capability asymmetry, contiguity,
#     durations, and the numbers of third states whose two ties to the pair have the same sign (n_same) or opposite signs (n_diff).

al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
  filter((defense == 1 | (Sys.getenv("ALLIANCE_DEF") != "defense" & (neutrality == 1 | nonaggression == 1 | entente == 1)))) %>%
  transmute(year, lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2), a = 1L) %>% distinct()
mi <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>%
  mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>% filter(hihost >= 3) %>% transmute(year, lo, hi, m = 1L) %>% distinct()
contig <- read_csv("raw/contiguity/DirectContiguity320/contdird.csv", show_col_types = FALSE) %>%
  transmute(year, lo = pmin(state1no, state2no), hi = pmax(state1no, state2no), contig = 1L) %>% distinct()
tt <- read_csv("data/tie_trade.csv.gz", show_col_types = FALSE) %>%
  transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), shl = trade_share_low, shh = trade_share_high, cap = cap_asymmetry, vol = total_dyadic_trade) %>%
  filter(year >= 1870, year <= 2011, !is.na(shl), !is.na(shh), !is.na(cap), shl <= 1, shh <= 1)
cat("trade dyad-years with valid shares:", nrow(tt), "\n")
d <- tt %>% left_join(al, by = c("year", "lo", "hi")) %>% left_join(mi, by = c("year", "lo", "hi")) %>% left_join(contig, by = c("year", "lo", "hi")) %>%
  mutate(a = replace_na(a, 0L), m = replace_na(m, 0L), contig = replace_na(contig, 0L), mx = pmax(shl, shh), mn = pmin(shl, shh), dyad = paste(lo, hi), lvol = log1p(vol), decade = factor(floor(year / 10)))
# next-year layers (only when the dyad is observed with valid trade data next year)
nx <- d %>% transmute(year = year - 1, lo, hi, a1 = a, m1 = m)
d <- d %>% left_join(nx, by = c("year", "lo", "hi")) %>% filter(!is.na(a1))
# own-layer durations: consecutive observed years with the same status
d <- d %>% arrange(dyad, year) %>% group_by(dyad) %>%
  mutate(ra = cumsum(c(TRUE, diff(year) != 1 | diff(a) != 0)), rm_ = cumsum(c(TRUE, diff(year) != 1 | diff(m) != 0))) %>%
  group_by(dyad, ra) %>% mutate(dur_a = row_number()) %>% group_by(dyad, rm_) %>% mutate(dur_m = row_number()) %>% ungroup()
# third-state structure from the full signed network of each year (dispute overrides alliance)
raw_edges <- bind_rows(al %>% transmute(year, lo, hi, s = 1L), mi %>% transmute(year, lo, hi, s = -1L)) %>%
  group_by(year, lo, hi) %>% summarise(s = ifelse(any(s == -1L), -1L, 1L), .groups = "drop")
cnt <- lapply(split(raw_edges, raw_edges$year), function(e) {
  ids <- sort(unique(c(e$lo, e$hi))); n <- length(ids); ix <- setNames(seq_len(n), ids)
  P <- matrix(0, n, n); N <- matrix(0, n, n)
  p <- e$s == 1L
  P[cbind(ix[as.character(e$lo[p])], ix[as.character(e$hi[p])])] <- 1; P <- pmax(P, t(P))
  N[cbind(ix[as.character(e$lo[!p])], ix[as.character(e$hi[!p])])] <- 1; N <- pmax(N, t(N))
  same <- P %*% P + N %*% N; diff <- P %*% N + N %*% P
  list(ids = ids, ix = ix, same = same, diff = diff)
})
get_cnt <- function(yr, lo, hi) {
  out_s <- rep(0, length(yr)); out_d <- rep(0, length(yr))
  for (y in unique(yr)) { k <- which(yr == y); c0 <- cnt[[as.character(y)]]; if (is.null(c0)) next
    i <- c0$ix[as.character(lo[k])]; j <- c0$ix[as.character(hi[k])]; ok <- !is.na(i) & !is.na(j)
    out_s[k[ok]] <- c0$same[cbind(i[ok], j[ok])]; out_d[k[ok]] <- c0$diff[cbind(i[ok], j[ok])] }
  list(s = out_s, d = out_d)
}
cc <- get_cnt(d$year, d$lo, d$hi); d$n_same <- cc$s; d$n_diff <- cc$d
zs <- function(x) as.numeric(scale(x))
# standardise once on the whole panel so that coefficients are comparable across cells
d <- d %>% mutate(z_mx = zs(mx), z_mn = zs(mn), z_cap = zs(cap), z_vol = zs(lvol), l_same = log1p(n_same), l_diff = log1p(n_diff))
run <- function(x, y, label, extra = "") {
  durv <- if (grepl("dispute layer", label)) "dur_m" else "dur_a"
  x$l_dur <- log(x[[durv]])
  m <- glm(as.formula(paste(y, "~ z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff", extra)), data = x, family = binomial)
  ct <- coeftest(m, vcov = vcovCL(m, cluster = x$dyad))
  cat(sprintf("%-62s larger %6.3f (se %.3f, p=%.4f) | smaller %6.3f | events %d of %d | dyads with event %d\n", label, ct["z_mx",1], ct["z_mx",2], ct["z_mx",4], ct["z_mn",1], sum(x[[y]]), nrow(x), n_distinct(x$dyad[x[[y]] == 1])))
  invisible(list(m = m, ct = ct))
}
