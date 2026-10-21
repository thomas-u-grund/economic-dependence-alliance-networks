# 38: Library. Builds the dyad-year panel for the hazard analyses from ANY signed-ties table
#     (year, lo, hi, sign), so alternative alliance/dispute definitions can be swapped in.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(igraph) })


pressure_counts <- function(st) {   # n_unbal, n_bal per (year, lo, hi): triads a sign flip would fix / break
  out <- lapply(sort(unique(st$year)), function(yr) {
    e <- st[st$year == yr, ]
    if (nrow(e) < 3) return(NULL)
    g <- graph_from_data_frame(data.frame(a = as.character(e$lo), b = as.character(e$hi)), directed = FALSE)
    tri <- igraph::triangles(g); if (length(tri) == 0) return(NULL)
    tm <- matrix(V(g)$name[tri], ncol = 3, byrow = TRUE)
    key <- function(x, y) paste(pmin(as.numeric(x), as.numeric(y)), pmax(as.numeric(x), as.numeric(y)))
    sgn <- setNames(e$sign, key(e$lo, e$hi))
    k12 <- key(tm[,1], tm[,2]); k13 <- key(tm[,1], tm[,3]); k23 <- key(tm[,2], tm[,3])
    bal <- (sgn[k12] * sgn[k13] * sgn[k23]) > 0
    tibble(key = c(k12, k13, k23), bal = rep(bal, 3)) %>% group_by(key) %>%
      summarise(n_unbal = sum(!bal), n_bal = sum(bal), .groups = "drop") %>% mutate(year = yr)
  }) %>% bind_rows()
  out %>% separate(key, c("lo", "hi"), sep = " ", convert = TRUE)
}

build_panel <- function(st, label = "") {
  st <- st %>% filter(year >= 1816, year <= 2012) %>% distinct(year, lo, hi, .keep_all = TRUE)
  ed <- pressure_counts(st)
  nx <- st %>% transmute(year = year - 1, lo, hi, sign_next = sign)
  tt <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
    transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), dep = pmax(trade_share_low, trade_share_high), asym = dep_asymmetry,
              shl = trade_share_low, shh = trade_share_high, cap = cap_asymmetry)
  d <- st %>% left_join(ed, by = c("year", "lo", "hi")) %>% left_join(nx, by = c("year", "lo", "hi")) %>%
    left_join(tt, by = c("year", "lo", "hi")) %>% filter(year < 2012) %>%
    mutate(n_unbal = replace_na(n_unbal, 0L), n_bal = replace_na(n_bal, 0L), dyad = paste(lo, hi),
           mut = pmin(shl, shh), mx = pmax(shl, shh),
           flip = !is.na(sign_next) & sign_next != sign, gone = is.na(sign_next),
           any_tri = (n_unbal + n_bal) > 0) %>%
    arrange(dyad, year) %>% group_by(dyad) %>% mutate(run = cumsum(c(TRUE, diff(year) != 1 | diff(sign) != 0))) %>%
    group_by(dyad, run) %>% mutate(dur = row_number()) %>% ungroup() %>%
    mutate(flip_i = as.integer(flip), gone_i = as.integer(gone), prior_spells = run - 1)
  if (label != "") cat(sprintf("[%s] dyad-years %d, disputes %d, alliances %d, flips %d, gone %d\n", label, nrow(d),
      sum(d$sign == -1), sum(d$sign == 1), sum(d$flip), sum(d$gone)))
  d
}

base_signed <- function() {
  # ALLIANCE_DEF=defense rebuilds the signed network with COW defence pacts only (Type I); otherwise all four COW commitment types.
  if (Sys.getenv("ALLIANCE_DEF") == "defense") return(alt_signed(defense_only = TRUE))
  read_csv("data/signed_ties.csv", show_col_types = FALSE) %>%
    filter(year >= 1816, year <= 2012) %>% transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), sign, multiplex)
}

alt_signed <- function(defense_only = FALSE, hihost_min = 3) {
  al <- read_csv("raw/alliances/version4.1_csv/alliance_v4.1_by_dyad_yearly.csv", show_col_types = FALSE) %>%
    mutate(lo = pmin(ccode1, ccode2), hi = pmax(ccode1, ccode2))
  al <- if (defense_only) al %>% filter(defense == 1) else al %>% filter(defense == 1 | neutrality == 1 | nonaggression == 1 | entente == 1)
  al <- al %>% transmute(year, lo, hi, sign = 1) %>% distinct()
  mi <- read_csv("raw/mids/dyadic_mid_4.03_update/dyadic_mid_4.03.csv", show_col_types = FALSE) %>%
    mutate(lo = pmin(statea, stateb), hi = pmax(statea, stateb)) %>% filter(hihost >= hihost_min) %>%
    transmute(year, lo, hi, sign = -1) %>% distinct()
  bind_rows(al, mi) %>% group_by(year, lo, hi) %>%
    summarise(multiplex = n_distinct(sign) > 1, sign = ifelse(any(sign == -1), -1, 1), .groups = "drop") %>%
    filter(year >= 1816, year <= 2012)
}
