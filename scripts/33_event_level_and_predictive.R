# ------------------------------------------------------------------
# 33: Event-level weighting, restricted-overlap, transition-type split,
#     out-of-sample predictive gain, and diversification heterogeneity.
#
#     (a) A single political episode (one dyad changing sign in a given
#         year) can appear in many triad-events. Re-estimate M4 giving
#         each unique (changed dyad, year, direction) total weight 1.
#     (b) Restrict to episodes appearing in few triads.
#     (c) Separate dispute->alliance (alliance-preserving dispute
#         endings) from alliance->dispute transitions.
#     (d) Out-of-sample gain from adding dependence asymmetry, with
#         folds that keep episodes / blocs together.
#     (e) Heterogeneity by trade diversification of the more exposed
#         state (partner-share HHI).
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(survival); library(sandwich); library(lmtest) })

set.seed(20261020)

h3 <- readRDS("results/h3_displacement.rds")
strict <- h3$strict %>% mutate(event_id = row_number())
tie_trade <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high))
tt <- tie_trade %>% select(year, ccode_low, ccode_high, trade_share_low, trade_share_high)
get_shares <- function(df, a, b_, suffix) {
  df %>% mutate(lo = as.character(pmin(as.numeric(.data[[a]]), as.numeric(.data[[b_]]))),
                hi = as.character(pmax(as.numeric(.data[[a]]), as.numeric(.data[[b_]])))) %>%
    left_join(tt, by = c("t1"="year","lo"="ccode_low","hi"="ccode_high")) %>%
    rename_with(~ paste0(., "_", suffix), c(trade_share_low, trade_share_high)) %>%
    rename(!!paste0("lo_", suffix) := lo, !!paste0("hi_", suffix) := hi)
}
strict <- strict %>% get_shares("node1","node2","12") %>% get_shares("node1","node3","13") %>% get_shares("node2","node3","23") %>%
  mutate(across(starts_with("trade_share_"), as.numeric))
for (s in c("12","13","23"))
  strict[[paste0("mutual_min_", s)]] <- pmin(strict[[paste0("trade_share_low_", s)]], strict[[paste0("trade_share_high_", s)]])

# diversification: HHI of each state's partner trade shares, by year
hhi <- bind_rows(
  tie_trade %>% transmute(year, ccode = ccode_low,  sh = trade_share_low),
  tie_trade %>% transmute(year, ccode = ccode_high, sh = trade_share_high)) %>%
  filter(!is.na(sh)) %>% group_by(year, ccode) %>% summarise(hhi = sum(sh^2), .groups = "drop")
for (s in c("12","13","23")) {
  lo <- strict[[paste0("lo_", s)]]; hi <- strict[[paste0("hi_", s)]]
  shl <- strict[[paste0("trade_share_low_", s)]]; shh <- strict[[paste0("trade_share_high_", s)]]
  exposed <- ifelse(shl >= shh, lo, hi)
  d <- tibble(year = strict$t1, ccode = exposed) %>% left_join(hhi, by = c("year","ccode"))
  strict[[paste0("hhi_exp_", s)]] <- d$hhi
  # concentration of the exposed state's trade among its OTHER partners (excludes the focal tie)
  sh_exp <- pmax(shl, shh)
  strict[[paste0("hhio_exp_", s)]] <- pmax(d$hhi - sh_exp^2, 0)
}

# unique political episode: the changed dyad x year x direction
dy <- function(a, b) paste(pmin(as.numeric(a), as.numeric(b)), pmax(as.numeric(a), as.numeric(b)), sep = "-")
strict <- strict %>% mutate(
  d12 = dy(node1,node2), d13 = dy(node1,node3), d23 = dy(node2,node3),
  ch_dyad = case_when(changed12 ~ d12, changed13 ~ d13, changed23 ~ d23),
  ch_sign = case_when(changed12 ~ sign12_t1, changed13 ~ sign13_t1, changed23 ~ sign23_t1),
  episode = paste(ch_dyad, t1, ch_sign, sep = "|")) %>%
  group_by(episode) %>% mutate(n_in_episode = n(), w_ep = 1 / n_in_episode) %>% ungroup()

cat("Triad-events:", nrow(strict), " Unique episodes:", n_distinct(strict$episode), "\n")
cat("Triad-events per episode: ")
print(summary(strict$n_in_episode))
cat("Share of events in episodes with >=10 triad-events:", mean(strict$n_in_episode >= 10), "\n")

long <- strict %>%
  select(event_id, triad_id, episode, w_ep, n_in_episode, t1, ch_sign,
         sign12_t1, sign13_t1, sign23_t1,
         ul12, ul13, ul23, emb12, emb13, emb23, dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         cap_asymmetry_12, cap_asymmetry_13, cap_asymmetry_23,
         mutual_min_12, mutual_min_13, mutual_min_23, hhi_exp_12, hhi_exp_13, hhi_exp_23, hhio_exp_12, hhio_exp_13, hhio_exp_23,
         changed12, changed13, changed23) %>%
  pivot_longer(cols = c(ul12,ul13,ul23,emb12,emb13,emb23,dep_asymmetry_12,dep_asymmetry_13,dep_asymmetry_23,
                        cap_asymmetry_12,cap_asymmetry_13,cap_asymmetry_23, mutual_min_12,mutual_min_13,mutual_min_23,
                        hhi_exp_12,hhi_exp_13,hhi_exp_23, hhio_exp_12,hhio_exp_13,hhio_exp_23, changed12,changed13,changed23, sign12_t1,sign13_t1,sign23_t1),
               names_to = c(".value","tie"),
               names_pattern = "(ul|emb|dep_asymmetry|cap_asymmetry|mutual_min|hhio_exp|hhi_exp|changed|sign)_?(?:t1_)?(12|13|23)(?:_t1)?") %>%
  rename(sign_t1 = sign)
# (sign columns are named signXX_t1; the pattern above maps them to 'sign')
long <- long %>% mutate(initial_negative = as.numeric(sign_t1 == -1),
  z_ul = as.numeric(scale(ul)), z_dep = as.numeric(scale(dep_asymmetry)), z_emb = as.numeric(scale(emb)),
  z_cap = as.numeric(scale(cap_asymmetry)), z_mut = as.numeric(scale(mutual_min))) %>%
  filter(complete.cases(z_ul, z_dep, z_emb, z_cap, z_mut, initial_negative))
long <- long %>% group_by(event_id) %>% filter(n() == 3) %>% ungroup()
cat("Events in analysis:", n_distinct(long$event_id), "\n")

RHS <- "z_ul + z_dep + initial_negative + z_mut + z_cap + z_emb"
fit <- function(data, label, w = NULL, cl = "triad_id", rhs = RHS) {
  surv <- Surv(rep(1, nrow(data)), data$changed)
  f <- as.formula(paste("surv ~", rhs, "+ strata(event_id)"))
  m <- if (is.null(w)) coxph(f, data = data, method = "efron") else coxph(f, data = data, weights = data[[w]], method = "efron")
  ct <- coeftest(m, vcov = vcovCL(m, cluster = data[[cl]]))
  b <- ct["z_dep", ]
  cat(sprintf("%-62s b=%6.3f se=%5.3f p=%6.4f  events=%d  (episodes=%d)\n", label, b[1], b[2], b[4],
              n_distinct(data$event_id), n_distinct(data$episode)))
  invisible(m)
}

sink("results/event_level_and_predictive_output.txt", split = TRUE)
cat("\n===== (a)-(c) EVENT-LEVEL WEIGHTING AND TRANSITION SPLIT (M4) =====\n")
fit(long, "M4 baseline (triad-clustered)")
fit(long, "M4 baseline (clustered by episode)", cl = "episode")
fit(long, "M4 episode-weighted: each episode total weight 1 (cl: episode)", w = "w_ep", cl = "episode")
for (k in c(1, 2, 5, 10))
  fit(filter(long, n_in_episode <= k), sprintf("M4 restricted to episodes in <=%d triad-events (cl: episode)", k), cl = "episode")
cat("\n-- by transition type --\n")
R2 <- "z_ul + z_dep + z_mut + z_cap + z_emb"   # starting sign omitted: it perfectly separates within a transition type
fit(filter(long, ch_sign == -1), "dispute -> alliance (dispute endings), unweighted [no start-sign term]", cl = "episode", rhs = R2)
fit(filter(long, ch_sign == -1), "dispute -> alliance, episode-weighted", w = "w_ep", cl = "episode", rhs = R2)
fit(filter(long, ch_sign == 1),  "alliance -> dispute, unweighted", cl = "episode", rhs = R2)
fit(filter(long, ch_sign == 1),  "alliance -> dispute, episode-weighted", w = "w_ep", cl = "episode", rhs = R2)

# (d) out-of-sample predictive gain: conditional log-lik + top-1 accuracy,
#     folds keep episodes together (5-fold) and leave-one-decade-out-of-bloc proxy: leave-one-year-block-out
cat("\n===== (d) OUT-OF-SAMPLE PREDICTIVE GAIN =====\n")
RHS_NULL <- "z_ul + initial_negative + z_mut + z_cap + z_emb"
score <- function(m, data) {
  b <- coef(m); X <- model.matrix(as.formula(paste("~ 0 +", gsub("\\+ strata.*", "", paste(names(b), collapse = "+")))), data = data)
  eta <- as.numeric(X %*% b)
  data %>% mutate(eta = eta) %>% group_by(event_id) %>%
    summarise(ll = sum(eta * changed) - log(sum(exp(eta))), top1 = as.numeric(changed[which.max(eta)] == 1), .groups = "drop")
}
cv <- function(group_var, K, label) {
  g <- unique(long[[group_var]]); folds <- sample(rep(1:K, length.out = length(g))); names(folds) <- g
  long$fold <- folds[as.character(long[[group_var]])]
  res <- lapply(1:K, function(k) {
    tr <- filter(long, fold != k); te <- filter(long, fold == k)
    mk <- function(rhs) coxph(as.formula(paste("Surv(rep(1,nrow(tr)), changed) ~", rhs, "+ strata(event_id)")), data = tr, method = "efron")
    m1 <- mk(RHS); m0 <- mk(RHS_NULL)
    bind_rows(score(m1, te) %>% mutate(model = "full"), score(m0, te) %>% mutate(model = "no_dep_asym"))
  }) %>% bind_rows()
  s <- res %>% group_by(model) %>% summarise(mean_ll = mean(ll), top1 = mean(top1), n = n(), .groups = "drop")
  w <- res %>% select(event_id, model, ll, top1) %>% pivot_wider(names_from = model, values_from = c(ll, top1)) %>%
    mutate(dll = ll_full - ll_no_dep_asym, dt = top1_full - top1_no_dep_asym)
  cat("\n", label, "\n"); print(as.data.frame(s))
  cat(sprintf("  mean per-event log-lik gain = %.4f (bootstrap 95%% CI %.4f to %.4f); top-1 accuracy gain = %.4f\n",
              mean(w$dll), quantile(replicate(2000, mean(sample(w$dll, replace = TRUE))), .025),
              quantile(replicate(2000, mean(sample(w$dll, replace = TRUE))), .975), mean(w$dt)))
}
cv("episode", 10, "10-fold CV, folds keep each political episode together")
long$decade <- paste0(floor(long$t1 / 10) * 10)
cv("decade", 8, "8-fold CV, folds keep each decade together (temporal/bloc-era block)")

# (e) diversification heterogeneity
cat("\n===== (e) HETEROGENEITY BY DIVERSIFICATION OF THE MORE EXPOSED STATE =====\n")
ld <- long %>% filter(!is.na(hhi_exp)) %>% mutate(z_hhi = as.numeric(scale(hhi_exp)),
                                                  concentrated = as.numeric(hhi_exp > median(hhi_exp)))
fit(ld, "M4, subsample with HHI available", cl = "episode")
fit(filter(ld, concentrated == 1), "M4, exposed state's trade concentrated (HHI > median)", cl = "episode")
fit(filter(ld, concentrated == 0), "M4, exposed state's trade diversified (HHI <= median)", cl = "episode")
cat("\n-- concentration among the exposed state's OTHER partners (focal tie excluded) --\n")
lo2 <- long %>% filter(!is.na(hhio_exp)) %>% mutate(conc2 = as.numeric(hhio_exp > median(hhio_exp)), z_hhio = as.numeric(scale(hhio_exp)))
fit(filter(lo2, conc2 == 1), "M4, other-partner trade concentrated (> median)", cl = "episode")
fit(filter(lo2, conc2 == 0), "M4, other-partner trade diversified (<= median)", cl = "episode")
m_int2 <- fit(lo2, "M4 + z_dep x z_hhio", cl = "episode", rhs = paste(RHS, "+ z_hhio + z_dep:z_hhio"))
print(coeftest(m_int2, vcov = vcovCL(m_int2, cluster = lo2$episode))[c("z_dep","z_dep:z_hhio"), ])
m_int <- fit(ld, "M4 + z_dep x z_hhi (coef 'z_dep' is at mean HHI)", cl = "episode", rhs = paste(RHS, "+ z_hhi + z_dep:z_hhi"))
ct <- coeftest(m_int, vcov = vcovCL(m_int, cluster = ld$episode)); print(ct)
sink()
cat("\n33_event_level_and_predictive.R complete.\n")

# ---- (f) recurrence: same dyad changing in different years is not independent ----
long <- long %>% left_join(strict %>% select(event_id, ch_dyad), by = "event_id")
sink("results/event_level_and_predictive_output.txt", append = TRUE, split = TRUE)
cat("\n===== (f) RECURRENCE: cluster by changed dyad =====\n")
cat("Distinct changed dyads:", n_distinct(long$ch_dyad[long$changed == 1]), " episodes:", n_distinct(long$episode), "\n")
rec <- long %>% filter(changed == 1) %>% distinct(episode, ch_dyad) %>% count(ch_dyad)
cat("Episodes per dyad: "); print(summary(rec$n)); cat("Dyads with >=3 episodes:", sum(rec$n >= 3), "\n")
long$ch_dyad_ev <- ave(as.character(long$ch_dyad), long$event_id, FUN = function(x) x[1])
fit(long, "M4 unweighted, clustered by changed dyad", cl = "ch_dyad")
fit(long, "M4 episode-weighted, clustered by changed dyad", w = "w_ep", cl = "ch_dyad")
sink()
