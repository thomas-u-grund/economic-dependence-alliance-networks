# ------------------------------------------------------------------
# 31: Joint "bend vs. break" model.
#
#     Theory: dependence asymmetry shapes the FORM of adjustment a tie
#     absorbs -- asymmetric ties are more likely to BEND (change sign)
#     and less likely to BREAK (alliance ends outright).
#
#     Design: every trade-valid triad-year in an unbalanced spell is a
#     risk set. Each of its three ties is at risk of two competing,
#     mutually exclusive transitions out of that year:
#         bend  : this tie changes sign (the single changed tie)
#         break : this tie, an alliance at t, disappears (single gone tie)
#     Stacking tie-years by transition type lets asymmetry carry a
#     type-specific coefficient; the bend-minus-break contrast is the
#     direct test of "bend rather than break". Cause-specific logits,
#     SEs clustered by triad (and by triad+dyad-sharing checked in 30).
#
#     Only single-tie transitions are used (as in scripts 05/18). Break
#     is defined only for ties that are alliances at t (eligibility).
# ------------------------------------------------------------------
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })



pyr <- readRDS("data/person_year_risk_table.rds")
tr  <- readRDS("data/triads_all_years.rds") %>%
  select(year, triad_id, sign12, sign13, sign23, tie12_emb, tie13_emb, tie23_emb)
signed_ties <- read_csv("data/signed_ties.csv", show_col_types = FALSE) %>%
  mutate(ccode_low = as.character(ccode_low), ccode_high = as.character(ccode_high)) %>%
  select(year, ccode_low, ccode_high) %>% mutate(present = TRUE)

# risk rows: trade-valid on all three ties
rr <- pyr %>% filter(complete.cases(dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23)) %>%
  mutate(row_id = row_number())
cat("Trade-valid triad-year risk rows:", nrow(rr), "\n")
print(table(rr$row_event))

# signs at t and t+1
rr <- rr %>%
  left_join(tr %>% rename(s12 = sign12, s13 = sign13, s23 = sign23,
                          emb12 = tie12_emb, emb13 = tie13_emb, emb23 = tie23_emb),
            by = c("triad_id", "year")) %>%
  left_join(tr %>% select(year, triad_id, n12 = sign12, n13 = sign13, n23 = sign23) %>%
              mutate(year = year - 1), by = c("triad_id", "year"))
# (n.. = signs at t+1, joined by shifting year back one)

present_next <- function(df, a, b, newcol) {
  df %>% mutate(lo = as.character(pmin(as.numeric(.data[[a]]), as.numeric(.data[[b]]))),
                hi = as.character(pmax(as.numeric(.data[[a]]), as.numeric(.data[[b]]))),
                yr1 = year + 1) %>%
    left_join(signed_ties, by = c("yr1" = "year", "lo" = "ccode_low", "hi" = "ccode_high")) %>%
    mutate(!!newcol := !is.na(present)) %>% select(-lo, -hi, -yr1, -present)
}
rr <- rr %>% present_next("node1", "node2", "pr12") %>%
  present_next("node1", "node3", "pr13") %>% present_next("node2", "node3", "pr23")

rr <- rr %>% mutate(
  ch12 = row_event == "flipped" & !is.na(n12) & s12 != n12,
  ch13 = row_event == "flipped" & !is.na(n13) & s13 != n13,
  ch23 = row_event == "flipped" & !is.na(n23) & s23 != n23,
  n_ch = ch12 + ch13 + ch23,
  gn12 = row_event == "dissolved" & !pr12, gn13 = row_event == "dissolved" & !pr13,
  gn23 = row_event == "dissolved" & !pr23, n_gn = gn12 + gn13 + gn23)

# keep: continues, single-change flips, single-gone dissolutions
rr <- rr %>% filter(row_event == "continues" | (row_event == "flipped" & n_ch == 1) |
                    (row_event == "dissolved" & n_gn == 1))

# unbalanced load per (year, edge), as in script 05
tr_all <- readRDS("data/triads_all_years.rds")
edge_long <- bind_rows(
  tr_all %>% transmute(year, balanced, lo = pmin(as.numeric(node1), as.numeric(node2)), hi = pmax(as.numeric(node1), as.numeric(node2))),
  tr_all %>% transmute(year, balanced, lo = pmin(as.numeric(node1), as.numeric(node3)), hi = pmax(as.numeric(node1), as.numeric(node3))),
  tr_all %>% transmute(year, balanced, lo = pmin(as.numeric(node2), as.numeric(node3)), hi = pmax(as.numeric(node2), as.numeric(node3)))) %>%
  group_by(year, lo, hi) %>% summarise(nub = sum(!balanced), .groups = "drop")
get_ul <- function(df, a, b, newcol) {
  df %>% mutate(lo = pmin(as.numeric(.data[[a]]), as.numeric(.data[[b]])),
                hi = pmax(as.numeric(.data[[a]]), as.numeric(.data[[b]]))) %>%
    left_join(edge_long, by = c("year", "lo", "hi")) %>%
    mutate(!!newcol := replace_na(nub - 1, 0)) %>% select(-lo, -hi, -nub)
}
rr <- rr %>% get_ul("node1", "node2", "ul12") %>% get_ul("node1", "node3", "ul13") %>% get_ul("node2", "node3", "ul23")

# long: one row per (triad-year, tie)
long <- rr %>%
  select(row_id, triad_id, year, duration, row_event,
         s12, s13, s23, emb12, emb13, emb23, ul12, ul13, ul23,
         dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
         ch12, ch13, ch23, gn12, gn13, gn23) %>%
  pivot_longer(cols = c(s12, s13, s23, emb12, emb13, emb23, ul12, ul13, ul23,
                        dep_asymmetry_12, dep_asymmetry_13, dep_asymmetry_23,
                        ch12, ch13, ch23, gn12, gn13, gn23),
               names_to = c(".value", "tie"),
               names_pattern = "(s|emb|ul|dep_asymmetry|ch|gn)_?(12|13|23)") %>%
  rename(sign_t = s, dep_asym = dep_asymmetry, bend = ch, brk = gn) %>%
  mutate(is_alliance = sign_t == 1, bend = as.integer(bend), brk = as.integer(brk))

cat("\nEvents in analysis sample: bend =", sum(long$bend), " break(all gone ties) =", sum(long$brk),
    " break(alliance gone) =", sum(long$brk[long$is_alliance]), "\n")

# standardize on the stacked tie-year sample (pooled across the two types)
z <- function(x) as.numeric(scale(x))
long <- long %>% mutate(z_dep = z(dep_asym), z_emb = z(emb), z_ul = z(ul), z_dur = z(log1p(duration)))

# ---- cause-specific logits ----------------------------------------
fit <- function(f, d) { m <- glm(f, data = d, family = binomial()); list(m = m, ct = coeftest(m, vcov = vcovCL(m, cluster = d$triad_id))) }

d_bend <- long                               # every tie is at risk of bending
d_brk  <- long %>% filter(is_alliance)       # only alliances can break
f_rhs  <- "z_dep + z_emb + z_ul + z_dur"
f_rhs_s <- paste(f_rhs, "+ is_alliance")

m_bend <- fit(as.formula(paste("bend ~", f_rhs_s)), d_bend)
m_brk  <- fit(as.formula(paste("brk ~",  f_rhs)),   d_brk)

# ---- joint: stack both types on the SAME eligible risk set (alliance ties) so the
#      contrast is not driven by who is eligible ----------------------
d_all <- long %>% filter(is_alliance)
stk <- bind_rows(d_all %>% mutate(type = "bend", y = bend), d_all %>% mutate(type = "break", y = brk)) %>%
  mutate(is_break = as.integer(type == "break"))
m_joint <- glm(y ~ is_break * (z_dep + z_emb + z_ul + z_dur), data = stk, family = binomial())
ct_joint <- coeftest(m_joint, vcov = vcovCL(m_joint, cluster = stk$triad_id))

# bend-on-alliance-ties only (like-for-like with break)
m_bend_all <- fit(as.formula(paste("bend ~", f_rhs)), d_all)

# genuinely non-overlap bends are handled by script 18; here, report share of bends that are
# on alliance ties at t (alliance->dispute) vs dispute ties at t (dispute->alliance)
cat("\nBends by starting sign: alliance ties =", sum(long$bend[long$is_alliance]),
    " | dispute ties =", sum(long$bend[!long$is_alliance]), "\n")

sink("results/joint_bend_break_output.txt")
cat("Trade-valid single-transition risk rows:", nrow(rr), "\n")
print(table(rr$row_event))
cat("\nBend events:", sum(long$bend), "(alliance->dispute:", sum(long$bend[long$is_alliance]),
    "| dispute->alliance:", sum(long$bend[!long$is_alliance]), ")\n")
cat("Break events (alliance ends):", sum(long$brk[long$is_alliance]), "\n")
cat("\n=== (1) BEND, all ties at risk (cause-specific logit, clustered by triad) ===\n"); print(m_bend$ct)
cat("\n=== (2) BREAK, alliance ties at risk ===\n"); print(m_brk$ct)
cat("\n=== (3) BEND restricted to alliance ties at t (like-for-like with break) ===\n"); print(m_bend_all$ct)
cat("\n=== (4) JOINT, alliance ties at risk of either; is_break x z_dep = contrast (break - bend) ===\n"); print(ct_joint)
b <- coef(m_joint); V <- vcovCL(m_joint, cluster = stk$triad_id)
k <- c("z_dep", "is_break:z_dep")
cat("\nbend coef (z_dep):", round(b["z_dep"], 3), " break coef:", round(b["z_dep"] + b["is_break:z_dep"], 3),
    "\ncontrast bend - break:", round(-b["is_break:z_dep"], 3), " SE:", round(sqrt(V["is_break:z_dep","is_break:z_dep"]), 3), "\n")
sink()
cat(readLines("results/joint_bend_break_output.txt"), sep = "\n")
saveRDS(list(m_bend = m_bend, m_brk = m_brk, m_bend_all = m_bend_all, m_joint = m_joint), "results/joint_bend_break.rds")

# ---- (5) triad fixed effects: compare tie-years within the same triad ----
suppressMessages(library(survival))
cl <- function(f, d) { m <- do.call(clogit, list(formula = f, data = d, method = "efron")); list(m = m, ct = coeftest(m, vcov = vcovCL(m, cluster = d$triad_id))) }
fe_bend_all <- cl(bend ~ z_dep + z_emb + z_ul + is_alliance + strata(triad_id), long)
fe_brk      <- cl(brk  ~ z_dep + z_emb + z_ul + strata(triad_id), d_brk)
fe_bend_alc <- cl(bend ~ z_dep + z_emb + z_ul + strata(triad_id), d_all)
stk_fe <- stk
fe_joint <- do.call(clogit, list(formula = y ~ is_break:z_dep + z_dep + z_emb + z_ul + is_break:z_emb + is_break:z_ul + is_break + strata(triad_id), data = stk_fe, method = "efron"))
ct_fe_joint <- coeftest(fe_joint, vcov = vcovCL(fe_joint, cluster = stk_fe$triad_id))
sink("results/joint_bend_break_output.txt", append = TRUE)
cat("\n\n######## TRIAD FIXED EFFECTS (conditional logit, strata = triad) ########\n")
cat("\n(5a) BEND, all ties\n"); print(fe_bend_all$ct)
cat("\n(5b) BREAK, alliance ties\n"); print(fe_brk$ct)
cat("\n(5c) BEND, alliance ties only\n"); print(fe_bend_alc$ct)
cat("\n(5d) JOINT, alliance ties; is_break:z_dep = contrast\n"); print(ct_fe_joint)
sink()
cat(readLines("results/joint_bend_break_output.txt")[-(1:75)], sep = "\n")
