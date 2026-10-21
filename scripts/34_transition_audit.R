# 34: Descriptive audit of how all unbalanced-triad spells end, with counts of
#     unique dyadic episodes (changed/gone dyad x year x prior sign) per type.
suppressMessages({ library(dplyr); library(readr); library(tidyr) })

sb <- readRDS("data/spell_bounds.rds")
tr <- readRDS("data/triads_all_years.rds") %>% select(year, triad_id, node1, node2, node3, sign12, sign13, sign23)
st <- read_csv("data/signed_ties.csv", show_col_types = FALSE) %>%
  mutate(lo = as.numeric(ccode_low), hi = as.numeric(ccode_high)) %>% select(year, lo, hi) %>% mutate(present = TRUE)
tt <- read_csv("data/tie_trade.csv", show_col_types = FALSE) %>% transmute(year, lo = as.numeric(ccode_low), hi = as.numeric(ccode_high), da = dep_asymmetry)

d <- sb %>% left_join(tr, by = c("triad_id", "t1" = "year")) %>%
  mutate(across(c(node1,node2,node3), as.numeric))
pr <- function(df, a, b, nm) df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]]), y1 = t1 + 1) %>%
  left_join(st, by = c("y1" = "year", "lo", "hi")) %>% mutate(!!paste0("pr", nm) := !is.na(present)) %>% select(-lo, -hi, -y1, -present)
ta <- function(df, a, b, nm) df %>% mutate(lo = pmin(.data[[a]], .data[[b]]), hi = pmax(.data[[a]], .data[[b]])) %>%
  left_join(tt, by = c("t1" = "year", "lo", "hi")) %>% mutate(!!paste0("da", nm) := da) %>% select(-lo, -hi, -da)
d <- d %>% pr("node1","node2","12") %>% pr("node1","node3","13") %>% pr("node2","node3","23") %>%
  ta("node1","node2","12") %>% ta("node1","node3","13") %>% ta("node2","node3","23")
# next-year sign for flips
tn <- tr %>% transmute(triad_id, t1 = year - 1, n12 = sign12, n13 = sign13, n23 = sign23)
d <- d %>% left_join(tn, by = c("triad_id", "t1"))
d <- d %>% mutate(
  ch12 = event == "flipped" & !is.na(n12) & sign12 != n12, ch13 = event == "flipped" & !is.na(n13) & sign13 != n13,
  ch23 = event == "flipped" & !is.na(n23) & sign23 != n23,
  gn12 = event == "dissolved" & !pr12, gn13 = event == "dissolved" & !pr13, gn23 = event == "dissolved" & !pr23,
  nch = ch12 + ch13 + ch23, ngn = gn12 + gn13 + gn23,
  type = case_when(
    event == "flipped" & nch == 1 & ((ch12 & sign12 == -1) | (ch13 & sign13 == -1) | (ch23 & sign23 == -1)) ~ "sign change: dispute -> alliance",
    event == "flipped" & nch == 1 ~ "sign change: alliance -> dispute",
    event == "flipped" ~ "sign change: multiple ties",
    event == "dissolved" & ngn == 1 & ((gn12 & sign12 == 1) | (gn13 & sign13 == 1) | (gn23 & sign23 == 1)) ~ "dissolution: alliance ends",
    event == "dissolved" & ngn == 1 ~ "dissolution: dispute ends (no alliance)",
    event == "dissolved" ~ "dissolution: multiple ties",
    TRUE ~ paste("other:", event)),
  gdy = case_when(ch12 | gn12 ~ paste(pmin(node1,node2), pmax(node1,node2), sep="-"),
                  ch13 | gn13 ~ paste(pmin(node1,node3), pmax(node1,node3), sep="-"),
                  ch23 | gn23 ~ paste(pmin(node2,node3), pmax(node2,node3), sep="-")),
  gsign = case_when(ch12 | gn12 ~ sign12, ch13 | gn13 ~ sign13, ch23 | gn23 ~ sign23),
  episode = ifelse(is.na(gdy), NA, paste(gdy, t1, gsign, sep = "|")),
  tradevalid = !is.na(da12) & !is.na(da13) & !is.na(da23))
sink("results/transition_audit_output.txt", split = TRUE)
cat("Spells:", nrow(d), "\n\n")
a <- d %>% group_by(type) %>% summarise(spells = n(), episodes = n_distinct(episode, na.rm = TRUE),
  spells_trade_valid = sum(tradevalid), episodes_trade_valid = n_distinct(episode[tradevalid], na.rm = TRUE), .groups = "drop") %>% arrange(desc(spells))
print(as.data.frame(a))
cat("\nAlliance-ends events, trade-valid: triads per episode\n")
al <- d %>% filter(type == "dissolution: alliance ends", tradevalid)
print(summary(as.numeric(table(al$episode)))); print(table(al$ngn))
cat("\nWithin alliance-ends (trade valid): eligible alliances = positive ties at t1. Gone alliance more/less asymmetric than the other alliance?\n")
al <- al %>% rowwise() %>% mutate(
  s = list(c(sign12, sign13, sign23)), a = list(c(da12, da13, da23)), g = list(c(gn12, gn13, gn23)),
  gone_a = a[g][1], oth = list(a[!g & s == 1]), other_a = ifelse(length(oth) >= 1, oth[1], NA)) %>% ungroup()
cat("events with exactly one other alliance:", sum(!is.na(al$other_a)), "\n")
dd <- al %>% filter(!is.na(other_a))
cat("gone alliance less asymmetric in", mean(dd$gone_a < dd$other_a), "of events; episode-level (mean within episode):",
    mean((dd %>% group_by(episode) %>% summarise(m = mean(gone_a < other_a)))$m), "\n")
ep <- dd %>% group_by(episode) %>% summarise(m = mean(gone_a < other_a), n = n())
cat("episodes:", nrow(ep), " episodes where gone alliance is less asymmetric (majority):", sum(ep$m > .5), " more:", sum(ep$m < .5), " tie:", sum(ep$m == .5), "\n")
print(binom.test(sum(ep$m > .5), sum(ep$m != .5)))
cat("\nDispute-ends-without-alliance: triads per episode\n"); x <- d %>% filter(type == "dissolution: dispute ends (no alliance)", tradevalid); print(summary(as.numeric(table(x$episode))))
sink()
