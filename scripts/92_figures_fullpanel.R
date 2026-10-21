# 92: Figure 2 (coefficient map, full panel) and Figure 3 (predicted annual onset probability by exposure percentile, contiguous vs non-contiguous). Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(stringr); library(ggplot2); library(patchwork); library(sandwich) })
BLUE <- "#2a78d6"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GREY <- "#8a8985"; GRID <- "#e6e5e1"
th <- theme_minimal(base_size = 10, base_family = "Helvetica") + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
  axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2), plot.title = element_text(colour = INK, face = "bold", size = 10, hjust = 0), legend.position = "top", legend.justification = "left", legend.text = element_text(colour = INK2), plot.margin = margin(4, 8, 4, 4))
# ---- Figure 2 data parsed from outputs
p86 <- readLines("results/population_selection_output.txt"); p86 <- p86[seq_len(grep("^-- full panel, larger share x contiguity", p86)[1] - 1)]; rows <- list()
for (l in p86[1:12]) { m <- str_match(l, "^(dispute onset|dispute ending|alliance formation|alliance ending) (baseline|\\+ vol \\+ decade FE)\\s+larger share\\s+(-?[0-9.]+) \\[\\s*(-?[0-9.]+),\\s*(-?[0-9.]+)\\].*events=(\\d+)")
  if (!is.na(m[1, 1])) rows[[length(rows) + 1]] <- tibble(row = m[1, 2], kind = ifelse(m[1, 3] == "baseline", "Baseline controls", "Adjusted"), est = as.numeric(m[1, 4]), ci_lo = as.numeric(m[1, 5]), ci_hi = as.numeric(m[1, 6]), events = as.integer(m[1, 7])) }
A <- bind_rows(rows) %>% mutate(row = recode(row, `dispute onset` = "Dispute begins", `alliance formation` = "Alliance forms", `dispute ending` = "Dispute ends", `alliance ending` = "Alliance ends"))
A$row <- factor(A$row, levels = rev(c("Dispute begins", "Alliance forms", "Dispute ends", "Alliance ends")))
sv <- readLines("results/fullpanel_severity_initiator_output.txt"); spec <- NA; rows <- list()
for (l in sv) { if (str_detect(l, "spec: baseline")) spec <- "Baseline controls"; if (str_detect(l, "spec: \\+ z_vol")) spec <- "Adjusted"
  m <- str_match(l, "^\\s+(y_\\w+)\\s+(-?[0-9.]+) \\[\\s*(-?[0-9.]+),\\s*(-?[0-9.]+)\\] p=[0-9.]+ events (\\d+)"); if (!is.na(m[1, 1]) && !is.na(spec) && m[1, 2] != "y_any") rows[[length(rows) + 1]] <- tibble(key = m[1, 2], kind = spec, est = as.numeric(m[1, 3]), ci_lo = as.numeric(m[1, 4]), ci_hi = as.numeric(m[1, 5]), events = as.integer(m[1, 6])) }
B <- bind_rows(rows) %>% mutate(row = recode(key, y_thr = "Threat only", y_dis = "Display of force", y_force = "Use of force or war", y_fat = "Any fatalities"))
B$row <- factor(B$row, levels = rev(c("Threat only", "Display of force", "Use of force or war", "Any fatalities")))
all <- bind_rows(A %>% mutate(panel = "A"), B %>% select(-key) %>% mutate(panel = "B")) %>% mutate(kind = factor(kind, levels = c("Baseline controls", "Adjusted")), thin = row == "Threat only")
lims <- c(-0.4, 0.45)
mk <- function(df, title, legend = FALSE) ggplot(df, aes(est, row, fill = kind, alpha = thin)) + geom_vline(xintercept = 0, colour = INK2, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.55, colour = BLUE, position = position_dodge(width = 0.55)) + geom_point(shape = 21, size = 2.6, colour = BLUE, stroke = 0.9, position = position_dodge(width = 0.55)) +
  geom_text(data = df %>% filter(kind == "Baseline controls"), aes(x = lims[2], y = row, label = paste0("n = ", trimws(format(events, big.mark = ",")))), inherit.aes = FALSE, hjust = 1, vjust = -0.1, size = 2.9, colour = INK2, family = "Helvetica") +
  scale_fill_manual(values = c("Baseline controls" = BLUE, "Adjusted" = "white"), name = NULL, guide = if (legend) guide_legend(override.aes = list(shape = 21, alpha = 1, colour = BLUE)) else "none") + scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.55), guide = "none") +
  scale_x_continuous(limits = lims, breaks = seq(-0.4, 0.4, 0.2)) + labs(x = NULL, y = NULL, title = title) + th
fig2 <- mk(all %>% filter(panel == "A"), "A  Transitions, all dyad-years", TRUE) / (mk(all %>% filter(panel == "B"), "B  Severity of dispute onsets, by highest hostility level reached") + labs(x = "Larger trade share: log-odds per standard deviation (95% interval)")) + plot_layout(heights = c(4, 4))
ggsave("results/fig_exposure_map_full.pdf", fig2, width = 6.4, height = 5.2, device = cairo_pdf); ggsave("results/fig_exposure_map_full.png", fig2, width = 6.4, height = 5.2, dpi = 300)
# ---- Figure 3: predicted onset probability by percentile of the larger share
source_lines <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", source_lines)[1]; eval(parse(text = source_lines[1:(k - 1)]))
MU <- mean(d$mx); SDV <- sd(d$mx); x <- d %>% filter(m == 0) %>% mutate(onset_m = as.integer(m1 == 1), l_dur = log(dur_m))
m <- glm(onset_m ~ z_mx * contig + z_mn + z_cap + l_dur + l_same + l_diff + z_vol, x, family = binomial); V <- vcovCL(m, cluster = x$dyad); b <- coef(m); ok <- names(b)[!is.na(b)]; V <- V[ok, ok]; b <- b[ok]
base <- c(z_mn = mean(x$z_mn), z_cap = mean(x$z_cap), l_dur = mean(x$l_dur), l_same = mean(x$l_same), l_diff = mean(x$l_diff), z_vol = mean(x$z_vol))
qs <- c(seq(0.05, 0.95, 0.05), 0.97, 0.98, 0.99, 0.995); grid <- expand.grid(q = qs, contig = c(0, 1))
pred <- lapply(seq_len(nrow(grid)), function(i) { share <- unname(quantile(x$mx, grid$q[i])); z <- (share - MU) / SDV; ct <- grid$contig[i]; v <- c(`(Intercept)` = 1, z_mx = z, contig = ct, base, `z_mx:contig` = z * ct)[names(b)]; eta <- sum(b * v); se <- sqrt(drop(t(v) %*% V %*% v))
  tibble(q = grid$q[i], share = share, contig = ct, p = 1000 * plogis(eta), lo = 1000 * plogis(eta - 1.96 * se), hi = 1000 * plogis(eta + 1.96 * se)) }) %>% bind_rows() %>% mutate(grp = ifelse(contig == 1, "Contiguous pairs", "Non-contiguous pairs"))
write_csv(pred, "results/predicted_onset_by_percentile.csv")
brk <- unname(quantile(x$mx, c(0.5, 0.9, 0.99)))
fig3 <- ggplot(pred, aes(share, p, colour = grp, fill = grp, linetype = grp)) + geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c("Non-contiguous pairs" = BLUE, "Contiguous pairs" = GREY), name = NULL) + scale_fill_manual(values = c("Non-contiguous pairs" = BLUE, "Contiguous pairs" = GREY), name = NULL) + scale_linetype_manual(values = c("Non-contiguous pairs" = "solid", "Contiguous pairs" = "22"), name = NULL) +
  scale_x_log10(breaks = brk, labels = c("Median\n(0.02%)", "90th percentile\n(2.9%)", "99th percentile\n(23.1%)")) + scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20), labels = c("0.5", "1", "2", "5", "10", "20")) +
  labs(x = "Larger trade share (log scale; percentile of dyad-years at risk)", y = "Predicted onsets per 1,000 dyad-years (log scale)") + th + theme(panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3))
ggsave("results/fig_predicted_onset.pdf", fig3, width = 6.4, height = 3.9, device = cairo_pdf); ggsave("results/fig_predicted_onset.png", fig3, width = 6.4, height = 3.9, dpi = 300)
cat(sprintf("figures written; shares at p50/p90/p99: %s\n", paste(round(brk, 4), collapse = " ")))
