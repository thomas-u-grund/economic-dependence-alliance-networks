# 18: Figure 2 (what exposure does: coefficient plot) and Figure 3 (embedding over time). Reads results/claims_table.csv, results/hostility_gradient_output.txt, results/era_bins_*.csv.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(stringr); library(ggplot2); library(patchwork) })

BLUE <- "#111111"; INK <- "#0b0b0b"; INK2 <- "#505050"; GREY <- "#8c8c8c"; GRID <- "#e6e6e6"
th <- theme_minimal(base_size = 10, base_family = "Helvetica") + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
  axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2), plot.title = element_text(colour = INK, face = "bold", size = 10, hjust = 0), legend.position = "top", legend.justification = "left",
  legend.text = element_text(colour = INK2), plot.margin = margin(4, 8, 4, 4))
# ---- Figure 2
ct <- read_csv("results/claims_table.csv", show_col_types = FALSE)
lab <- c("dispute onset" = "Dispute begins", "alliance formation" = "Alliance forms", "dispute ending" = "Dispute ends", "alliance ending" = "Alliance ends", "dispute ending, unbalanced (+,+,-) triads" = "Dispute ends,\ninside unbalanced triads")
pa <- ct %>% filter(spec %in% c("baseline", "+ volume + decade FE", "spec D (controls, allied)", "spec D + pact FE")) %>%
  mutate(row = lab[outcome], kind = ifelse(spec %in% c("baseline", "spec D (controls, allied)"), "Baseline controls", "Adjusted"), n = paste0(events)) %>% select(row, kind, est, ci_lo, ci_hi, events)
pa$row <- factor(pa$row, levels = rev(unname(lab))); pa$panel <- "A  Transitions in embedded pairs"
hl <- readLines("results/hostility_gradient_output.txt"); spec <- NA; out <- list()
for (l in hl) { if (str_detect(l, "spec: baseline")) spec <- "Baseline controls"; if (str_detect(l, "spec: \\+ z_vol")) spec <- "Adjusted"
  m <- str_match(l, "^\\s+(y_\\w+)\\s+est\\s+(-?[0-9.]+) .*CI \\[\\s*(-?[0-9.]+),\\s*(-?[0-9.]+)\\].*events (\\d+)")
  if (!is.na(m[1, 1]) && !is.na(spec) && m[1, 2] != "y_any") out[[length(out) + 1]] <- tibble(key = m[1, 2], kind = spec, est = as.numeric(m[1, 3]), ci_lo = as.numeric(m[1, 4]), ci_hi = as.numeric(m[1, 5]), events = as.integer(m[1, 6])) }
sev <- bind_rows(out) %>% mutate(row = recode(key, y_thr = "Threat only", y_dis = "Display of force", y_force = "Use of force or war", y_fat = "Any fatalities"))
sev$row <- factor(sev$row, levels = rev(c("Threat only", "Display of force", "Use of force or war", "Any fatalities"))); sev$panel <- "B  Severity of dispute onsets"
all <- bind_rows(pa, sev %>% select(row, kind, est, ci_lo, ci_hi, events, panel)) %>% mutate(kind = factor(kind, levels = c("Baseline controls", "Adjusted")), thin = row == "Threat only")
lims <- c(-1.0, 0.6)
mk <- function(df, title, legend = FALSE) ggplot(df, aes(est, row, shape = kind, fill = kind, alpha = thin)) +
  geom_vline(xintercept = 0, colour = INK2, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0, linewidth = 0.55, colour = BLUE, position = position_dodge(width = 0.55)) +
  geom_point(size = 2.6, colour = BLUE, stroke = 0.9, position = position_dodge(width = 0.55)) +
  geom_text(data = df %>% filter(kind == "Baseline controls"), aes(x = lims[2], y = row, label = paste0("n = ", events)), inherit.aes = FALSE, hjust = 1, vjust = -0.1, size = 2.9, colour = INK2, family = "Helvetica") +
  scale_shape_manual(values = c("Baseline controls" = 21, "Adjusted" = 21), guide = "none") +
  scale_fill_manual(values = c("Baseline controls" = BLUE, "Adjusted" = "white"), name = NULL, labels = c("Baseline controls", "Adjusted"), guide = if (legend) guide_legend(override.aes = list(shape = 21, alpha = 1, colour = BLUE)) else "none") +
  scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.45), guide = "none") + scale_x_continuous(limits = lims, breaks = seq(-1, 0.4, 0.2)) + labs(x = NULL, y = NULL, title = title) + th
fa <- mk(all %>% filter(str_starts(panel, "A")), "A  Transitions in embedded pairs", TRUE); fb <- mk(all %>% filter(str_starts(panel, "B")), "B  Severity of dispute onsets, by highest hostility level reached") + labs(x = "Larger trade share: log-odds per standard deviation (95% interval)")
fig2 <- fa / fb + plot_layout(heights = c(5, 4))
ggsave("results/fig_exposure_map.pdf", fig2, width = 6.4, height = 5.6, device = cairo_pdf); ggsave("results/fig_exposure_map.png", fig2, width = 6.4, height = 5.6, dpi = 300)
# ---- Figure 3
r <- read_csv("results/era_bins_rates.csv", show_col_types = FALSE); sh <- read_csv("results/era_bins_shares.csv", show_col_types = FALSE)
lv <- c("1870-89", "1890-1909", "1910-29", "1930-49", "1950-69", "1970-89", "1990-2010"); r$per <- factor(r$per, lv); sh$per <- factor(sh$per, lv); r$grp <- ifelse(r$tri, "Embedded pairs", "Other pairs")
f3a <- ggplot(r, aes(per, rate, group = grp, colour = grp, shape = grp, linetype = grp)) + geom_linerange(aes(ymin = lo, ymax = hi), linewidth = 0.5, position = position_dodge(width = 0.25), show.legend = FALSE) +
  geom_line(linewidth = 0.7, position = position_dodge(width = 0.25)) + geom_point(size = 2.4, position = position_dodge(width = 0.25), fill = "white", stroke = 0.9) +
  scale_colour_manual(values = c("Embedded pairs" = BLUE, "Other pairs" = GREY), name = NULL) + scale_shape_manual(values = c("Embedded pairs" = 19, "Other pairs" = 22), name = NULL) + scale_linetype_manual(values = c("Embedded pairs" = "solid", "Other pairs" = "22"), name = NULL) +
  scale_y_log10(breaks = c(1, 3, 10, 30), labels = c("1", "3", "10", "30")) + labs(x = NULL, y = "Onsets per 1,000 dyad-years (log scale)", title = "A  Dispute onset rate") + th + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3))
s <- sh %>% pivot_longer(c(pct_dyad_years_embedded, pct_onsets_embedded), names_to = "k", values_to = "pct") %>% mutate(k = recode(k, pct_dyad_years_embedded = "Share of dyad-years", pct_onsets_embedded = "Share of onsets"))
f3b <- ggplot(s, aes(per, pct, group = k, colour = k, shape = k, linetype = k)) + geom_line(linewidth = 0.7) + geom_point(size = 2.4, fill = "white", stroke = 0.9) +
  scale_colour_manual(values = c("Share of onsets" = BLUE, "Share of dyad-years" = GREY), name = NULL) + scale_shape_manual(values = c("Share of onsets" = 19, "Share of dyad-years" = 22), name = NULL) + scale_linetype_manual(values = c("Share of onsets" = "solid", "Share of dyad-years" = "22"), name = NULL) +
  scale_y_continuous(limits = c(0, 70), breaks = seq(0, 60, 20)) + labs(x = NULL, y = "Per cent that are embedded", title = "B  Embedded pairs' share") + th + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3))
fig3 <- f3a / f3b + plot_layout(heights = c(5, 3.4))
ggsave("results/fig_eras.pdf", fig3, width = 6.4, height = 5.4, device = cairo_pdf); ggsave("results/fig_eras.png", fig3, width = 6.4, height = 5.4, dpi = 300)
cat("figures written\n")
