# 102: Geographic descriptives (raw onset rates by exposure tier, contiguity and distance) and map of the highly exposed non-contiguous pairs. Run with ALLIANCE_DEF=defense.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(ggplot2); library(patchwork); library(sf); library(cshapes); library(rnaturalearth) })
src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]; eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(onset_m = as.integer(m1 == 1)) %>% filter(year >= 1886, year <= 2011)
yrs <- sort(unique(d$year))
cap <- bind_rows(lapply(yrs, function(y) { m <- distmatrix(as.Date(paste0(y, "-01-01")), type = "capdist", useGW = FALSE); ids <- as.numeric(rownames(m)); ix <- which(upper.tri(m), arr.ind = TRUE)
  tibble(year = y, lo = pmin(ids[ix[,1]], ids[ix[,2]]), hi = pmax(ids[ix[,1]], ids[ix[,2]]), capdist = m[ix]) }))
x <- d %>% left_join(cap, by = c("year", "lo", "hi")) %>% filter(m == 0, !is.na(capdist))
BLUE <- "#2a78d6"; ORANGE <- "#c8651b"; INK <- "#0b0b0b"; INK2 <- "#52514e"; GRID <- "#e6e5e1"
th <- theme_minimal(base_size = 10, base_family = "Helvetica") + theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
  axis.text = element_text(colour = INK2), axis.title = element_text(colour = INK2), plot.title = element_text(colour = INK, face = "bold", size = 10, hjust = 0), legend.position = "top", legend.justification = "left", plot.margin = margin(4, 8, 4, 4))
bt <- function(k, n) { r <- binom.test(k, n)$conf.int; c(1000 * k / n, 1000 * r[1], 1000 * r[2]) }
sink("results/geography_descriptives_output.txt", split = TRUE)
cat(sprintf("At-risk dyad-years with capital distance (1886-2010): %d; onsets %d; contiguous %d (%d onsets)\n", nrow(x), sum(x$onset_m), sum(x$contig == 1), sum(x$onset_m[x$contig == 1])))
q <- quantile(x$mx, c(.5, .9, .99)); cat(sprintf("Larger-share cut points (panel-wide): p50 %.4f, p90 %.4f, p99 %.4f\n", q[1], q[2], q[3]))
x <- x %>% mutate(tier = cut(mx, c(-Inf, q, Inf), labels = c("Below median", "Median to 90th", "90th to 99th", "Top 1%")), geo = ifelse(contig == 1, "Contiguous", "Not contiguous"))
A <- x %>% group_by(geo, tier) %>% summarise(n = n(), k = sum(onset_m), .groups = "drop") %>% rowwise() %>% mutate(r = list(bt(k, n)), rate = r[1], lo = r[2], hi = r[3]) %>% ungroup()
cat("\n(a) Onsets per 1,000 dyad-years at risk by exposure tier\n"); print(as.data.frame(A %>% select(geo, tier, n, k, rate, lo, hi)), digits = 3)
B <- x %>% filter(contig == 0) %>% mutate(band = cut(capdist, c(0, 2000, 5000, 9000, Inf), labels = c("<2,000", "2,000-5,000", "5,000-9,000", ">9,000")), expo = ifelse(mx > q[2], "Above 90th percentile", "At or below 90th percentile")) %>%
  group_by(band, expo) %>% summarise(n = n(), k = sum(onset_m), .groups = "drop") %>% rowwise() %>% mutate(r = list(bt(k, n)), rate = r[1], lo = r[2], hi = r[3]) %>% ungroup()
cat("\n(b) Non-contiguous pairs: onsets per 1,000 by capital-distance band (km) and exposure\n"); print(as.data.frame(B %>% select(band, expo, n, k, rate, lo, hi)), digits = 3)
pa <- ggplot(A, aes(tier, rate, colour = geo, group = geo)) + geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.4), size = 0.35) + scale_colour_manual(values = c(Contiguous = ORANGE, `Not contiguous` = BLUE), name = NULL) +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50), labels = function(v) format(v, trim = TRUE, drop0trailing = TRUE)) + labs(x = "Larger trade share, percentile band", y = "Onsets per 1,000 dyad-years", title = "A  Exposure tier and contiguity") + th
pb <- ggplot(B, aes(band, rate, colour = expo, group = expo)) + geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(0.4), size = 0.35) + scale_colour_manual(values = c(`At or below 90th percentile` = "#8a8985", `Above 90th percentile` = BLUE), name = NULL) +
  scale_y_log10(breaks = c(0.5, 1, 2, 5), labels = function(v) format(v, trim = TRUE, drop0trailing = TRUE)) + labs(x = "Distance between capitals (km), non-contiguous pairs", y = "Onsets per 1,000 dyad-years", title = "B  Distance and exposure") + th
ggsave("results/fig_geography_descriptives.png", pa / pb, width = 6.4, height = 6, dpi = 300, bg = "white"); ggsave("results/fig_geography_descriptives.pdf", pa / pb, width = 6.4, height = 6, device = cairo_pdf)
# ---- tail and map
nc <- x %>% filter(contig == 0); t99 <- quantile(nc$mx, .99); tail <- nc %>% filter(mx > t99)
cat(sprintf("\n(c) Top 1%% of the larger share among non-contiguous at-risk dyad-years: > %.4f; %d dyad-years in %d dyads; onsets %d\n", t99, nrow(tail), n_distinct(tail$dyad), sum(tail$onset_m)))
cat(sprintf("Median capital distance: tail %.0f km; other non-contiguous %.0f km; tail above 5,000 km: %.0f%%; other: %.0f%%\n", median(tail$capdist), median(nc$capdist[nc$mx <= t99]), 100 * mean(tail$capdist > 5000), 100 * mean(nc$capdist[nc$mx <= t99] > 5000)))
caps <- as.data.frame(cshp(date = NA, useGW = FALSE)) %>% arrange(cowcode, desc(end)) %>% group_by(cowcode) %>% slice(1) %>% ungroup() %>% select(cowcode, country_name, caplong, caplat)
ed <- tail %>% group_by(lo, hi) %>% summarise(ny = n(), on = sum(onset_m), .groups = "drop") %>% left_join(caps, by = c("lo" = "cowcode")) %>% rename(n1 = country_name, x1 = caplong, y1 = caplat) %>% left_join(caps, by = c("hi" = "cowcode")) %>% rename(n2 = country_name, x2 = caplong, y2 = caplat) %>% filter(!is.na(x1), !is.na(x2))
cat(sprintf("Dyads drawn: %d of %d\n", nrow(ed), n_distinct(tail$dyad)))
st <- bind_rows(tail %>% transmute(cowcode = lo), tail %>% transmute(cowcode = hi)) %>% count(cowcode, name = "ny") %>% left_join(caps, by = "cowcode") %>% arrange(desc(ny)); cat("\nStates by tail dyad-years:\n"); print(head(as.data.frame(st %>% select(country_name, ny)), 12))
PROJ <- "+proj=robin +lon_0=0"
gc <- function(x1, y1, x2, y2, n = 60) { r <- pi / 180; v <- function(lo, la) c(cos(la * r) * cos(lo * r), cos(la * r) * sin(lo * r), sin(la * r)); a <- v(x1, y1); b <- v(x2, y2); w <- acos(pmin(1, sum(a * b)))
  t <- seq(0, 1, length.out = n); P <- sapply(t, function(u) if (w < 1e-8) a else (sin((1 - u) * w) * a + sin(u * w) * b) / sin(w)); lon <- atan2(P[2, ], P[1, ]) / r; lat <- asin(pmax(-1, pmin(1, P[3, ]))) / r
  brk <- c(0, which(abs(diff(lon)) > 180), length(lon)); st_multilinestring(lapply(seq_len(length(brk) - 1), function(i) cbind(lon[(brk[i] + 1):brk[i + 1]], lat[(brk[i] + 1):brk[i + 1]]))) }
lines <- st_transform(st_sfc(lapply(seq_len(nrow(ed)), function(i) gc(ed$x1[i], ed$y1[i], ed$x2[i], ed$y2[i])), crs = 4326), PROJ)
el <- st_sf(ed, geometry = lines) %>% mutate(grp = ifelse(on > 0, "At least one onset", "No onset"))
wd <- ne_countries(scale = 110, returnclass = "sf") %>% filter(!name %in% "Antarctica") %>% st_transform(PROJ)
pts <- st_transform(st_as_sf(st %>% filter(!is.na(caplong)), coords = c("caplong", "caplat"), crs = 4326), PROJ)
lab <- pts %>% filter(country_name %in% c("United States of America", "United Kingdom", "Japan")); lab$nm <- recode(lab$country_name, `United States of America` = "United States")
mp <- ggplot() + geom_sf(data = wd, fill = "#f1f0ec", colour = "#d6d4cd", linewidth = 0.15) + geom_sf(data = el %>% filter(grp == "No onset"), aes(linewidth = ny), colour = "#8fb0de", alpha = 0.35) +
  geom_sf(data = el %>% filter(grp != "No onset"), aes(linewidth = ny), colour = ORANGE, alpha = 0.85) + geom_sf(data = pts, aes(size = ny), colour = INK, alpha = 0.7, shape = 16) +
  geom_sf_text(data = lab, aes(label = nm), size = 2.6, colour = INK, nudge_y = 900000, fontface = "bold", family = "Helvetica") +
  scale_linewidth(range = c(0.15, 1.4), guide = "none") + scale_size(range = c(0.6, 4.5), guide = "none") + theme_void(base_family = "Helvetica") + theme(plot.margin = margin(2, 2, 2, 2), plot.background = element_rect(fill = "white", colour = NA))
ggsave("results/fig_geography_map.png", mp, width = 6.8, height = 3.6, dpi = 300, bg = "white"); ggsave("results/fig_geography_map.pdf", mp, width = 6.8, height = 3.8, device = cairo_pdf)
sink()
