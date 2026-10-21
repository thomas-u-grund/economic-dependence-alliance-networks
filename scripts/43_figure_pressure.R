# 43: Schematic of common and split partners around a focal pair A-C in a dispute (two panels).
setwd("/Users/tgrund/FILES_NEW/RESEARCH/2026 Trade")
pos <- "#2166ac"; neg <- "#b2182b"
panel <- function(title, third, s1, s2, tint, verdict, why, expect) {
  plot(NA, xlim = c(-1.0, 5.0), ylim = c(-4.3, 3.9), axes = FALSE, xlab = "", ylab = "", asp = 1)
  P <- list(A = c(0, 0), C = c(4, 0), X = c(2, 3))
  polygon(rbind(P$A, P$C, P$X), col = tint, border = NA)
  seg <- function(a, b, s, lwd = 3.4) { u <- P[[b]] - P[[a]]; u <- u / sqrt(sum(u^2)); r <- 0.42; p0 <- P[[a]] + r * u; p1 <- P[[b]] - r * u
    segments(p0[1], p0[2], p1[1], p1[2], lty = ifelse(s > 0, 1, 2), lwd = lwd, col = ifelse(s > 0, pos, neg), lend = 1) }
  seg("A", "C", -1, 4.6); seg("A", "X", s1); seg("X", "C", s2)
  for (n in c("A", "C", "X")) { points(P[[n]][1], P[[n]][2], pch = 21, bg = "white", cex = 6.2, lwd = 2.4); text(P[[n]][1], P[[n]][2], ifelse(n == "X", third, n), font = 2, cex = 1.35) }
  text(2, -0.62, "focal pair in a dispute", col = neg, cex = 1.0, font = 3)
  text(2, 4.05, title, font = 2, cex = 1.5, xpd = NA)
  text(2, -1.55, verdict, cex = 1.2, font = 2, xpd = NA)
  text(2, -2.25, why, cex = 1.0, xpd = NA)
  segments(-0.8, -2.85, 4.8, -2.85, col = "grey70", lwd = 1)
  text(2, -3.45, expect, cex = 1.05, font = 3, xpd = NA)
}
draw <- function() {
  par(mfrow = c(1, 2), mar = c(0.3, 0.3, 0.6, 0.3), family = "sans")
  panel("Common partner", "B", 1, 1, "#e8eef7", "Triad (+, +, −): unbalanced", "A and C are both allied with B,\nso a dispute between them is inconsistent", "Balance theory expects the dispute to end")
  panel("Split partner", "D", 1, -1, "#f7e8ea", "Triad (+, −, −): balanced", "D is allied with A and in dispute with C,\nso a dispute between A and C is consistent", "Balance theory expects disputes to begin")
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0)); plot.new()
  legend("bottom", c("Alliance (+)", "Dispute (−)"), lty = c(1, 2), lwd = 3, col = c(pos, neg), bty = "n", horiz = TRUE, cex = 1.05, inset = -0.005, xpd = NA)
}
pdf("results/fig_pressure.pdf", width = 8.4, height = 4.9); draw(); dev.off()
png("results/fig_pressure.png", width = 1680, height = 980, res = 200); draw(); dev.off()
