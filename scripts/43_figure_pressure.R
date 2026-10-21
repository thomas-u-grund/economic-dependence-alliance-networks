# 43: Schematic of structural pressure on one tie (focal dispute A-C): one unbalanced and one balanced triad.

draw <- function() {
  par(mar = c(0.5, 0.5, 0.5, 0.5), family = "sans")
  plot(NA, xlim = c(-1.2, 9.2), ylim = c(-3.4, 3.4), axes = FALSE, xlab = "", ylab = "", asp = 1)
  P <- list(A = c(0, 0), C = c(4, 0), B = c(2, 2.3), D = c(2, -2.3))
  polygon(rbind(P$A, P$C, P$B), col = "#fde0dd", border = NA)
  polygon(rbind(P$A, P$C, P$D), col = "#deebf7", border = NA)
  seg <- function(a, b, pos = TRUE, lwd = 2.6) segments(P[[a]][1], P[[a]][2], P[[b]][1], P[[b]][2], lty = ifelse(pos, 1, 2), lwd = lwd, col = ifelse(pos, "#2166ac", "#b2182b"))
  seg("A", "B"); seg("B", "C"); seg("A", "D"); seg("A", "C", FALSE, 3.6); seg("C", "D", FALSE)
  for (n in names(P)) { points(P[[n]][1], P[[n]][2], pch = 21, bg = "white", cex = 4.6, lwd = 2); text(P[[n]][1], P[[n]][2], n, font = 2, cex = 1.1) }
  text(2, 0.32, "focal pair (dispute)", col = "#b2182b", cex = .85, font = 3)
  text(5.0, 1.45, "B is a common partner\n(ally of both)\nTriad A-B-C: (+, +, \u2212), unbalanced\nA dispute between A and C is\ninconsistent with it", cex = .85, col = "#67000d", adj = 0)
  text(5.0, -1.55, "D is a split partner\n(ally of A, adversary of C)\nTriad A-C-D: (+, \u2212, \u2212), balanced\nA dispute between A and C is\nconsistent with it", cex = .85, col = "#08306b", adj = 0)
  arrows(4.9, 1.0, 2.9, 0.9, length = .08, col = "#67000d"); arrows(4.9, -1.1, 2.9, -0.9, length = .08, col = "#08306b")
  legend("topleft", c("Alliance (+)", "Dispute (−)"), lty = c(1, 2), lwd = 2.4, col = c("#2166ac", "#b2182b"), bty = "n", cex = .85)
}
pdf("results/fig_pressure.pdf", width = 6.8, height = 5.2); draw(); dev.off()
png("results/fig_pressure.png", width = 1360, height = 1040, res = 200); draw(); dev.off()
