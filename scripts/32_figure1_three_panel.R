# ------------------------------------------------------------------
# 32: Figure 1 -- three panels, one inconsistent triad
#     (A-B alliance, A-C dispute, B-C alliance).
#
#     A: the allocation problem (structural balance alone: three equal
#        candidate sites of adjustment).
#     B: BEND (H1a) -- the highly asymmetric A-B tie is the most likely
#        tie to change sign.
#     C: BREAK (H1b) -- of the two eligible ties (alliances), the highly
#        asymmetric A-B tie is the LESS likely to end outright. A-C is a
#        dispute and is not eligible to be "the alliance that ends".
#
#     Illustrative values only; grayscale; base R. Edge width is
#     constant across panels so the only thing that varies is the
#     predicted likelihood (arrow weight), not an implied mechanism.
# ------------------------------------------------------------------


A <- c(0, 1); B <- c(-0.87, -0.5); C <- c(0.87, -0.5)
nodes <- list(A = A, B = B, C = C)
centroid <- (A + B + C) / 3
sgn  <- c(AB = 1, AC = -1, BC = 1)
ends <- list(AB = list(A, B), AC = list(A, C), BC = list(B, C))
asym <- c(AB = 0.57, AC = 0.15, BC = 0.29)

base <- function(title, sub) {
  plot(NA, xlim = c(-2.3, 2.3), ylim = c(-1.8, 2.0), axes = FALSE, xlab = "", ylab = "", asp = 1)
  mtext(title, side = 3, line = 2.6, cex = 1.45, font = 2)
  mtext(sub, side = 3, line = 0.9, cex = 0.95, font = 3)
}
draw_edges <- function(hl = NULL) {
  for (tn in names(sgn)) {
    p1 <- ends[[tn]][[1]]; p2 <- ends[[tn]][[2]]
    lines(c(p1[1], p2[1]), c(p1[2], p2[2]), col = if (!is.null(hl) && tn == hl) "black" else "grey55",
          lty = if (sgn[[tn]] > 0) "solid" else "22", lwd = if (!is.null(hl) && tn == hl) 3.4 else 2.2, lend = 1)
  }
  for (nm in names(nodes)) {
    p <- nodes[[nm]]
    points(p[1], p[2], pch = 21, bg = "white", col = "black", cex = 3.0, lwd = 1.5)
    text(p[1], p[2], nm, cex = 1.2, font = 2)
  }
}
outdir <- function(tn) {
  mid <- (ends[[tn]][[1]] + ends[[tn]][[2]]) / 2
  d <- mid - centroid; list(mid = mid, dir = d / sqrt(sum(d^2)))
}
arrow_at <- function(tn, strength, label, lab_font = 1, lab_col = "grey30", lab_cex = 1.05, show_asym = FALSE) {
  o <- outdir(tn)
  start <- o$mid + o$dir * 0.40
  tip   <- o$mid + o$dir * (0.62 + 0.45 * strength)
  arrows(start[1], start[2], tip[1], tip[2], lwd = 1.2 + 6 * strength, length = 0.06 + 0.08 * strength,
         angle = 20, col = grey(0.75 - 0.75 * strength), lend = 1)
  lp <- o$mid + o$dir * 1.28
  text(lp[1], lp[2], label, cex = lab_cex, font = lab_font, col = lab_col)
  if (show_asym) {
    vp <- o$mid + o$dir * 0.2
    text(vp[1], vp[2], sprintf("asym. %.2f", asym[[tn]]), cex = 0.8, col = "grey20")
  }
}

panel_a <- function() {
  base("A. The allocation problem", "Imbalance creates pressure to adjust, not where it lands")
  draw_edges()
  for (tn in names(sgn)) arrow_at(tn, 0.4, c(AB = "A-B changes", AC = "A-C changes", BC = "B-C changes")[[tn]])
}
panel_b <- function() {
  base("B. Bend (H1a)", "Asymmetric ties are more likely to change sign")
  draw_edges("AB")
  arrow_at("AB", 1.0, "A-B changes sign", 2, "black", 1.25, TRUE)
  arrow_at("AC", 0.15, "A-C", 1, "grey50", 1.0, TRUE)
  arrow_at("BC", 0.15, "B-C", 1, "grey50", 1.0, TRUE)
}
panel_c <- function() {
  base("C. Break (H1b)", "Asymmetric alliances are less likely to end outright")
  draw_edges("AB")
  arrow_at("AB", 0.15, "A-B ends (less likely)", 1, "grey40", 1.05, TRUE)
  arrow_at("BC", 0.9, "B-C ends (more likely)", 2, "black", 1.2, TRUE)
  o <- outdir("AC")
  text(o$mid[1] + o$dir[1] * 0.75, o$mid[2] + o$dir[2] * 0.75, "A-C: a dispute,\nnot eligible to end\nas an alliance", cex = 0.85, col = "grey40")
}
render <- function() {
  layout(matrix(1:3, nrow = 1))
  par(mar = c(1, 0.5, 6, 0.5), oma = c(3, 0, 2, 0))
  panel_a(); panel_b(); panel_c()
  mtext("Dependence asymmetry channels adjustment: asymmetric relationships are more likely to bend than to break.",
        side = 1, outer = TRUE, cex = 1.1, line = 1.2, font = 2)
}
png("results/fig1_triad_diagram_three_panel.png", width = 4200, height = 1700, res = 260); render(); dev.off()
pdf("results/fig1_triad_diagram_three_panel.pdf", width = 15.5, height = 6.3); render(); dev.off()
cat("32_figure1_three_panel.R complete.\n")
