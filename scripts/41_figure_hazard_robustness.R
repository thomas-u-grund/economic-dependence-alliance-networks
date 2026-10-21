# 41: Coefficient plot of the larger-share (maximum exposure) effect on sign change across specifications.
# Estimates/SEs transcribed from results/hazard_*_output.txt (scripts 35, 36, 39) under the defence-pact definition; the four-type row comes from a run without ALLIANCE_DEF.

s <- data.frame(
 spec = c("Baseline (defence pacts)","All four treaty types","+ Bloc and decade FE","+ Lagged 1y, bloc, decade FE","Larger share lagged 3 years","Larger share lagged 5 years",
          "Complementary log-log","Cox hazard","+ Prior spells of dyad","Dyadic-robust (shared-state) SEs",
          "Hostility level >= 4","Change within 3 years","Dyad fixed effects","Dyad and decade FE"),
 dB = c(.405,.315,.477,.413,.301,.215,.198,.174,.428,.405,.378,.403,.323,.429),
 dS = c(.118,.090,.117,.122,.142,.143,.141,.083,.120,.120,.143,.145,.232,.264),
 aB = c(.347,.317,.368,.372,.356,.320,.334,.336,.153,.347,.371,.392,.172,.118),
 aS = c(.056,.049,.056,.056,.058,.053,.054,.053,.059,.043,.065,.063,.105,.108))
write.csv(s, "results/figure_hazard_robustness_data.csv", row.names = FALSE)
draw <- function() {
  par(mar = c(4, 15.5, 2.2, 1), mgp = c(2.2, .7, 0), family = "sans")
  n <- nrow(s); y <- rev(seq_len(n))
  plot(NA, xlim = c(-0.3, 1.0), ylim = c(0.5, n + 1.6), yaxt = "n", xlab = "Log-odds per SD of the larger share (95% CI)", ylab = "", bty = "l")
  abline(v = 0, lty = 3, col = "grey40"); abline(h = y, col = "grey92")
  axis(2, at = y, labels = s$spec, las = 1, tick = FALSE, cex.axis = .85)
  off <- .17
  segments(s$dB - 1.96*s$dS, y + off, s$dB + 1.96*s$dS, y + off, col = "#b2182b", lwd = 1.6); points(s$dB, y + off, pch = 19, col = "#b2182b", cex = .9)
  segments(s$aB - 1.96*s$aS, y - off, s$aB + 1.96*s$aS, y - off, col = "#2166ac", lwd = 1.6); points(s$aB, y - off, pch = 15, col = "#2166ac", cex = .9)
  legend("top", horiz = TRUE, c("Dispute sign change (ending between allies)", "Alliance sign change (dispute begins)"), pch = c(19, 15), col = c("#b2182b", "#2166ac"), bty = "n", cex = .8)
}
pdf("results/fig_hazard_robustness.pdf", width = 7.4, height = 5.4); draw(); dev.off()
png("results/fig_hazard_robustness.png", width = 1480, height = 1080, res = 200); draw(); dev.off()
