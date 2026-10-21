# 51: Bend vs break for alliances: dispute begins inside the alliance (A->B) vs alliance ends (A->D). Stacked cause-specific logit contrast.
suppressMessages({ library(dplyr); library(readr); library(tidyr); library(sandwich); library(lmtest) })

src <- readLines("scripts/49_layered_hazards.R"); k <- grep("^sink\\(", src)[1]
eval(parse(text = src[1:(k - 1)]))
d <- d %>% mutate(l_dur = log(dur_a))
A <- d %>% filter(a == 1, m == 0) %>% mutate(bend = as.integer(a1 == 1 & m1 == 1), brk = as.integer(a1 == 0 & m1 == 0)) %>% filter(!(a1 == 0 & m1 == 1))
cat("alliance-years, no dispute:", nrow(A), " bend events:", sum(A$bend), " break events:", sum(A$brk), "\n")
stack <- function(x, rhs) {
  xb <- x %>% mutate(y = bend, out = "bend"); xk <- x %>% mutate(y = brk, out = "brk")
  s <- bind_rows(xb, xk) %>% mutate(ob = as.integer(out == "bend"))
  vars <- strsplit(rhs, " \\+ ")[[1]]
  f <- as.formula(paste("y ~ ob +", paste(vars, collapse = " + "), "+", paste0("ob:", vars, collapse = " + ")))
  m <- glm(f, data = s, family = binomial); V <- vcovCL(m, cluster = s$dyad); b <- coef(m)
  nm <- grep("^(ob:)?z_mx$|^z_mx:ob$|^ob:z_mx$", names(b), value = TRUE)
  ib <- "ob:z_mx"; if (!ib %in% names(b)) ib <- "z_mx:ob"
  cat(sprintf("   break %6.3f (se %.3f) | bend %6.3f (se %.3f) | contrast bend-break %6.3f (se %.3f, p=%.4f)\n", b["z_mx"], sqrt(V["z_mx","z_mx"]), b["z_mx"] + b[ib],
      sqrt(V["z_mx","z_mx"] + V[ib,ib] + 2 * V["z_mx", ib]), b[ib], sqrt(V[ib, ib]), 2 * pnorm(-abs(b[ib] / sqrt(V[ib, ib])))))
}
sink("results/bend_not_break_output.txt", split = TRUE)
base <- "z_mx + z_mn + z_cap + contig + l_dur + l_same + l_diff"
cat("baseline\n"); stack(A, base)
cat("+ volume + decade FE\n"); stack(A, paste(base, "+ z_vol + decade"))
cat("contiguous only\n"); stack(A %>% filter(contig == 1), "z_mx + z_mn + z_cap + l_dur + l_same + l_diff")
cat("1946 onward\n"); stack(A %>% filter(year >= 1946), base)
sink()
