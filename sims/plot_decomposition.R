#!/usr/bin/env Rscript
# WHY two slices -- invariant map + perturbed-data film strips + truncation sets.
#   INVARIANT MAP (between vs within plane): the observed data can be moved two ways.
#     SLIDE (known sigma): a HORIZONTAL line -- within ||P1X|| fixed, vary between ||P0X||=phi.
#     ROTATE (unknown sigma): the CIRCLE  between^2 + within^2 = T  (total energy fixed).
#       "TRADE" is defined here = sliding along that circle: between up  <=>  within down.
#   Then each motion is shown on REAL perturbed data (check=clusters preserved, cross=broken)
#   and its TRUNCATION SET (the values the selective test conditions on): phi-interval for
#   the slide, R/theta-interval for the rotation.
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(ggplot2); library(gridExtra); library(grid) })
source("sims/unknown_var_union.R")
.lloyd <- KmeansInference:::.kmeans_estimation_fixed_init
.presv <- KmeansInference:::.kmeans_preserves_observed_clusters
GRN <- "#1A9850"; RED <- "#D6604D"; BLK <- km_col[["oracle"]]; BLU <- km_col[["studentized"]]

set.seed(5); npc <- 16
X <- rbind(cbind(rnorm(npc, -1.15, 0.52), rnorm(npc, 0, 0.62)),
           cbind(rnorm(npc,  1.15, 0.52), rnorm(npc, 0, 0.62)))
est <- kmeans_estimation(X, 2, 10, 2021, verbose = FALSE); cl <- est$final_cluster; init <- est$random_init_obs
P0 <- fun_P0(cl, 1, 2); P1 <- fun_P1(cl, 1, 2); P0X <- P0 %*% X; P1X <- P1 %*% X; P2X <- X - P0X - P1X
a0 <- sqrt(sum(P0X^2)); a1 <- sqrt(sum(P1X^2)); Tt <- a0^2 + a1^2; m <- sum(cl == 1) + sum(cl == 2)
U0 <- P0X / a0; U1 <- P1X / a1; th_obs <- atan2(a0, a1); rT <- sqrt(Tt)
pres <- function(Xp) { fr <- tryCatch(.lloyd(Xp, 2, init, 10, tol_eps = 1e-6, verbose = FALSE), error = function(e) NULL)
  if (is.null(fr)) FALSE else isTRUE(.presv(cl, fr$final_cluster, 1, 2)) }
Rof <- function(th) (m - 2) * tan(th)^2

phi_s <- c(0.30, 0.65, 1.00, 1.55, 2.10); phi_v <- phi_s * a0
phi_F <- lapply(phi_s, function(s) s * P0X + P1X + P2X); phi_keep <- sapply(phi_F, pres)
th_s  <- pmin(pmax(th_obs + c(-0.45, -0.22, 0, 0.24, 0.48), 0.12), 1.40); R_v <- Rof(th_s)
th_F  <- lapply(th_s, function(th) P2X + rT * (sin(th) * U0 + cos(th) * U1)); th_keep <- sapply(th_F, pres)
lim   <- max(abs(do.call(rbind, c(phi_F, th_F)))) * 1.05

frame <- function(Xp, is_obs, rail) {
  col <- if (pres(Xp)) GRN else RED
  ggplot(data.frame(x = Xp[,1], y = Xp[,2], cl = factor(cl)), aes(x, y, colour = cl)) +
    geom_point(size = 1.3) + scale_colour_manual(values = c("1" = "#4477AA", "2" = "#B2456E"), guide = "none") +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) + theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2),
          panel.border = element_rect(colour = if (is_obs) rail else col, fill = NA, linewidth = if (is_obs) 1.7 else 1.0),
          plot.background = element_rect(fill = if (is_obs) "#FFFCEF" else "white", colour = NA))
}
cap  <- function(txt, keep) textGrob(txt, gp = gpar(fontsize = 8.5, col = if (keep) GRN else RED, fontface = "bold"))
rlab <- function(txt, col) textGrob(txt, gp = gpar(fontface = "bold", fontsize = 8.5, col = col), rot = 90)
strip <- function(frames, rowlab, rcol) arrangeGrob(grobs = c(list(rlab(rowlab, rcol)),
                  mapply(function(Xp, i) frame(Xp, i == 3, rcol), frames, seq_along(frames), SIMPLIFY = FALSE)),
                  nrow = 1, widths = c(0.13, rep(1, 5)))
capstrip <- function(labs, keeps) arrangeGrob(grobs = c(list(nullGrob()), mapply(cap, labs, keeps, SIMPLIFY = FALSE)),
                  nrow = 1, widths = c(0.13, rep(1, 5)))

## truncation panels (a value-axis with the preserved interval shaded + frames placed)
trunc_panel <- function(grid, keepg, obs, ivl_fun, frame_x, frame_keep, frame_lab, xlim, xlab, ttl, prelab) {
  rr <- which(diff(c(FALSE, keepg, FALSE)) != 0); st <- rr[seq(1, length(rr), 2)]; en <- rr[seq(2, length(rr), 2)] - 1
  iv <- data.frame(lo = grid[st], hi = grid[en])
  fr <- data.frame(x = frame_x, lab = frame_lab, keep = frame_keep)
  ggplot() +
    geom_rect(data = iv, aes(xmin = lo, xmax = hi, ymin = 0, ymax = 1), fill = "#CFEBD6", colour = GRN, linewidth = 0.4) +
    annotate("text", x = mean(c(iv$lo[1], iv$hi[1])), y = 0.5, label = prelab, size = 2.8, colour = GRN) +
    geom_vline(xintercept = obs, linetype = 2, colour = "black") +
    annotate("text", x = obs, y = 1.2, label = "observed", size = 2.7, fontface = "bold") +
    geom_point(data = fr, aes(x, -0.18, colour = keep), size = 2.6) +
    geom_text(data = fr, aes(x, -0.58, label = lab), size = 2.4) +
    scale_colour_manual(values = c(`TRUE` = GRN, `FALSE` = RED), guide = "none") +
    scale_x_continuous(limits = xlim) + scale_y_continuous(limits = c(-0.8, 1.35)) +
    labs(x = xlab, y = NULL, title = ttl) +
    theme_km() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                       panel.grid = element_blank(), plot.title = element_text(size = 8.5))
}
sg <- seq(0.02, 2.6, length.out = 320); pint_phi <- trunc_panel(
  sg * a0, sapply(sg, function(s) pres(s*P0X + P1X + P2X)), a0, NULL, phi_v, phi_keep, sprintf("phi=%.1f", phi_v),
  c(0, 20), expression(paste(phi, "  =  between magnitude  ", P[0], "X")),
  "Slide truncation: phi where clustering holds (compare observed phi to chi here)",
  "clustering preserved")
tg <- seq(0.06, pi/2 - 0.04, length.out = 320); pint_R <- trunc_panel(
  tg, sapply(tg, function(th) pres(P2X + rT*(sin(th)*U0 + cos(th)*U1))), th_obs, NULL, th_s, th_keep, sprintf("R=%.0f", R_v),
  c(0, pi/2), expression(paste("angle ", theta, "    (R = (m-2)(between/within)"^2, "  = the F statistic)")),
  "Rotation truncation: the between/within RATIO R where clustering holds (compare R to F here)",
  "clustering preserved")

## ---- INVARIANT MAP: between vs within plane (line vs circle; defines TRADE) --
circ <- data.frame(th = seq(0.07, 1.5, length.out = 250)); circ$x <- rT*sin(circ$th); circ$y <- rT*cos(circ$th)
pmap <- ggplot() +
  geom_path(data = circ, aes(x, y), colour = BLU, linewidth = 1.1) +
  annotate("segment", x = 0.25*a0, xend = 2.15*a0, y = a1, yend = a1, colour = BLK, linewidth = 1.1,
           arrow = arrow(ends = "both", length = unit(0.14, "cm"))) +
  geom_point(data = data.frame(x = phi_v, y = a1, k = phi_keep), aes(x, y, colour = k), size = 2.5) +
  geom_point(data = data.frame(x = rT*sin(th_s), y = rT*cos(th_s), k = th_keep), aes(x, y, colour = k), size = 2.5) +
  annotate("point", x = a0, y = a1, size = 3.4, colour = "black") +
  annotate("text", x = a0 + 0.25, y = a1 + 0.45, label = "observed", size = 2.8, fontface = "bold", hjust = 0) +
  annotate("text", x = 9.3, y = a1 - 0.75, label = "SLIDE (known sigma):\nvary the DISTANCE only,\nwithin-spread fixed", size = 2.9, colour = BLK, hjust = 0, lineheight = 0.92) +
  annotate("text", x = 7.6, y = 8.5, label = "ROTATE (unknown sigma):\nshift the between/within RATIO,\ntotal  b^2 + w^2 = T  fixed", size = 2.9, colour = BLU, hjust = 0, lineheight = 0.92) +
  scale_colour_manual(values = c(`TRUE` = GRN, `FALSE` = RED), guide = "none") +
  coord_equal(xlim = c(0, 16), ylim = c(0, 9.6), expand = FALSE) +
  labs(x = expression(paste("between   ", P[0], "X")), y = expression(paste("within   ", P[1], "X")),
       title = "Invariant: slide = line (within fixed),  rotate = circle (T fixed)") +
  theme_km() + theme(panel.grid.major = element_blank(), plot.title = element_text(size = 10))

g <- arrangeGrob(
  pmap,
  strip(phi_F, "KNOWN sigma\nSLIDE: vary separation", BLK), capstrip(sprintf("phi=%.1f %s", phi_v, ifelse(phi_keep,"✓","✗")), phi_keep), pint_phi,
  strip(th_F,  "UNKNOWN sigma\nROTATE: vary ratio", BLU), capstrip(sprintf("R=%.0f %s", R_v, ifelse(th_keep,"✓","✗")), th_keep), pint_R,
  ncol = 1, heights = c(2.3, 1, 0.12, 0.6, 1, 0.12, 0.6),
  top = textGrob("Perturb the data, keep the clustering:  SLIDE (within fixed) vs ROTATE (total T fixed)",
                 gp = gpar(fontface = "bold", fontsize = 12)))
ggsave_km(g, "sims/results/decomposition", width = 12, height = 11.5)
cat(sprintf("Wrote decomposition. T=%.0f (invariant)  phi-trunc & R-trunc panels added; trade defined on the map.\n", Tt))
