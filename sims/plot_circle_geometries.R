#!/usr/bin/env Rscript
# EXPERIMENT: hang real data geometries on the two motions, side by side.
#   LEFT  ROTATE (unknown sigma): geometries on the CIRCLE between^2+within^2=T.
#         Round the arc the clusters MORPH: overlapping blob (high within) -> two
#         tight far-apart groups (high between). Only the between/within ratio moves.
#   RIGHT SLIDE (known sigma): geometries on a horizontal LINE (within fixed). The
#         cluster shape stays the same; the clusters just slide apart.
# Green frame = k-means still returns the two clusters; red = it breaks; ring = observed.
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(ggplot2); library(gridExtra); library(grid) })
source("sims/unknown_var_union.R")
.lloyd <- KmeansInference:::.kmeans_estimation_fixed_init
.presv <- KmeansInference:::.kmeans_preserves_observed_clusters
GRN <- "#1A9850"; RED <- "#D6604D"; BLU <- km_col[["studentized"]]; BLK <- km_col[["oracle"]]

set.seed(5); npc <- 16
X <- rbind(cbind(rnorm(npc, -1.15, 0.52), rnorm(npc, 0, 0.62)),
           cbind(rnorm(npc,  1.15, 0.52), rnorm(npc, 0, 0.62)))
est <- kmeans_estimation(X, 2, 10, 2021, verbose = FALSE); cl <- est$final_cluster; init <- est$random_init_obs
P0 <- fun_P0(cl, 1, 2); P1 <- fun_P1(cl, 1, 2); P0X <- P0 %*% X; P1X <- P1 %*% X; P2X <- X - P0X - P1X
a0 <- sqrt(sum(P0X^2)); a1 <- sqrt(sum(P1X^2)); Tt <- a0^2 + a1^2; rT <- sqrt(Tt)
U0 <- P0X / a0; U1 <- P1X / a1; th_obs <- atan2(a0, a1)
pres <- function(Xp) { fr <- tryCatch(.lloyd(Xp, 2, init, 10, tol_eps = 1e-6, verbose = FALSE), error = function(e) NULL)
  if (is.null(fr)) FALSE else isTRUE(.presv(cl, fr$final_cluster, 1, 2)) }

# big thumbnail at an anchor + a thin connector to the anchor's point on the curve
necklace <- function(params, curvept, thumbpt, perturb, SC) {
  pts <- list(); box <- list(); seg <- list()
  for (k in seq_along(params)) {
    Xk <- perturb(params[k]); Xc <- scale(Xk, center = TRUE, scale = FALSE); kp <- pres(Xk)
    cpt <- curvept(params[k]); tpt <- thumbpt(params[k])
    px <- tpt[1] + SC*Xc[,1]; py <- tpt[2] + SC*Xc[,2]; pad <- 0.35
    pts[[k]] <- data.frame(x = px, y = py, cl = factor(cl), k = k)
    box[[k]] <- data.frame(xmin = min(px)-pad, xmax = max(px)+pad, ymin = min(py)-pad, ymax = max(py)+pad, keep = kp)
    seg[[k]] <- data.frame(x = cpt[1], y = cpt[2], xend = tpt[1], yend = tpt[2], keep = kp)
  }
  list(pts = do.call(rbind, pts), box = do.call(rbind, box), seg = do.call(rbind, seg))
}
gpanel <- function(nk, curve, anchor_obs, ttl, lab, lab_xy, xl, yl) {
  ggplot() +
    geom_path(data = curve, aes(x, y), colour = curve$col[1], linewidth = 1.2) +
    geom_segment(data = nk$seg, aes(x = x, y = y, xend = xend, yend = yend, colour = keep), linewidth = 0.4) +
    geom_point(data = nk$seg, aes(x, y, colour = keep), size = 1.6) +
    geom_rect(data = nk$box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, colour = keep), fill = "white", linewidth = 0.8) +
    geom_point(data = nk$pts, aes(x, y, fill = cl), shape = 21, colour = "white", stroke = 0.15, size = 2.6) +
    annotate("point", x = anchor_obs[1], y = anchor_obs[2], shape = 1, size = 8, stroke = 1.3, colour = "black") +
    annotate("text", x = lab_xy[1], y = lab_xy[2], label = lab, size = 3.2, colour = curve$col[1], hjust = 0, lineheight = 0.92) +
    scale_fill_manual(values = c("1" = "#4477AA", "2" = "#B2456E"), guide = "none") +
    scale_colour_manual(values = c(`TRUE` = GRN, `FALSE` = RED), guide = "none") +
    coord_equal(xlim = xl, ylim = yl, expand = FALSE) +
    labs(x = expression(paste("between   ", P[0], "X")), y = expression(paste("within   ", P[1], "X")), title = ttl) +
    theme_km() + theme(panel.grid.major = element_blank(), plot.title = element_text(size = 11))
}

## LEFT: ROTATE -- thumbnails pushed OUTSIDE the circle, connected
thv <- seq(0.30, 1.40, length.out = 5); OFF <- 3.6; SCc <- 0.62
nkC <- necklace(thv, function(th) c(rT*sin(th), rT*cos(th)),
                function(th) c((rT+OFF)*sin(th), (rT+OFF)*cos(th)),
                function(th) P2X + rT*(sin(th)*U0 + cos(th)*U1), SCc)
circ <- data.frame(a = seq(0.05, 1.55, length.out = 200)); circ$x <- rT*sin(circ$a); circ$y <- rT*cos(circ$a); circ$col <- BLU
pL <- gpanel(nkC, circ, c(rT*sin(th_obs), rT*cos(th_obs)),
             "ROTATE (unknown sigma): fix total, change the ratio", "shape MORPHS\n(spread -> tight)", c(0.3, 13.6),
             xl = c(-0.5, 14.5), yl = c(-0.5, 14.5))
# the INVARIANT: every frame is at the same radius sqrt(T) from the origin
spk <- data.frame(cx = rT*sin(thv), cy = rT*cos(thv))
pL <- pL +
  geom_segment(data = spk, aes(x = 0, y = 0, xend = cx, yend = cy), colour = km_ref, linetype = "22", linewidth = 0.35) +
  annotate("text", x = 0.3, y = 2.6, label = "INVARIANT: every frame is the\nsame radius sqrt(T) from the origin\n(between² + within² = T)",
           size = 2.95, colour = km_ref, hjust = 0, lineheight = 0.92)

## RIGHT: SLIDE -- thumbnails ABOVE the line, connected
phv <- seq(3.5, 13.0, length.out = 5); OFF2 <- 3.3; SCl <- 0.62
nkL <- necklace(phv, function(p) c(p, a1), function(p) c(p, a1 + OFF2),
                function(p) (p/a0)*P0X + P1X + P2X, SCl)
line <- data.frame(x = c(2, 14.5), y = a1, col = BLK)
pR <- gpanel(nkL, line, c(a0, a1),
             "SLIDE (known sigma): fix within, change the distance", "same SHAPE,\njust slides apart", c(2.2, 10.2),
             xl = c(1, 15.5), yl = c(-3.5, 11.5))

g <- arrangeGrob(pL, pR, ncol = 2,
                 top = textGrob("Two motions on real data:  rotate = same energy, different ratio   vs   slide = same shape, different distance",
                                gp = gpar(fontface = "bold", fontsize = 12)))
ggsave_km(g, "sims/results/circle_geometries", width = 12.5, height = 6.4)
cat(sprintf("Wrote circle_geometries (2-panel, bigger thumbnails). T=%.0f.\n", Tt))
