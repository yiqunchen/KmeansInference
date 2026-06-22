#!/usr/bin/env Rscript
# theta_toy -- the SMALLEST possible example: 4 points (2 per cluster).
#   (1 point per cluster has zero within-scatter -> nothing to rotate; need >=2.)
#   Same three views as theta_matching but you can watch every single point:
#   A feature space (each of the 4 points rotates to its own image, row i->row i),
#   B polar (both datasets on the same circle r=sqrt(T)), C invariance (b^2+w^2=T).
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(ggplot2); library(gridExtra); library(grid) })
source("sims/unknown_var_union.R")
CL1 <- "#4477AA"; CL2 <- "#B2456E"; BETW <- "#D55E00"; WITH <- "#0072B2"

## four hand-placed points: 2 left (cluster 1), 2 right (cluster 2)
X  <- rbind(c(-1.2, 0.7), c(-0.8, -0.5), c(1.1, 0.6), c(1.3, -0.4)); cl <- c(1, 1, 2, 2)
P0 <- fun_P0(cl, 1, 2); P1 <- fun_P1(cl, 1, 2); P0X <- P0 %*% X; P1X <- P1 %*% X; P2X <- X - P0X - P1X
a0 <- sqrt(sum(P0X^2)); a1 <- sqrt(sum(P1X^2)); Tt <- a0^2 + a1^2; rT <- sqrt(Tt)
U0 <- P0X / a0; U1 <- P1X / a1; th_obs <- atan2(a0, a1)
Xof <- function(th) P2X + rT*(sin(th)*U0 + cos(th)*U1)
thA <- 0.60; thB <- 1.30

####################### A: feature space, 4 points ##########################
tg   <- seq(thA, thB, length.out = 60)
arc1 <- function(i) { xy <- t(sapply(tg, function(t) Xof(t)[i, ]))
  data.frame(x = xy[,1], y = xy[,2], i = i, cl = factor(cl[i])) }
arcs <- do.call(rbind, lapply(1:4, arc1))
A <- data.frame(x = Xof(thA)[,1], y = Xof(thA)[,2], cl = factor(cl), lab = 1:4)
B <- data.frame(x = Xof(thB)[,1], y = Xof(thB)[,2], cl = factor(cl), lab = 1:4)

pA <- ggplot() +
  geom_path(data = arcs, aes(x, y, group = i, colour = cl), linewidth = 1.1,
            arrow = arrow(length = unit(0.18, "cm"), type = "closed"), alpha = 0.9) +
  geom_point(data = A, aes(x, y, colour = cl), shape = 1,  size = 5.5, stroke = 1.6) +
  geom_point(data = B, aes(x, y, fill = cl), shape = 21, size = 5.0, stroke = 0.3, colour = "white") +
  geom_text(data = A, aes(x, y, label = lab), size = 3.0, fontface = "bold", colour = "grey20", nudge_x = -0.17, nudge_y = 0.17) +
  geom_text(data = B, aes(x, y, label = lab), size = 3.0, fontface = "bold", colour = "white") +
  scale_colour_manual(values = c("1" = CL1, "2" = CL2), guide = "none") +
  scale_fill_manual(values   = c("1" = CL1, "2" = CL2), guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.4, size = 3.0, colour = "grey20", lineheight = 1.0,
           label = "open = A (theta=0.60)\nfilled = B (theta=1.30)\narrow = each point rotates to itself") +
  coord_equal(clip = "off") +
  labs(x = expression(x[1]), y = expression(x[2]), title = "A.  4 points, matched 1-to-1") +
  theme_km() + theme(plot.title = element_text(size = 10.5), panel.grid.minor = element_blank())

####################### B: polar (r, theta) ##################################
gA   <- seq(0, pi/2, length.out = 140)
ring <- function(r) data.frame(x = r*sin(gA), y = r*cos(gA), r = r)
guides <- do.call(rbind, lapply(seq(0.5, 2.5, by = 0.5), ring)); main <- ring(rT)
radl <- do.call(rbind, lapply(seq(0, pi/2, by = pi/12), function(a) data.frame(x = c(0, (rT+0.4)*sin(a)), y = c(0, (rT+0.4)*cos(a)), a = a)))
ptsAB <- data.frame(x = c(rT*sin(thA), rT*sin(thB)), y = c(rT*cos(thA), rT*cos(thB)))
pB <- ggplot() +
  geom_line(data = radl, aes(x, y, group = a), colour = "grey90", linewidth = 0.3) +
  geom_path(data = guides, aes(x, y, group = r), colour = "grey90", linewidth = 0.3) +
  geom_path(data = main, aes(x, y), colour = "grey45", linewidth = 1.2) +
  geom_segment(aes(x = 0, y = 0, xend = rT*sin(thA), yend = rT*cos(thA)), colour = "grey25", linewidth = 0.6) +
  geom_segment(aes(x = 0, y = 0, xend = rT*sin(thB), yend = rT*cos(thB)), colour = "grey25", linewidth = 0.6) +
  geom_point(data = ptsAB[1,], aes(x, y), shape = 1,  size = 4.6, stroke = 1.5, colour = "black") +
  geom_point(data = ptsAB[2,], aes(x, y), shape = 16, size = 4.0, colour = "black") +
  annotate("text", x = rT*sin(thA)-0.12, y = rT*cos(thA)+0.18, label = "A", size = 4, fontface = "bold", hjust = 1) +
  annotate("text", x = rT*sin(thB)+0.16, y = rT*cos(thB)+0.06, label = "B", size = 4, fontface = "bold", hjust = 0) +
  annotate("text", x = (rT+0.28)*sin(0.93), y = (rT+0.28)*cos(0.93), label = sprintf("r = sqrt(T) = %.2f", rT), angle = -38, size = 3.0, colour = "grey30") +
  coord_equal(xlim = c(0, rT+0.6), ylim = c(0, rT+0.6), expand = FALSE) +
  labs(x = "between (radial)", y = "within (radial)", title = "B.  polar: same radius, angle moves") +
  theme_km() + theme(plot.title = element_text(size = 10.5), panel.grid = element_blank())

####################### C: the invariance ####################################
ener <- data.frame(d = rep(c("A", "B"), each = 2),
                   comp = factor(rep(c("between^2", "within^2"), 2), levels = c("within^2", "between^2")),
                   val  = c(Tt*sin(thA)^2, Tt*cos(thA)^2, Tt*sin(thB)^2, Tt*cos(thB)^2))
pC <- ggplot(ener, aes(d, val, fill = comp)) +
  geom_col(width = 0.62, colour = "white") +
  geom_text(aes(label = sprintf("%.1f", val)), position = position_stack(vjust = 0.5), size = 3.2, colour = "white", fontface = "bold") +
  geom_hline(yintercept = Tt, linetype = 2, colour = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = Tt + 0.45, label = sprintf("total = T = %.1f  (invariant)", Tt), fontface = "bold", size = 3.1) +
  scale_fill_manual(values = c("between^2" = BETW, "within^2" = WITH),
                    labels = c("between^2" = expression(between^2), "within^2" = expression(within^2)), name = NULL) +
  scale_x_discrete(labels = c("A" = "A\n(theta=0.60)", "B" = "B\n(theta=1.30)")) +
  coord_cartesian(ylim = c(0, Tt*1.18)) +
  labs(x = NULL, y = "energy", title = "C.  invariance: split moves, total = T") +
  theme_km() + theme(plot.title = element_text(size = 10.5), legend.position = "top",
                     legend.margin = margin(0, 0, -4, 0), panel.grid.major.x = element_blank())

g <- arrangeGrob(pA, pB, pC, ncol = 3, widths = c(1.3, 1.0, 0.92),
                 top = textGrob("Smallest example -- 4 points (2 per cluster): one angle theta rotates each point to itself; between^2+within^2 = T stays fixed",
                                gp = gpar(fontface = "bold", fontsize = 12)))
ggsave_km(g, "sims/results/theta_toy", width = 15.5, height = 5.6)
cat(sprintf("Wrote theta_toy. a0=%.2f a1=%.2f T=%.2f th_obs=%.2f | A:%.1f+%.1f  B:%.1f+%.1f\n",
            a0, a1, Tt, th_obs, Tt*sin(thA)^2, Tt*cos(thA)^2, Tt*sin(thB)^2, Tt*cos(thB)^2))
