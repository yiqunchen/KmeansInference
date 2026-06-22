#!/usr/bin/env Rscript
# theta_matching -- two datasets on the spectrum, three views of the SAME rotation.
#   A  FEATURE SPACE : overlay x(thetaA) (open) and x(thetaB) (filled); each point i
#                      is joined to its own image by the arc it rides (row i -> row i).
#   B  POLAR (r,theta): the (between,within) summary in polar form -- both datasets sit
#                      on the SAME circle r = sqrt(T); only the angle theta differs.
#   C  INVARIANCE    : between^2 + within^2 = T for BOTH -- the split moves, total stays.
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(ggplot2); library(gridExtra); library(grid) })
source("sims/unknown_var_union.R")
CL1 <- "#4477AA"; CL2 <- "#B2456E"; BETW <- "#D55E00"; WITH <- "#0072B2"

set.seed(5); npc <- 16
X <- rbind(cbind(rnorm(npc, -1.15, 0.52), rnorm(npc, 0, 0.62)),
           cbind(rnorm(npc,  1.15, 0.52), rnorm(npc, 0, 0.62)))
est <- kmeans_estimation(X, 2, 10, 2021, verbose = FALSE); cl <- est$final_cluster
P0 <- fun_P0(cl, 1, 2); P1 <- fun_P1(cl, 1, 2); P0X <- P0 %*% X; P1X <- P1 %*% X; P2X <- X - P0X - P1X
a0 <- sqrt(sum(P0X^2)); a1 <- sqrt(sum(P1X^2)); Tt <- a0^2 + a1^2; rT <- sqrt(Tt)
U0 <- P0X / a0; U1 <- P1X / a1
Xof <- function(th) P2X + rT*(sin(th)*U0 + cos(th)*U1)        # the dataset at angle theta
thA <- 0.62; thB <- 1.24

####################### A: feature space, 1-to-1 match #######################
tg   <- seq(thA, thB, length.out = 50)
arc1 <- function(i) { xy <- t(sapply(tg, function(t) Xof(t)[i, ]))
  data.frame(x = xy[,1], y = xy[,2], i = i, cl = factor(cl[i])) }
ell <- rT * sqrt(rowSums(U0^2) + rowSums(U1^2)); i1 <- which(cl == 1); i2 <- which(cl == 2)
hl  <- c(i1[order(ell[i1], decreasing = TRUE)][1:2], i2[order(ell[i2], decreasing = TRUE)][1:2])
allarc <- do.call(rbind, lapply(setdiff(seq_len(nrow(X)), hl), arc1))
hlarc  <- do.call(rbind, lapply(hl, arc1))
A   <- data.frame(x = Xof(thA)[,1],  y = Xof(thA)[,2],  cl = factor(cl))
B   <- data.frame(x = Xof(thB)[,1],  y = Xof(thB)[,2],  cl = factor(cl))
hlS <- data.frame(x = Xof(thA)[hl,1], y = Xof(thA)[hl,2], lab = seq_along(hl))

pA <- ggplot() +
  geom_path(data = allarc, aes(x, y, group = i, colour = cl), linewidth = 0.35, alpha = 0.28) +
  geom_point(data = A, aes(x, y, colour = cl), shape = 1,  size = 1.9, stroke = 0.7, alpha = 0.55) +
  geom_point(data = B, aes(x, y, fill   = cl), shape = 21, size = 2.4, stroke = 0.2, colour = "white", alpha = 0.9) +
  geom_path(data = hlarc, aes(x, y, group = i, colour = cl), linewidth = 1.0,
            arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
  geom_point(data = hlS, aes(x, y), shape = 1, size = 3.6, stroke = 1.1, colour = "grey15") +
  geom_text(data = hlS, aes(x, y, label = lab), nudge_x = -0.22, nudge_y = 0.20, size = 3.3, fontface = "bold", colour = "grey15") +
  scale_colour_manual(values = c("1" = CL1, "2" = CL2), guide = "none") +
  scale_fill_manual(values   = c("1" = CL1, "2" = CL2), guide = "none") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 2.9, colour = "grey15", lineheight = 1.0,
           label = "open = A (theta=0.62)\nfilled = B (theta=1.24)\narrow = point i rotates to itself") +
  coord_equal(clip = "off") +
  labs(x = expression(x[1]), y = expression(x[2]), title = "A.  feature space: matched 1-to-1") +
  theme_km() + theme(plot.title = element_text(size = 10.5), panel.grid.minor = element_blank())

####################### B: polar (r, theta) ##################################
gA   <- seq(0, pi/2, length.out = 140)
ring <- function(r) data.frame(x = r*sin(gA), y = r*cos(gA), r = r)
guides <- do.call(rbind, lapply(c(2, 4, 6, 8), ring)); main <- ring(rT)
radl <- do.call(rbind, lapply(seq(0, pi/2, by = pi/12), function(a) data.frame(x = c(0, 9.2*sin(a)), y = c(0, 9.2*cos(a)), a = a)))
ptsAB <- data.frame(x = c(rT*sin(thA), rT*sin(thB)), y = c(rT*cos(thA), rT*cos(thB)))
pB <- ggplot() +
  geom_line(data = radl, aes(x, y, group = a), colour = "grey90", linewidth = 0.3) +
  geom_path(data = guides, aes(x, y, group = r), colour = "grey90", linewidth = 0.3) +
  geom_path(data = main, aes(x, y), colour = "grey45", linewidth = 1.2) +
  geom_segment(aes(x = 0, y = 0, xend = rT*sin(thA), yend = rT*cos(thA)), colour = "grey25", linewidth = 0.6) +
  geom_segment(aes(x = 0, y = 0, xend = rT*sin(thB), yend = rT*cos(thB)), colour = "grey25", linewidth = 0.6) +
  geom_point(data = ptsAB[1,], aes(x, y), shape = 1,  size = 4.6, stroke = 1.5, colour = "black") +
  geom_point(data = ptsAB[2,], aes(x, y), shape = 16, size = 4.0, colour = "black") +
  annotate("text", x = rT*sin(thA)-0.3, y = rT*cos(thA)+0.5, label = "A", size = 4, fontface = "bold", hjust = 1) +
  annotate("text", x = rT*sin(thB)+0.5, y = rT*cos(thB)+0.2, label = "B", size = 4, fontface = "bold", hjust = 0) +
  annotate("text", x = (rT+0.7)*sin(0.93), y = (rT+0.7)*cos(0.93), label = "r = sqrt(T) = 8.26", angle = -38, size = 3.0, colour = "grey30") +
  annotate("text", x = 1.5*sin(thA), y = 1.5*cos(thA)+0.1, label = "theta[A]", parse = TRUE, size = 3.1, colour = "grey25") +
  annotate("text", x = 2.0*sin(thB), y = 2.0*cos(thB), label = "theta[B]", parse = TRUE, size = 3.1, colour = "grey25") +
  coord_equal(xlim = c(0, 9.6), ylim = c(0, 9.6), expand = FALSE) +
  labs(x = "between  (radial)", y = "within  (radial)",
       title = "B.  polar (r, theta): same radius, angle moves") +
  theme_km() + theme(plot.title = element_text(size = 10.5), panel.grid = element_blank())

####################### C: the invariance ####################################
ener <- data.frame(d = rep(c("A", "B"), each = 2),
                   comp = factor(rep(c("between^2", "within^2"), 2), levels = c("within^2", "between^2")),
                   val  = c(Tt*sin(thA)^2, Tt*cos(thA)^2, Tt*sin(thB)^2, Tt*cos(thB)^2))
pC <- ggplot(ener, aes(d, val, fill = comp)) +
  geom_col(width = 0.62, colour = "white") +
  geom_text(aes(label = sprintf("%.0f", val)), position = position_stack(vjust = 0.5), size = 3.2, colour = "white", fontface = "bold") +
  geom_hline(yintercept = Tt, linetype = 2, colour = "black", linewidth = 0.5) +
  annotate("text", x = 1.5, y = Tt + 5, label = "total = T = 68  (invariant)", fontface = "bold", size = 3.1) +
  scale_fill_manual(values = c("between^2" = BETW, "within^2" = WITH),
                    labels = c("between^2" = expression(between^2), "within^2" = expression(within^2)), name = NULL) +
  scale_x_discrete(labels = c("A" = "A\n(theta=0.62)", "B" = "B\n(theta=1.24)")) +
  coord_cartesian(ylim = c(0, Tt*1.18)) +
  labs(x = NULL, y = "energy", title = "C.  invariance: split moves, total = T") +
  theme_km() + theme(plot.title = element_text(size = 10.5), legend.position = "top",
                     legend.margin = margin(0, 0, -4, 0), panel.grid.major.x = element_blank())

g <- arrangeGrob(pA, pB, pC, ncol = 3, widths = c(1.35, 1.0, 0.92),
                 top = textGrob("Two datasets, one rotation by theta:  matched 1-to-1  ->  same circle r=sqrt(T)  ->  between^2+within^2=T fixed",
                                gp = gpar(fontface = "bold", fontsize = 12)))
ggsave_km(g, "sims/results/theta_matching", width = 15.5, height = 5.6)
cat(sprintf("Wrote theta_matching (3 panels). A:%.0f+%.0f  B:%.0f+%.0f  =T=%.0f\n",
            Tt*sin(thA)^2, Tt*cos(thA)^2, Tt*sin(thB)^2, Tt*cos(thB)^2, Tt))
