#!/usr/bin/env Rscript
# theta_sphere -- the rotation as a MERIDIAN on a sphere.
#   Fix total energy T and the nuisance. A dataset is then pinned by two angles:
#     latitude  = theta            -> the between/within split  (R = (m-2)tan^2 theta)
#     longitude = within-pattern   -> which way the within-cluster scatter points
#   => datasets of fixed T live on a sphere of radius sqrt(T). POLE = all separation
#   (within=0), EQUATOR = clusters merged (between=0). R depends ONLY on latitude.
#   Our R-fiber test fixes the longitude and sweeps theta -> ONE meridian (A,B on it).
#   (Schematic 2-sphere; the longitude sphere is higher-dim in general, but the
#    latitude=theta / R=(m-2)tan^2 theta story is exact.)
source("sims/house_style.R")
suppressMessages({ library(ggplot2) })
BLU <- "#1F5AA6"

rT <- 2.48                      # radius sqrt(T) from the 4-point toy
az <- 0.42; el <- 0.40          # viewing azimuth / elevation
P <- function(lat, lon) {       # orthographic globe projection
  X <- rT*cos(lat)*cos(lon); Y <- rT*cos(lat)*sin(lon); Z <- rT*sin(lat)
  data.frame(sx = -sin(az)*X + cos(az)*Y,
             sy = -sin(el)*cos(az)*X - sin(el)*sin(az)*Y + cos(el)*Z,
             d  =  cos(el)*cos(az)*X + cos(el)*sin(az)*Y + sin(el)*Z)
}
latc <- function(L) { lon <- seq(0, 2*pi, length.out = 160); p <- P(rep(L, 160), lon); p$g <- paste0("la", L); p }
lonc <- function(G) { lat <- seq(-pi/2, pi/2, length.out = 160); p <- P(lat, rep(G, 160)); p$g <- paste0("lo", G); p }
wire <- rbind(do.call(rbind, lapply(seq(-60, 60, 30)*pi/180, latc)),
              do.call(rbind, lapply(seq(0, 315, 45)*pi/180, lonc)))
wire$face <- ifelse(wire$d >= 0, "front", "back")

phi0 <- az + 0.62                                             # our meridian, front-facing
mer  <- P(seq(0, pi/2, length.out = 120), rep(phi0, 120))     # R-fiber: theta 0 -> pi/2
thA <- 0.60; thB <- 1.30; thO <- 1.09
ptA <- P(thA, phi0); ptB <- P(thB, phi0); ptO <- P(thO, phi0)
pole <- P(pi/2, 0); eqpt <- P(0, phi0)
equ  <- latc(0)                                               # equator highlighted
othr <- P(seq(0, pi/2, length.out = 120), rep(phi0 - 1.4, 120))  # ONE other within-pattern meridian

p <- ggplot() +
  geom_path(data = wire[wire$face == "back", ],  aes(sx, sy, group = g), colour = "grey90", linewidth = 0.3, linetype = "22") +
  geom_path(data = wire[wire$face == "front", ], aes(sx, sy, group = g), colour = "grey78", linewidth = 0.35) +
  geom_path(data = equ[equ$d >= 0, ], aes(sx, sy), colour = "grey35", linewidth = 1.0) +
  geom_path(data = equ[equ$d <  0, ], aes(sx, sy), colour = "grey65", linewidth = 0.7, linetype = "22") +
  geom_path(data = othr, aes(sx, sy), colour = "#9970AB", linewidth = 1.0) +
  # our meridian = the R-fiber
  geom_path(data = mer, aes(sx, sy), colour = BLU, linewidth = 1.9,
            arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
  geom_point(data = pole, aes(sx, sy), size = 3.4, colour = "black") +
  geom_point(data = eqpt, aes(sx, sy), size = 2.4, colour = "grey35") +
  geom_point(data = ptO, aes(sx, sy), shape = 1, size = 5.2, stroke = 1.3, colour = "black") +
  geom_point(data = ptA, aes(sx, sy), size = 3.6, colour = BLU) +
  geom_point(data = ptB, aes(sx, sy), size = 3.6, colour = BLU) +
  annotate("text", x = pole$sx, y = pole$sy + 0.30, label = "POLE:  all separation  (within = 0,  theta = 90)", size = 3.2, fontface = "bold", hjust = 0.5) +
  annotate("text", x = eqpt$sx + 0.20, y = eqpt$sy - 0.40, label = "EQUATOR:  clusters merged\n(between = 0,  theta = 0)", size = 3.0, hjust = 0.5, lineheight = 0.95, colour = "grey25") +
  annotate("text", x = ptA$sx - 0.20, y = ptA$sy + 0.02, label = "A", size = 4.2, fontface = "bold", colour = BLU, hjust = 1) +
  annotate("text", x = ptB$sx - 0.20, y = ptB$sy + 0.02, label = "B", size = 4.2, fontface = "bold", colour = BLU, hjust = 1) +
  annotate("text", x = ptO$sx + 0.18, y = ptO$sy + 0.04, label = "observed", size = 2.7, hjust = 0) +
  annotate("text", x = mer$sx[70] + 0.30, y = mer$sy[70] + 0.95, label = "our R-fiber test = ONE meridian:\nsweep latitude theta  (R = (m-2)tan^2 theta,\ndepends ONLY on latitude)", size = 3.0, colour = BLU, hjust = 0, lineheight = 0.95) +
  annotate("text", x = othr$sx[70] - 0.15, y = othr$sy[70], label = "another within-\npattern (a different\nlongitude) -- conditioned\nAWAY by our test", size = 2.7, colour = "#7B4FA0", hjust = 1, lineheight = 0.95) +
  annotate("text", x = 0, y = -rT - 0.35, label = "longitude = which within-scatter pattern (held fixed)      latitude = theta (the between/within split, what we test)", size = 2.9, colour = "grey40", hjust = 0.5) +
  coord_equal(clip = "off") +
  labs(title = "Fixed energy T -> datasets live on a sphere of radius sqrt(T);  our test walks ONE meridian (latitude = theta)") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(4, 0, 10, 0)),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(8, 95, 28, 70))
ggsave_km(p, "sims/results/theta_sphere", width = 9.2, height = 7.8)
cat("Wrote theta_sphere.\n")
