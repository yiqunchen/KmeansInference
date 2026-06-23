# ---------------------------------------------------------------------------
# plot_conditioning_set.R -- the signature "condition on less" figure, following
# Chen-Witten (k-means) Fig 3 + Chen-Jewell-Witten (GFL) Fig 2(d) / Fig 1(d).
#
# Top: the data perturbed along the contrast, x'(phi), at three phi (a,b,c). The
#   two tested clusters CHANGE outside the conditioning set (a), are PRESERVED at a
#   smaller separation the union reaches but the path does not (b), and at the
#   observed value (c). Colours are aligned across panels so a changed partition
#   visibly recolours.
# Bottom: the set of phi preserving the two clusters. The PATH test conditions on
#   the single interval at the observed statistic; the UNION test on the whole set,
#   which reaches BELOW the observed statistic (dashed) -- exactly why
#   p_union << p_path. Truncated-chi null underneath.
# ---------------------------------------------------------------------------
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(intervals); library(ggplot2)
                   library(dplyr); library(gridExtra); library(grid) })

SEED <- 61; D <- 2.7; NPER <- 12; Q <- 2; K <- 3; C1 <- 1; C2 <- 2
set.seed(SEED)
mu <- rbind(c(0,0), c(D,0), c(D/2, D*sqrt(3)/2))
X  <- do.call(rbind, lapply(1:3, function(g) matrix(rnorm(NPER*Q),NPER,Q) + matrix(mu[g,],NPER,Q,byrow=TRUE)))
fit <- kmeans_inference_union(X, K, C1, C2, sig = 1, seed = SEED)
cl0 <- fit$final_cluster; n1 <- sum(cl0==C1); n2 <- sum(cl0==C2)

nu <- ifelse(cl0==C1, 1/n1, ifelse(cl0==C2, -1/n2, 0))
XTnu <- as.numeric(t(X) %*% nu); phi_obs <- sqrt(sum(XTnu^2)); dir <- XTnu/phi_obs; nu_n2 <- sum(nu^2)
xprime <- function(phi) X + ((phi - phi_obs)/nu_n2) * outer(nu, dir)
recluster <- function(phi) kmeans_estimation(xprime(phi), K, seed = SEED)$final_cluster
agree <- function(cl) sum(cl == cl0) / length(cl0)
align <- function(cl) { tab <- table(cl0, cl); map <- rownames(tab)[apply(tab, 2, which.max)]
                        as.integer(map[match(cl, as.integer(colnames(tab)))]) }

iu <- as.matrix(fit$interval_union); ip <- as.matrix(fit$interval_path)
ius <- iu[order(iu[,1]), , drop = FALSE]; union_lo <- min(ius[,1]); path_lo <- ip[1,1]
sf <- fit$scale_factor; ts <- fit$test_stat
cat(sprintf("seed %d: stat=%.3f p_path=%.4f p_union=%.4f | union %s\n", SEED, ts,
            fit$pval_path, fit$pval_union, paste(round(range(iu),2), collapse="-")))

# phi (a) changed: scan below the union for the most-changed partition
cand <- seq(max(0.4, union_lo - 1.8), union_lo - 0.2, length.out = 24)
phi_a <- cand[which.min(sapply(cand, function(p) agree(recluster(p))))]
phi_b <- mean(c(union_lo, path_lo))           # (b) preserved, in union but below path
phi_c <- phi_obs                               # (c) observed
phis  <- c(phi_a, phi_b, phi_c); tags <- letters[1:3]

# ---- top row: three perturbed-data scatters ------------------------------
clr <- c("1" = "#CC6699", "2" = "#1F5AA6", "3" = "#E69F00")   # pink / blue / orange
scatter <- function(phi, tg) {
  cl <- recluster(phi); a <- agree(cl); Xp <- xprime(phi)
  df <- data.frame(x = Xp[,1], y = Xp[,2], cl = factor(align(cl)))
  changed <- a < 0.9
  pre <- sprintf("(%s)  ", tg)
  ttl <- if (abs(phi - phi_obs) < 1e-6) bquote(bold(.(pre)) * phi[obs] == .(sprintf("%.2f", phi))) else
                                        bquote(bold(.(pre)) * phi == .(sprintf("%.2f", phi)))
  ggplot(df, aes(x, y, colour = cl)) +
    geom_point(size = 1.8, alpha = 0.9) +
    scale_colour_manual(values = clr, guide = "none") +
    labs(title = ttl, subtitle = if (changed) "clusters 1,2 change" else "clusters 1,2 preserved",
         x = NULL, y = NULL) +
    theme_km(12) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          plot.subtitle = element_text(face = "bold", size = rel(0.85),
                          colour = if (changed) km_col[["sample"]] else km_col[["studentized"]]))
}
tops <- Map(scatter, phis, tags)

# ---- bottom: the conditioning set + truncated-chi null --------------------
xlo <- min(phi_a, union_lo) - 0.2; xhi <- max(ts, max(iu)) + 0.25
gx <- seq(xlo, xhi, length.out = 700); dens <- function(x){ z <- x^2/sf; dchisq(z, df=Q)*(2*x/sf) }
dd <- data.frame(x = gx, d = dens(gx))
seg_u <- as.data.frame(iu); seg_p <- as.data.frame(ip); names(seg_u) <- names(seg_p) <- c("lo","hi")
ymax <- max(dd$d); yu <- ymax*1.20; yp <- ymax*1.05
mk <- data.frame(x = phis, tg = paste0("(", tags, ")"))
bottom <- ggplot() +
  geom_area(data = dd, aes(x, d), fill = "grey88") +
  geom_line(data = dd, aes(x, d), colour = "grey55", linewidth = 0.4) +
  geom_segment(data = seg_u, aes(x=lo, xend=hi, y=yu, yend=yu), colour = km_col[["studentized"]], linewidth = 5, lineend = "butt") +
  geom_segment(data = seg_p, aes(x=lo, xend=hi, y=yp, yend=yp), colour = "black", linewidth = 5, lineend = "butt") +
  geom_vline(xintercept = ts, linetype = "dashed", linewidth = 0.6, colour = "black") +
  annotate("text", x = ts, y = yu*1.07, label = "observed statistic", hjust = -0.03, size = 3.1, fontface = "bold") +
  annotate("text", x = min(iu), y = yu, label = "Union set  ", hjust = 1.03, size = 3.2, fontface = "bold", colour = km_col[["studentized"]]) +
  annotate("text", x = ip[1,1], y = yp, label = "Path set  ", hjust = 1.03, size = 3.2, fontface = "bold") +
  geom_point(data = mk, aes(x = x, y = 0), shape = 17, size = 2.6) +
  geom_text(data = mk, aes(x = x, y = ymax*0.13, label = tg), size = 3.1, fontface = "bold") +
  labs(x = expression(phi == "perturbation along the contrast " * (group("||", X^T * nu, "||"))), y = "truncated null") +
  scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(xlo, xhi)) +
  theme_km(12) + theme(axis.title.y = element_text(size = rel(0.85)))

g <- arrangeGrob(arrangeGrob(grobs = tops, ncol = 3), bottom, ncol = 1, heights = c(1.05, 1))
ggsave_km(g, "sims/results/conditioning_set", width = 11, height = 6.6)
cat("Wrote conditioning_set.{pdf,png}\n")
