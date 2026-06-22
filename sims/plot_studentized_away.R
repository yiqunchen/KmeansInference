#!/usr/bin/env Rscript
# "The union gain is studentized away" -- the SIDEDNESS mechanism.
#
# Each test reports a one-sided truncated p-value
#     p = P(stat >= observed | stat in S),
# and the UNION enlarges S.  An enlargement lowers a one-sided p-value ONLY if the
# extra set-mass lands BELOW the observed statistic (it then inflates the
# denominator alone).  Both observed statistics here sit deep in their MARGINAL
# tail -- so this is NOT a "bulk vs tail" story.  The difference is the SIDE the
# union's extra room falls on:
#
#   Known variance (chi): other Lloyd paths tolerate SMALLER separation, so the
#       union extends the truncation interval DOWNWARD (below observed) -- pure
#       denominator -> p_union << p_path.
#   Unknown var (F): k-means MINIMIZED the within-cluster scatter, pinning the
#       studentized R against its low-R merge boundary; the only extra room is at
#       LARGER R (above observed) -> added to numerator AND denominator -> no gain.
#
# We visualize the CONDITIONAL reference density restricted to (and normalized
# over) each union set -- computed in log-space so the ~1e-60 absolute masses do
# not underflow.  Area below the observed line = the denominator the union adds.
#
# Output: sims/results/studentized_away.{pdf,png}
source("sims/house_style.R")
suppressMessages({ library(KmeansInference); library(ggplot2)
                   library(gridExtra); library(grid) })
source("sims/kmeans_union_unknownvar.R")

set_seed_data <- function(seed, q = 10, n_per = 30, delta = 4) {
  set.seed(seed)
  mu <- rbind(c(0, 0), c(delta, 0), c(delta/2, delta*sqrt(3)/2))
  mu <- cbind(mu, matrix(0, 3, q - 2))
  tc <- rep(1:3, each = n_per)
  matrix(rnorm(3*n_per*q), 3*n_per, q) + mu[tc, ]
}
in_set <- function(x, M) Reduce(`|`, lapply(seq_len(nrow(M)),
                       function(i) x >= M[i,1] & x <= M[i,2]), rep(FALSE, length(x)))
# conditional reference density on the union set, normalized to integrate to 1,
# evaluated on a grid -- log-space so 1e-60 masses don't underflow. logf = log ref.
cond_density <- function(grid, logf, M) {
  keep <- in_set(grid, M); f <- rep(0, length(grid))
  lw <- logf[keep]; lw <- lw - max(lw); w <- exp(lw)           # relative weights O(1)
  dx <- median(diff(grid)); f[keep] <- w / (sum(w) * dx)        # proper density
  f
}

## representative rep: clear known-sigma gain, clean R-fiber solve (message is universal)
q <- 10; pick <- NULL
for (s in c(3, 1, 6, 8, 10, 11, 4, 9)) {
  X  <- set_seed_data(s, q = q)
  ft <- tryCatch(kmeans_inference_union(X, 3, 1, 3, sig = 1, seed = 2021), error = function(e) NULL)
  if (is.null(ft) || !is.finite(ft$pval_union) || !is.finite(ft$pval_path)) next
  if (ft$pval_path < 0.6 * ft$pval_union) next                  # want a real chi gain
  rf <- tryCatch(kmeans_union_unknownvar(X, 3, 1, 3, seed = 2021, n_theta = 400), error = function(e) NULL)
  if (!is.null(rf) && is.finite(rf$p_union) && is.finite(rf$p_path)) { pick <- list(s=s, X=X, fit=ft, rf=rf); break }
}
if (is.null(pick)) stop("no representative rep")
fit <- pick$fit; rf <- pick$rf; X <- pick$X
cat(sprintf("rep seed=%d | KNOWN p_path=%.4f->p_union=%.4f (x%.0f) | F p_path=%.4f->p_union=%.4f (x%.2f) R_obs=%.1f\n",
            pick$s, fit$pval_path, fit$pval_union, fit$pval_path/fit$pval_union,
            rf$p_path, rf$p_union, rf$p_path/rf$p_union, rf$R_obs))

# GRAMMAR EXCEPTION (documented): FILL here encodes REGIONS of a truncation set
# (path set / union gain / union extra), not a method or variance series -- there
# is no inference series to confuse, so semantic region fills are allowed. Rose =
# path set, blue = the union's extra mass, grey = union mass that does not help.
side_lab <- c(path = "path set", below = "union gain (below obs)", above = "union extra (above obs, no help)")
side_fill <- c("path set" = "#C48A97",
               "union gain (below obs)" = "#1F5AA6",
               "union extra (above obs, no help)" = "#BBBBBB")
classify <- function(grid, dens, Mpath, obs) {
  cl <- rep(NA_character_, length(grid))
  cl[dens > 0 & in_set(grid, Mpath)] <- "path set"
  ex <- dens > 0 & is.na(cl)
  cl[ex & grid <  obs] <- "union gain (below obs)"
  cl[ex & grid >= obs] <- "union extra (above obs, no help)"
  factor(cl, levels = side_lab)
}

## ---- Panel A : known variance (chi^2_q over the phi-ray) ---------------------
sc <- fit$scale_factor; uobs <- fit$test_stat^2 / sc
to_u <- function(M) { m <- cbind(M[,1]^2/sc, M[,2]^2/sc); m }
U_un <- to_u(as.matrix(fit$interval_union)); U_pa <- to_u(as.matrix(fit$interval_path))
gA <- seq(min(U_un) - 4, max(U_un) + 4, length.out = 3000)
dA <- cond_density(gA, dchisq(gA, q, log = TRUE), U_un)
datA <- data.frame(x = gA, d = dA, side = classify(gA, dA, U_pa, uobs))
pA <- ggplot(subset(datA, d > 0), aes(x, d, fill = side)) +
  geom_col(width = diff(gA)[1], position = "identity") +
  geom_vline(xintercept = uobs, colour = km_ref, linetype = 2, linewidth = 0.6) +
  annotate("text", x = uobs, y = max(dA)*0.96, label = "observed", hjust = 1.06, size = 3, colour = km_ref) +
  scale_fill_manual(values = side_fill, name = NULL, drop = FALSE) +
  labs(x = expression(paste("test coordinate  ", phi^2, "/", sigma^2, "   (ref ", chi[q]^2, ")")),
       y = "conditional reference density",
       title = sprintf("Known variance: union extends DOWN  (p: %.3f -> %.4f)",
                       fit$pval_path, fit$pval_union)) +
  theme_km() + theme(legend.position = "none")

## ---- Panel B : unknown variance (F_{q,(m-2)q} over the R-fiber) --------------
d1 <- q; d2 <- (nrow(X) - 2) * q
R_un0 <- rf$intervals_r
Rcap <- 3 * rf$R_obs                                          # F-mass past 3*R_obs is negligible
R_un <- R_un0[R_un0[,1] < Rcap, , drop = FALSE]              # drop arcs entirely beyond the cap
R_un[R_un[,2] > Rcap, 2] <- Rcap                              # truncate the runaway upper arc
cont <- R_un[,1] <= rf$R_obs & rf$R_obs <= R_un[,2]
R_pa <- if (any(cont)) R_un[cont, , drop = FALSE] else R_un[which.min(abs(rowMeans(R_un) - rf$R_obs)), , drop = FALSE]
gB <- seq(min(R_un) - 1, Rcap, length.out = 3000)
dB <- cond_density(gB, df(gB, d1, d2, log = TRUE), R_un)
datB <- data.frame(x = gB, d = dB, side = classify(gB, dB, R_pa, rf$R_obs))
pB <- ggplot(subset(datB, d > 0), aes(x, d, fill = side)) +
  geom_col(width = diff(gB)[1], position = "identity") +
  geom_vline(xintercept = rf$R_obs, colour = km_ref, linetype = 2, linewidth = 0.6) +
  annotate("text", x = rf$R_obs, y = max(dB)*0.96, label = "observed R", hjust = -0.06, size = 3, colour = km_ref) +
  scale_fill_manual(values = side_fill, name = NULL, drop = FALSE) +
  labs(x = expression(paste("studentized  R   (ref ", F[paste(q, ",", (m-2)*q)], ")")),
       y = "conditional reference density",
       title = sprintf("Unknown variance: room all ABOVE obs  (p: %.3f -> %.3f, R_obs=%.0f)",
                       rf$p_path, rf$p_union, rf$R_obs)) +
  theme_km() + theme(legend.position = "none")

leg  <- km_get_legend(pA + theme(legend.position = "bottom"))
body <- arrangeGrob(pA, pB, ncol = 2,
                    top = textGrob("The union helps only when it adds mass BELOW the observed statistic",
                                   gp = gpar(fontface = "bold", fontsize = 12.5)))
g <- arrangeGrob(body, leg, ncol = 1, heights = c(10, 1))
ggsave_km(g, "sims/results/studentized_away", width = 11.5, height = 5.2)
cat("Wrote sims/results/studentized_away.{pdf,png}\n")
