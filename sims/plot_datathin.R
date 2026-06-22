#!/usr/bin/env Rscript
# plot_datathin.R -- data thinning vs selective k-means inference, 3 metrics, in
# the UNIFIED GRAMMAR. COLOUR = sigma-handling (km_col); LINETYPE = paradigm
# (selective = solid, data thinning = dashed -- the path/union slot, reused for
# the two camps in this comparison figure).
#   A  Type I QQ vs Uniform (global null)         -- valid tests track the 45-deg line.
#   B  detection probability vs delta             -- full-data clustering vs thinned splits.
#   C  conditional power vs delta                 -- reject | the pair was recovered.
# Reads sims/results/datathin_compare.rds (from exp_datathin.R).
source("sims/house_style.R")
suppressMessages({ library(ggplot2); library(gridExtra); library(grid) })

d <- readRDS("sims/results/datathin_compare.rds")
nullc <- d[d$is_null, ]; sigc <- d[!d$is_null, ]
DELTAS <- sort(unique(sigc$delta)); A <- 0.05
wilson <- function(x, n) { if (n == 0) return(c(NA, NA)); p <- x/n; z <- 1.96
  ctr <- (p + z^2/(2*n))/(1 + z^2/n); hw <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/(1 + z^2/n); c(ctr-hw, ctr+hw) }

## grammar: colour = treatment (km_col), linetype = paradigm
trt_lab <- c(naive = "Naive (invalid)", oracle = "Known / oracle sigma",
             studentized = "Studentized-F (unknown sigma)", med = "sigma-MED", sample = "sigma-sample")
par_lt  <- c(selective = 1, thinning = 2)
par_lab <- c(selective = "Selective (fully conditional)", thinning = "Data thinning")
trtlev  <- names(trt_lab); parlev <- names(par_lt)

## method registry: key, p-col, detect-col, treatment(colour), paradigm(linetype)
M <- list(
  list(k="naive", p="p_naive",       det="sel_detected",   trt="naive",       par="selective"),
  list(k="union", p="p_union",       det="sel_detected",   trt="oracle",      par="selective"),
  list(k="studF", p="p_studF",       det="sel_detected",   trt="studentized", par="selective"),
  list(k="thinO", p="p_thin_oracle", det="thinO_detected", trt="oracle",      par="thinning"),
  list(k="thinM", p="p_thin_med",    det="thinM_detected", trt="med",         par="thinning"),
  list(k="thinS", p="p_thin_samp",   det="thinS_detected", trt="sample",      par="thinning"))
fct <- function(df) { df$trt <- factor(df$trt, trtlev); df$par <- factor(df$par, parlev); df }

## ---- Panel A: Type I QQ (global null) --------------------------------------
qqd <- do.call(rbind, lapply(M, function(m) { p <- sort(nullc[[m$p]][is.finite(nullc[[m$p]])])
  if (length(p) < 20) return(NULL)
  data.frame(theo = ppoints(length(p)), emp = p, trt = m$trt, par = m$par, k = m$k) }))
pA <- ggplot(fct(qqd), aes(theo, emp, colour = trt, linetype = par, group = k)) +
  geom_abline(slope = 1, colour = km_ref, alpha = 0.6, linewidth = 0.4) +
  geom_step(linewidth = 0.85) +
  scale_colour_manual(values = km_col, labels = trt_lab, name = NULL, drop = FALSE) +
  scale_linetype_manual(values = par_lt, labels = par_lab, name = NULL, drop = FALSE) +
  scale_x_continuous(breaks = c(0, .5, 1)) + scale_y_continuous(breaks = c(0, .5, 1)) +
  labs(x = "Uniform quantile", y = "p-value (global null)", title = "A.  Type I calibration") +
  theme_km() + theme(aspect.ratio = 1, legend.position = "none", plot.title = element_text(size = 11))

## ---- Panel B: detection probability vs delta -------------------------------
# detection depends on the CLUSTERING source: full data (selective, all methods share
# it) vs each thinned split (oracle / MED / sample sigma into the split).
detsrc <- list(list(k="sel", det="sel_detected", trt="oracle", par="selective"),
               list(k="thinO", det="thinO_detected", trt="oracle", par="thinning"),
               list(k="thinM", det="thinM_detected", trt="med",    par="thinning"),
               list(k="thinS", det="thinS_detected", trt="sample", par="thinning"))
detdf <- do.call(rbind, lapply(detsrc, function(s) do.call(rbind, lapply(DELTAS, function(dl) {
  cc <- sigc[sigc$delta == dl, ]; det <- cc[[s$det]]; n <- length(det); x <- sum(det); ci <- wilson(x, n)
  data.frame(delta = dl, rate = x/n, lo = ci[1], hi = ci[2], trt = s$trt, par = s$par, k = s$k) }))))
pB <- ggplot(fct(detdf), aes(delta, rate, colour = trt, linetype = par, group = k)) +
  geom_errorbar(aes(ymin = pmax(0, lo), ymax = pmin(1, hi)), width = 0.12, linewidth = 0.4, linetype = 1) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.7) +
  scale_colour_manual(values = km_col, guide = "none", drop = FALSE) +
  scale_linetype_manual(values = par_lt, guide = "none", drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "detection probability", title = "B.  Detection probability") +
  theme_km() + theme(legend.position = "none", plot.title = element_text(size = 11))

## ---- Panel C: conditional power vs delta -----------------------------------
cpkeys <- c("union", "studF", "thinO", "thinM", "thinS")
cpdf <- do.call(rbind, lapply(M[sapply(M, function(m) m$k %in% cpkeys)], function(m)
  do.call(rbind, lapply(DELTAS, function(dl) {
    cc <- sigc[sigc$delta == dl, ]; det <- cc[[m$det]]; p <- cc[[m$p]]
    sel <- det & is.finite(p); n <- sum(sel); x <- sum(p[sel] <= A); ci <- wilson(x, n)
    data.frame(delta = dl, rate = if (n) x/n else NA, lo = ci[1], hi = ci[2], trt = m$trt, par = m$par, k = m$k) }))))
pC <- ggplot(fct(cpdf), aes(delta, rate, colour = trt, linetype = par, group = k)) +
  geom_hline(yintercept = A, colour = km_ref, alpha = 0.6, linewidth = 0.4) +
  geom_errorbar(aes(ymin = pmax(0, lo), ymax = pmin(1, hi)), width = 0.12, linewidth = 0.4, linetype = 1) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.7) +
  scale_colour_manual(values = km_col, guide = "none", drop = FALSE) +
  scale_linetype_manual(values = par_lt, guide = "none", drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "conditional power  (reject | detected)",
       title = "C.  Conditional power") +
  theme_km() + theme(legend.position = "none", plot.title = element_text(size = 11))

leg  <- km_get_legend(pA + theme(legend.position = "bottom"))
body <- arrangeGrob(pA, pB, pC, ncol = 3, widths = c(1.0, 1.0, 1.0))
g <- arrangeGrob(body, leg, ncol = 1, heights = c(10, 1.3),
                 top = textGrob("Data thinning vs selective k-means inference   (k=3, q=2, n=30, sigma=1)",
                                gp = gpar(fontface = "bold", fontsize = 12)))
ggsave_km(g, "sims/results/datathin", width = 15.5, height = 5.4)
cat("Wrote sims/results/datathin.{pdf,png}\n")
