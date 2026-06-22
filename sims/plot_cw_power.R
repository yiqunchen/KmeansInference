#!/usr/bin/env Rscript
# Fig 2 -- Power by dimension q (n=90, sigma=1), facet_grid(variance ~ q):
#   rows = {known, unknown} variance; union vs path; 95% binomial CI bars.
# Fig 3 -- Detection probability + conditional power vs delta (Chen & Witten
#   Fig. 5 estimands), by q, union vs path.
# Both read the VALID columns (known: rej_union/rej_path; unknown R-fiber:
# rej_union_rfib/rej_path_rfib). Robust to partial data (drops cells < 40 reps).
source("sims/house_style.R")
suppressMessages({ library(gridExtra); library(grid) })
DIR <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sims/results/sweep_cw"
s   <- read.csv(file.path(DIR, "summary.csv"), stringsAsFactors = FALSE)
pw  <- s[grepl("cw_power", s$label) & s$n_valid >= 40, ]
if (nrow(pw) == 0) stop("no power cells with >=40 valid reps yet")
se  <- function(p, n) sqrt(pmax(p * (1 - p), 0) / n)
qlev <- paste0("q = ", c(2, 10, 50))
qf <- function(q) factor(paste0("q = ", q), levels = qlev)

## ---- Fig 2: power vs delta, facet_grid(variance ~ q) -----------------------
## GRAMMAR: COLOUR = treatment (known row -> oracle black, unknown row ->
## studentized blue); LINETYPE = conditioning (path solid, union dashed); q on
## facet. The treatment colour is constant within a row, so the per-row colour
## doubles as the variance-handling cue while linetype carries path-vs-union.
mk <- function(uy, py, v, treat) rbind(
  data.frame(q = pw$q, delta = pw$delta, n = pw$n_valid, rej = pw[[uy]], cond = "union", treat = treat, variance = v),
  data.frame(q = pw$q, delta = pw$delta, n = pw$n_valid, rej = pw[[py]], cond = "path",  treat = treat, variance = v))
## KNOWN-sigma only: the union's power gain is a known-variance result. The
## unknown-variance comparisons (which require a different fiber and so cannot be
## stacked against the known-sigma curves without inviting a false "studentized >
## known" cross-fiber read) live in repro_yunbarber + cw_path_known_vs_unknown.
d <- mk("rej_union", "rej_path", "known variance", "oracle")
d$q <- qf(d$q)
d$cond <- factor(d$cond, levels = c("path", "union"))
d$treat <- factor(d$treat, levels = c("oracle", "studentized"))
d$lo <- pmax(0, d$rej - 1.96 * se(d$rej, d$n)); d$hi <- pmin(1, d$rej + 1.96 * se(d$rej, d$n))
p_pow <- ggplot(d, aes(delta, rej, colour = treat, linetype = cond)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, linewidth = 0.4, alpha = 0.7) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.9) +
  facet_wrap(~ q, nrow = 1) +
  scale_km_colour() + scale_km_linetype() +
  scale_x_continuous(breaks = seq(2, 8, 2)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "power",
       title = expression(paste("Union vs path power by dimension q  (known ", sigma, ", n=90, 95% CI)"))) +
  theme_km() + theme(strip.background = element_blank(),
                     strip.text = element_text(face = "bold"), panel.spacing = unit(0.6, "lines"))
ggsave_km(p_pow, file.path(DIR, "../cw_power"), width = 11, height = 3.6)

## ---- Fig 3: detection probability + conditional power (known sigma) ---------
## with 95% binomial CIs: detection ~ n_valid; conditional power ~ n_detected.
sebin <- function(p, n) sqrt(pmax(p*(1-p), 0)/pmax(n, 1))
det <- pw[pw$delta <= 8, ]; det$qf <- qf(det$q)
det$dlo <- pmax(0, det$detect_prob - 1.96*sebin(det$detect_prob, det$n_valid))
det$dhi <- pmin(1, det$detect_prob + 1.96*sebin(det$detect_prob, det$n_valid))
## p_det: detection probability has NO treatment/conditioning series -- a single
## neutral series faceted by q (q must NOT be a colour). Reference grey is for
## lines only, so detection uses the neutral data ink (black) with q on facet.
p_det <- ggplot(det, aes(delta, detect_prob)) +
  geom_errorbar(aes(ymin = dlo, ymax = dhi), width = 0.18, linewidth = 0.4, alpha = 0.7) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  facet_wrap(~ qf, nrow = 1) +
  scale_x_continuous(breaks = seq(1, 8, 2)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "detection probability",
       title = "Detection (95% CI)") + theme_km()
## p_cp: known-sigma conditional power -> COLOUR = treatment (oracle black),
## LINETYPE = conditioning (path solid, union dashed), q on facet.
cp <- rbind(
  data.frame(delta = det$delta, q = det$qf, cp = det$cpow_union, nd = det$n_detected, cond = "union"),
  data.frame(delta = det$delta, q = det$qf, cp = det$cpow_path,  nd = det$n_detected, cond = "path"))
cp$cond <- factor(cp$cond, levels = c("path", "union")); cp$treat <- factor("oracle", levels = "oracle")
cp <- cp[is.finite(cp$cp), ]
cp$lo <- pmax(0, cp$cp - 1.96*sebin(cp$cp, cp$nd)); cp$hi <- pmin(1, cp$cp + 1.96*sebin(cp$cp, cp$nd))
p_cp <- ggplot(cp, aes(delta, cp, colour = treat, linetype = cond)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.4, alpha = 0.6) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  facet_wrap(~ q, nrow = 1) +
  scale_km_colour() + scale_km_linetype() +
  scale_x_continuous(breaks = seq(1, 8, 2)) + scale_y_continuous(limits = c(0, 1)) +
  labs(x = expression(paste("separation  ", delta)), y = "conditional power",
       title = "Conditional power, known sigma (union vs path, 95% CI)") + theme_km()
fig3 <- arrangeGrob(p_det, p_cp, ncol = 2)
ggsave_km(fig3, file.path(DIR, "../cw_detection"), width = 12, height = 4.6)
cat("Wrote cw_power.{pdf,png} + cw_detection.{pdf,png} (in", dirname(DIR), ")\n")
