#!/usr/bin/env Rscript
# Union vs path power, KNOWN vs UNKNOWN variance (house style, no subtitle).
source("sims/house_style.R")

res <- readRDS("sims/results/unknown_var_proto.rds")
alpha <- 0.05
rej <- function(p) mean(p <= alpha, na.rm = TRUE)

agg <- do.call(rbind, lapply(split(res, res$delta), function(d) data.frame(
  delta = d$delta[1],
  union_known = rej(d$p_union_known), path_known = rej(d$p_path_known),
  union_unknown = rej(d$p_union_unk), path_unknown = rej(d$p_path_unk))))

long <- with(agg, rbind(
  data.frame(delta, rej = union_known,   method = "union", variance = "known variance"),
  data.frame(delta, rej = path_known,    method = "path",  variance = "known variance"),
  data.frame(delta, rej = union_unknown, method = "union", variance = "unknown variance"),
  data.frame(delta, rej = path_unknown,  method = "path",  variance = "unknown variance")))
long$method   <- factor(long$method,   levels = c("union", "path"))
long$variance <- factor(long$variance, levels = c("known variance", "unknown variance"))

p <- ggplot(long, aes(delta, rej, colour = method, linetype = variance, shape = variance)) +
  geom_hline(yintercept = alpha, linetype = 2, colour = km_ref, alpha = 0.6) +
  geom_line(linewidth = 1.0) + geom_point(size = 2.8) +
  scale_colour_manual(values = km_pal, labels = km_lab) +
  scale_linetype_manual(values = c("known variance" = 1, "unknown variance" = 2)) +
  scale_shape_manual(values = c("known variance" = 16, "unknown variance" = 17)) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0.02)) +
  guides(colour = guide_legend(order = 1), linetype = guide_legend(order = 2),
         shape = guide_legend(order = 2)) +
  labs(x = expression(paste("separation  ", delta)), y = "rejection rate",
       title = "Power under known vs unknown variance") +
  theme_km()

ggsave_km(p, "sims/results/unknown_var", width = 7, height = 5)
cat("Wrote sims/results/unknown_var.{pdf,png}\n")
