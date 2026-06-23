# ---------------------------------------------------------------------------
# plot_power_geometries.R -- known-sigma union vs path power across cluster
# geometries (Yun-Barber Fig 4 style). Power vs separation delta, faceted by
# geometry; union dashed, path solid; 95% CI. Reads power_geometries.rds.
# ---------------------------------------------------------------------------
source("sims/house_style.R")
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

df <- as.data.frame(readRDS("sims/results/power_geometries.rds"))

summ <- df %>%
  pivot_longer(c(p_path, p_union), names_to = "cond", values_to = "p") %>%
  mutate(cond = factor(sub("^p_", "", cond), levels = c("path", "union"))) %>%
  group_by(geom, delta, cond) %>%
  summarise(n = sum(is.finite(p)), pow = mean(p < 0.05, na.rm = TRUE), .groups = "drop") %>%
  mutate(se = sqrt(pmax(pow*(1-pow), 0)/pmax(n,1)),
         lo = pmax(pow - 1.96*se, 0), hi = pmin(pow + 1.96*se, 1),
         geom = factor(geom, levels = c("K=2 horizontal", "K=3 triangle", "K=3 collinear")))

fig <- ggplot(summ, aes(delta, pow, linetype = cond, group = cond)) +
  geom_line(linewidth = 0.85, colour = km_col[["oracle"]]) +
  geom_point(size = 1.4, colour = km_col[["oracle"]]) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.4, colour = km_col[["oracle"]]) +
  facet_wrap(~ geom, nrow = 1) +
  scale_linetype_manual(values = c(path = 1, union = 2), labels = c(path = "Path", union = "Union"), name = NULL) +
  labs(x = expression(separation~~delta), y = "power") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_km() + theme(legend.key.width = unit(1.4, "cm"))

ggsave_km(fig, "sims/results/power_geometries", width = 11, height = 3.9)
cat("Wrote power_geometries.{pdf,png}\n")
