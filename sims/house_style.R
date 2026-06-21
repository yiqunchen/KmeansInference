# ---------------------------------------------------------------------------
# house_style.R -- shared ggplot aesthetics, matching ../distributional-ppi
# (src/distributional_ppi/plotting.py: apply_publication_style()).
#
#   * despined (no top/right), light dashed grid, bold axis labels/ticks
#   * legend has no frame; multi-panel figures use ONE shared bottom legend
#   * ONE concise bold title per panel; NO subtitles
#   * export PNG @ 400 dpi + PDF with selectable text, tight margins
#
# Semantic palette (proposed = blue, baseline = rose, naive = orange):
#   union  #1F5AA6   (their "augmented"/proposed)
#   path   #C48A97   (their "labeled only"/baseline)
#   naive  #E69F00   (their "naive proxy")
# ---------------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(gridExtra); library(grid) })

km_pal <- c(union = "#1F5AA6", path = "#C48A97", naive = "#E69F00")
km_ref <- "#4D4D4D"                                  # reference / null lines
km_lab <- c(union = "Union (proposed)", path = "Path (original)",
            naive = "Naive (invalid)")

theme_km <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title  = element_text(face = "bold"),
      axis.text   = element_text(face = "bold", colour = "black"),
      plot.title  = element_text(face = "bold", hjust = 0, size = rel(1.0)),
      axis.line   = element_line(linewidth = 0.7, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.6, colour = "black"),
      panel.grid.major = element_line(linewidth = 0.3, linetype = "dashed",
                                      colour = adjustcolor(km_ref, alpha.f = 0.22)),
      panel.grid.minor = element_blank(),
      legend.position   = "bottom",
      legend.title      = element_blank(),
      legend.key        = element_blank(),
      legend.background = element_blank(),
      legend.margin     = margin(2, 2, 2, 2),
      plot.title.position = "plot",
      plot.margin = margin(6, 10, 6, 6))
}

# extract the legend grob from a ggplot (for a single shared legend)
km_get_legend <- function(p) {
  g <- ggplotGrob(p + theme(legend.position = "bottom"))
  idx <- which(vapply(g$grobs, function(x) x$name, "") == "guide-box")
  if (length(idx) == 0) return(nullGrob())
  g$grobs[[idx[1]]]
}

# arrange panels (legend.position='none' each) under ONE shared bottom legend
km_panels <- function(plots, legend, ncol = length(plots), heights = c(10, 1)) {
  body <- arrangeGrob(grobs = plots, ncol = ncol)
  arrangeGrob(body, legend, ncol = 1, heights = heights)
}

# save a grob/plot to PDF (selectable text) + PNG @ >=300 dpi
ggsave_km <- function(plot, prefix, width = 7, height = 4.6, dpi = 400) {
  ggsave(paste0(prefix, ".pdf"), plot, width = width, height = height,
         device = "pdf", useDingbats = FALSE)
  ggsave(paste0(prefix, ".png"), plot, width = width, height = height, dpi = dpi)
  invisible(prefix)
}
