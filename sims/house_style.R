# ---------------------------------------------------------------------------
# house_style.R -- shared ggplot aesthetics + the UNIFIED VISUAL GRAMMAR.
# (Matches ../distributional-ppi: despined, light dashed grid, bold labels,
#  PNG@400dpi + selectable-text PDF, ONE shared bottom legend, NO subtitles.)
#
# ====================  UNIFIED GRAMMAR (apply everywhere)  ===================
#   COLOUR   = inference treatment / variance handling  (km_col)
#   LINETYPE = conditioning level: path = solid, union = dashed  (km_cond)
#   FACET    = the swept dimension (q, sigma, geometry) -- NEVER colour
#   SHAPE    = avoid on lines (redundant with linetype/colour). EXCEPTION on POINT
#              geoms where linetype cannot apply: dot-plots & real-data panel B use
#              shape for path=open(1)/union=filled(16), and real-data scatters use
#              shape for the cluster/species label. These are the only shape uses.
#   grey #4D4D4D (km_ref) = REFERENCE LINES ONLY (45-deg, null, zero, alpha);
#                           never a data series.
#
# "Union" is the conditioning concept (dashed). The studentized union is computed
# on the R-fiber, the known-sigma union on the phi-ray -- the fiber is implied by
# the COLOUR (known vs studentized), so "R-fiber" is a methods footnote, not a
# label. path solid / union dashed is FIXED across every figure.
# ---------------------------------------------------------------------------
suppressMessages({ library(ggplot2); library(gridExtra); library(grid) })

km_ref <- "#4D4D4D"                                  # reference / null lines ONLY

## COLOUR = treatment (variance handling). One hex per role, used paper-wide.
km_col <- c(
  naive       = "#E69F00",   # ordinary non-selective p-value (invalid) -- motivation/Type-I only
  oracle      = "#000000",   # known / true sigma  (also colour for known-sigma selective tests)
  studentized = "#1F5AA6",   # exact unknown-sigma (studentized-F / R-fiber) -- blue, the proposal
  med         = "#009E73",   # sigma_MED plug-in -- green
  sample      = "#D55E00")   # sigma_all / sample plug-in -- vermillion
km_col_lab <- c(                                     # ASCII "sigma" -> renders in PNG and PDF
  naive = "Naive (invalid)", oracle = "Known / oracle sigma",
  studentized = "Studentized-F (unknown sigma)",
  med = "sigma-MED plug-in", sample = "sigma-sample plug-in")

## LINETYPE = conditioning level. path solid, union dashed -- ALWAYS.
km_cond     <- c(path = 1, union = 2)
km_cond_lab <- c(path = "Path", union = "Union")

## convenience scales (drop-in)
scale_km_colour   <- function(...) scale_colour_manual(values = km_col, labels = km_col_lab, name = NULL, ...)
scale_km_fill     <- function(...) scale_fill_manual(values = km_col, labels = km_col_lab, name = NULL, ...)
scale_km_linetype <- function(...) scale_linetype_manual(values = km_cond, labels = km_cond_lab, name = NULL, ...)

## DEPRECATED (pre-grammar) -- kept so any un-refactored script still runs.
## Do NOT use in new/refactored figures: blue now means studentized, not union.
km_pal <- c(union = "#1F5AA6", path = "#C48A97", naive = "#E69F00")
km_lab <- c(union = "Union (proposed)", path = "Path (original)", naive = "Naive (invalid)")

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
      legend.key.width  = unit(1.5, "cm"),          # wide enough to SHOW the path/union dash
      legend.background = element_blank(),
      legend.margin     = margin(2, 2, 2, 2),
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold"),
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
