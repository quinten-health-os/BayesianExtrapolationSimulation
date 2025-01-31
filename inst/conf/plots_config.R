# Set up constants
WIDTH <- 468
LINEWIDTH <- 0.5
ADULT_BLUE <- c(0, 0, 128, 255 * 0.8) / 255 #' 000080cc'
CHILD_RED <- c(128, 0, 0, 255 * 0.8) / 255 #' 800000cc'
cmap <- viridis::viridis(255)
markers_list <- c("o", "s", "D", "^", "v", "p", "*", "+", "x")
markersize <- 1
relative_error_cap_width <- 0.05
dpi <- 300
text_size <- 14
small_text_size <- 12
width <- 11
height <- 8.5
plots_to_latex <<- FALSE
textwidth <- 426

font <- "CMU Serif"

# Set up style
STYLE <- "paper"
if (STYLE == "paper") {
  tex_fonts <- list(
    # Use LaTeX to write all text
    text.usetex = TRUE,
    font.family = "serif",
    font.serif = c("Computer Modern Roman"),
    # Computer modern roman 10 cmr10
    # Use 10pt font in plots, to match 10pt font in document
    axes.labelsize = 10,
    axes.titlesize = 10,
    font.size = 10,
    # Make the legend/label fonts a little smaller
    legend.fontsize = 8,
    xtick.labelsize = 8,
    ytick.labelsize = 8
  )
} else if (STYLE == "presentation") {
  tex_fonts <- list(
    # Use LaTeX to write all text
    text.usetex = TRUE,
    font.family = "serif",
    font.serif = c("Computer Modern Roman"),
    # [CMU_SERIF_FONT_PATH],
    # Use 10pt font in plots, to match 10pt font in document
    axes.labelsize = 14,
    axes.titlesize = 14,
    font.size = 14,
    # Make the legend/label fonts a little smaller
    legend.fontsize = 14,
    xtick.labelsize = 14,
    ytick.labelsize = 14
  )
}

library(ggplot2)
ggplot2::theme_set(theme_minimal() + theme(text = element_text(family = "serif")))
options(plot.title = expression(paste(bold(""))))
options(plot.subtitle = expression(paste(bold(""))))
options(plot.caption = expression(paste(bold(""))))
options(plot.margin = unit(c(1, 1, 1, 1), "lines"))
options(title.padding = unit(4, "lines"))
options(axis.title.y = expression(atop(bold(""), "")))
options(axis.title.x = expression(atop(bold(""), "")))
options(legend.position = "bottom")
options(legend.title = expression(paste(bold(""))))
options(legend.text = expression(paste(bold(""))))
options(legend.key.size = unit(1, "lines"))
options(legend.key.height = unit(1, "lines"))
options(legend.key.width = unit(1, "lines"))
options(legend.background = element_rect(fill = "white", color = NA))
options(strip.background = element_rect(fill = "white", color = NA))
options(strip.text = expression(paste(bold(""))))
options(strip.text.x = expression(paste(bold(""))))
options(strip.text.y = expression(paste(bold(""))))
options(panel.background = element_rect(fill = "white", color = NA))
options(panel.grid.major = element_line(
  color = "grey90",
  size = 0.25,
  linetype = "solid"
))
options(panel.grid.minor = element_line(
  color = "grey95",
  size = 0.15,
  linetype = "solid"
))
options(panel.border = element_rect(fill = NA, color = "black"))
options(panel.spacing = unit(1, "lines"))
options(axis.line = element_line(
  color = "black",
  size = 0.5,
  linetype = "solid"
))
options(axis.ticks = element_line(
  color = "black",
  size = 0.5,
  linetype = "solid"
))
options(axis.ticks.length = unit(0.2, "lines"))
options(axis.ticks.margin = unit(0.2, "lines"))
options(axis.text = element_text(family = "serif", color = "black"))
options(axis.text.x = element_text(
  angle = 45,
  hjust = 1,
  vjust = 1,
  size = 8
))
options(axis.text.y = element_text(
  hjust = 1,
  vjust = 1,
  size = 8
))
options(plot.title = element_text(
  family = "serif",
  face = "bold",
  size = 14,
  hjust = 0.5,
  margin = margin(t = 10)
))
options(plot.subtitle = element_text(
  family = "serif",
  face = "bold",
  size = 12,
  hjust = 0.5,
  margin = margin(t = 5)
))
options(plot.caption = element_text(
  family = "serif",
  face = "italic",
  size = 8,
  hjust = 0
))

# Set up x-variables
xvars <- list(
  control_drift = list(name = "control_drift", label = "Control drift", fmt = "o"),
  drift = list(name = "drift", label = "Drift in treatment effect", fmt = "-"),
  method_param = list(name = "method_param", label = "Method parameter", fmt = "o"),
  TIE = list(name = "TIE", label = "TIE rate", fmt = "o"),
  sample_size = list(name = "sample_size", label = "Target study sample size per arm", fmt = "o")
)
