# List of non-standard packages used in project scripts
required_pkgs <- c(
  "akima",         # interpolation for horizontal slices
  "fields",        # image.plot, color bars, etc.
  "RColorBrewer",  # colour palettes
  "latex2exp"      # LaTeX-style labels in flyby_plots.R
)

# Install any that are missing
installed <- rownames(installed.packages())
to_install <- setdiff(required_pkgs, installed)

if (length(to_install) > 0) {
  install.packages(to_install)
} else {
  message("All required packages are already installed.")
}
