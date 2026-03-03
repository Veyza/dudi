# Shared helper: draw a filled-contour plot on a rectangular grid
.plot_plane <- function(z, x_range, y_range, levels, cols,
                        filename, main, key_labels = NULL,
                        xlab = "km", ylab = "km", draw_circle = FALSE) {
  nx <- nrow(z); ny <- ncol(z)
  x_seq <- seq(x_range[1], x_range[2], length.out = nx)
  y_seq <- seq(y_range[1], y_range[2], length.out = ny)

  png(filename, width = 900, height = 800)
  par(oma = c(2.0, 3.0, 0.0, 0.5), mar = c(4, 5, 4, 1.5), cex.axis = 2)

  graphics::filled.contour(
    x = x_seq, y = y_seq, z = z,
    levels = levels,
    col    = cols,
    plot.title = {
      par(cex.main = 2.5, cex.lab = 2.5)
      title(main = main, xlab = xlab, ylab = ylab)
    },
    plot.axes = {
      axis(1); axis(2)
      if (draw_circle) {
        symbols(x = -3, y = 0, circles = 252, inches = FALSE, add = TRUE,
                bg = "darkgrey", fg = "darkgrey")
      }
    },
    key.axes  = {
      if (!is.null(key_labels)) {
        axis(4, at = levels, labels = key_labels, cex.axis = 1.5)
      } else {
        axis(4, cex.axis = 1.5)
      }
    }
  )

  dev.off()
}

## 1) Composition: salt_rich / (salt_rich + salt_poor) -----------------------

composition_in_plane <- function(close = FALSE) {
  # Plot salt-rich fraction in the E0 vertical plane

  library(RColorBrewer)

  if (close) {
    salt_poor_file <- "./results/salt_poor_0_plane_close.dat"
    salt_rich_file <- "./results/salt_rich_0_plane_close.dat"
    x_range <- c(-150, 150); y_range <- c(237, 372)
    out_file <- "./results/type3_proportion_close.png"
  } else {
    salt_poor_file <- "./results/salt_poor_0_plane.dat"
    salt_rich_file <- "./results/salt_rich_0_plane.dat"
    x_range <- c(-450, 450); y_range <- c(0, 900)
    out_file <- "./results/type3_proportion.png"
  }

  sp <- as.matrix(read.table(salt_poor_file))
  sr <- as.matrix(read.table(salt_rich_file))
  s_total <- sp + sr

  frac <- matrix(0, nrow = nrow(sp), ncol = ncol(sp))
  mask <- s_total > 1e-22
  frac[mask] <- sr[mask] / s_total[mask] * 100.0

  if (close) {
    lev <- c(5, 10, 17, seq(30, 90, by = 10))
  } else {
    lev <- seq(5, 75, by = 5)
  }
  nbin <- length(lev) - 1
  cols <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(nbin)

  .plot_plane(
    z           = frac,
    x_range     = x_range,
    y_range     = y_range,
    levels      = lev,
    cols        = cols,
    filename    = out_file,
    main        = "salt-rich dust fraction",
    key_labels  = lev,
    xlab        = "km",
    ylab        = "km",
    draw_circle = TRUE
  )
}

## 2) Number density: salt_poor + salt_rich ----------------------------------

number_density_in_plane <- function(close = FALSE) {
  # Plot total number density in the E0 vertical plane

  library(RColorBrewer)

  if (close) {
    salt_poor_file <- "./results/salt_poor_0_plane_close.dat"
    salt_rich_file <- "./results/salt_rich_0_plane_close.dat"
    x_range <- c(-150, 150); y_range <- c(237, 372)
    out_file <- "./results/number_density_E0_close.png"
  } else {
    salt_poor_file <- "./results/salt_poor_0_plane.dat"
    salt_rich_file <- "./results/salt_rich_0_plane.dat"
    x_range <- c(-450, 450); y_range <- c(0, 900)
    out_file <- "./results/number_density_E0.png"
  }

  sp <- as.matrix(read.table(salt_poor_file))
  sr <- as.matrix(read.table(salt_rich_file))
  s  <- sp + sr

  # Your preferred levels; plotted on log10 scale, labels in linear units
  lev_lin <- c(0.01, 0.1, 0.3, 1, 10, 100, 200, 300, 500,
               1000, 2000, 3000, 4000)
  min_lev <- min(lev_lin)

  s_pos <- pmax(s, min_lev * 1e-2)
  z     <- log10(s_pos)
  lev   <- log10(lev_lin)

  nbin <- length(lev) - 1
  cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(nbin))

  .plot_plane(
    z           = z,
    x_range     = x_range,
    y_range     = y_range,
    levels      = lev,
    cols        = cols,
    filename    = out_file,
    main        = "number density in E0 flyby plane",
    key_labels  = lev_lin,
    xlab        = "km",
    ylab        = "km",
    draw_circle = TRUE
  )
}

## 3) Dust-to-gas ratio: dust / gas ------------------------------------------

dust2gas_in_plane <- function(close = FALSE) {
  # Plot dust/gas ratio using vertical_structure outputs:
  # dust: fnum = 0.1 (dens_in_0_plane_dust*)
  # gas : fnum = 0.2 (dens_in_0_plane_gas*)

  library(RColorBrewer)

  if (close) {
    dust_file <- "./results/dens_in_0_plane_dust_close.dat"
    gas_file  <- "./results/dens_in_0_plane_gas_close.dat"
    x_range <- c(-150, 150); y_range <- c(237, 372)
    out_file <- "./results/dust2gas_E0_close.png"
  } else {
    dust_file <- "./results/dens_in_0_plane_dust.dat"
    gas_file  <- "./results/dens_in_0_plane_gas.dat"
    x_range <- c(-450, 450); y_range <- c(0, 900)
    out_file <- "./results/dust2gas_E0.png"
  }

  dust <- as.matrix(read.table(dust_file))
  gas  <- as.matrix(read.table(gas_file))

  # ratio only where gas > 0, NA elsewhere
  ratio <- matrix(NA_real_, nrow = nrow(dust), ncol = ncol(dust))
  mask  <- gas > 1e-22
  ratio[mask] <- dust[mask] / gas[mask]

  # User‑specified dust/gas levels (linear units)
  lev_lin <- c(4e-4, 0.00081, 0.0018, 0.0039, 0.008,
               0.014, 0.023, 0.036, 0.055, 0.08,
               0.11, 0.145, 0.19, 0.25, 0.5, 0.8, 1.0)
  min_lev <- min(lev_lin)

  # Log10 for plotting, clip very small positive values
  ratio_pos <- pmax(ratio, min_lev * 1e-2, na.rm = FALSE)
  z   <- log10(ratio_pos)
  lev <- log10(lev_lin)

  # Build grid
  nx <- nrow(z); ny <- ncol(z)
  x_seq <- seq(x_range[1], x_range[2], length.out = nx)
  y_seq <- seq(y_range[1], y_range[2], length.out = ny)

  # Inverted discrete palette (high ratio = light, low ratio = dark)
  nbin <- length(lev) - 1
  cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(nbin))

  png(out_file, width = 900, height = 800)
  par(oma = c(2.0, 3.0, 0.0, 0.5), mar = c(4, 5, 4, 1.5), cex.axis = 2)

  # Light‑grey background where gas = 0 (NA in z)
  plot.new()
  plot.window(xlim = x_range, ylim = y_range)
  rect(x_range[1], y_range[1], x_range[2], y_range[2],
       col = "lightgrey", border = NA)

  # Filled contours for ratio (only where z is finite)
  graphics::filled.contour(
    x = x_seq, y = y_seq, z = z,
    levels = lev,
    col    = cols,
    plot.title = {
      par(cex.main = 2.5, cex.lab = 2.5)
      title(main = "dust / gas ratio in E0 flyby plane",
            xlab = "km", ylab = "km")
    },
    plot.axes = {
      axis(1); axis(2)
      symbols(x = 0, y = 0, circles = 252, inches = FALSE, add = TRUE,
              bg = "darkgrey", fg = "darkgrey")
    },
    key.axes = {
      axis(4, at = lev, labels = lev_lin, cex.axis = 1.5)
    }
  )

  dev.off()
}
