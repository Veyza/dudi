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
  # Plot salt-rich fraction in the vertical structure plane

  library(RColorBrewer)

  if (close) {
    salt_poor_file <- "./results/salt_poor_0_plane_a-panel.dat"
    salt_rich_file <- "./results/salt_rich_0_plane_a-panel.dat"
    x_range <- c(-150, 150); y_range <- c(237, 372)
    out_file <- "./results/type3_proportion_a-panel.png"
  } else {
    salt_poor_file <- "./results/salt_poor_0_plane_b-panel.dat"
    salt_rich_file <- "./results/salt_rich_0_plane_b-panel.dat"
    x_range <- c(-450, 450); y_range <- c(0, 900)
    out_file <- "./results/type3_proportion_b-panel.png"
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

## 2) Densities in plane ------------------------------------------------------

density_in_plane <- function(quantity = c("number", "mass", "gas"),
                             close = FALSE) {
  # Plot number, mass, or gas density in the vertical structure plane

  library(RColorBrewer)
  quantity <- match.arg(quantity)

  if (close) {
    x_range <- c(-150, 150); y_range <- c(237, 372)
    if (quantity == "number") {
      salt_poor_file <- "./results/salt_poor_0_plane_a-panel.dat"
      salt_rich_file <- "./results/salt_rich_0_plane_a-panel.dat"
      out_file <- "./results/number_density_a-panel.png"
    } else if (quantity == "mass") {
      dens_file <- "./results/dens_in_0_plane_dust_a-panel.dat"
      out_file  <- "./results/mass_density_a-panel.png"
    } else { # gas
      dens_file <- "./results/dens_in_0_plane_gas_a-panel.dat"
      out_file  <- "./results/gas_density_a-panel.png"
    }
  } else {
    x_range <- c(-450, 450); y_range <- c(0, 900)
    if (quantity == "number") {
      salt_poor_file <- "./results/salt_poor_0_plane_b-panel.dat"
      salt_rich_file <- "./results/salt_rich_0_plane_b-panel.dat"
      out_file <- "./results/number_density_b-panel.png"
    } else if (quantity == "mass") {
      dens_file <- "./results/dens_in_0_plane_dust_b-panel.dat"
      out_file  <- "./results/mass_density_b-panel.png"
    } else { # gas
      dens_file <- "./results/dens_in_0_plane_gas_b-panel.dat"
      out_file  <- "./results/gas_density_b-panel.png"
    }
  }

  if (quantity == "number") {
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

    main_title <- "number density in vertical structure plane"
    key_lbls   <- lev_lin
  } else {
    s <- as.matrix(read.table(dens_file))

    s[s <= 0] <- NA
    z <- log10(s)
    rng <- range(z, finite = TRUE)

    lev <- pretty(rng, n = 12)
    key_lbls <- signif(10^lev, digits = 3)

    main_title <- if (quantity == "mass") then
      "mass density in vertical structure plane"
    else
      "gas density in vertical structure plane"
    end if
  }

  nbin <- length(lev) - 1
  cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, "YlGnBu"))(nbin))

  .plot_plane(
    z           = z,
    x_range     = x_range,
    y_range     = y_range,
    levels      = lev,
    cols        = cols,
    filename    = out_file,
    main        = main_title,
    key_labels  = key_lbls,
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
    dust_file <- "./results/dens_in_0_plane_dust_a-panel.dat"
    gas_file  <- "./results/dens_in_0_plane_gas_a-panel.dat"
    x_range <- c(-150, 150); y_range <- c(237, 372)
    out_file <- "./results/dust2gas_a-panel.png"
  } else {
    dust_file <- "./results/dens_in_0_plane_dust_b-panel.dat"
    gas_file  <- "./results/dens_in_0_plane_gas_b-panel.dat"
    x_range <- c(-450, 450); y_range <- c(0, 900)
    out_file <- "./results/dust2gas_b-panel.png"
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
               0.11, 0.145, 0.19, 0.25)
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
      title(main = "dust / gas ratio in vertical structure plane",
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
