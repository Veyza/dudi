plot_horizontal_structure <- function(
  altitudes_km = c(25, 50, 75, 100),
  results_dir = "./results",
  input_dir   = "./input_data_files"
) {
  # Required packages
  if (!requireNamespace("akima", quietly = TRUE)) {
    stop("Please install package 'akima' (install.packages('akima')).")
  }
  if (!requireNamespace("fields", quietly = TRUE)) {
    stop("Please install package 'fields' (install.packages('fields')).")
  }
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Please install package 'RColorBrewer' (install.packages('RColorBrewer')).")
  }

  library(akima)
  library(fields)
  library(RColorBrewer)

  # ------------------------------------------------------------------
  # Helper: read one altitude file, return data.frame(lat, lon, dens_g)
  # ------------------------------------------------------------------
  read_surface_file <- function(alt_km) {
    fname <- file.path(results_dir,
                       sprintf("surface_dust_mass_density_%dkm.dat", alt_km))
    dat <- read.table(fname)
    if (ncol(dat) < 3) {
      stop("File ", fname, " must have at least 3 columns: lat lon density")
    }
    colnames(dat)[1:3] <- c("lat_deg", "lon_deg", "dens_kg_m3")
    dat$dens_g_m3 <- dat$dens_kg_m3 * 1e3  # kg/m^3 -> g/m^3
    dat
  }

  # ------------------------------------------------------------------
  # Helper: read jet and diffuse source positions for overlay
  # vertical_jets.dat and diffuse_sources.dat:
  #   header line (ignored), then alphaM betaM zeta eta (degrees)
  # We treat (alphaM, betaM) as latitude, longitude in degrees.
  # ------------------------------------------------------------------
  read_sources <- function() {
    read_src <- function(fname) {
      f <- file.path(input_dir, fname)
      con <- file(f, "r")
      on.exit(close(con))
      header <- readLines(con, n = 1)
      dat <- read.table(con)
      if (ncol(dat) < 2) stop("File ", f, " must have at least 2 columns.")
      colnames(dat)[1:2] <- c("alpha_deg", "beta_deg")
      # interpret as latitude/longitude in degrees
      data.frame(
        lat_deg = dat$alpha_deg,
        lon_deg = dat$beta_deg
      )
    }
    jets    <- read_src("vertical_jets.dat")
    diffuse <- read_src("diffuse_sources.dat")
    list(jets = jets, diffuse = diffuse)
  }

  src <- read_sources()

  # ------------------------------------------------------------------
  # Projection from latitude/longitude to local Cartesian coordinates
  # centred at the south pole. We neglect curvature and simply use the
  # Cartesian components on a sphere of radius R_km, then drop the z
  # coordinate: south pole (lat = -90 deg) maps to (0, 0).
  # ------------------------------------------------------------------
  latlon_to_xy <- function(lat_deg, lon_deg, R_km) {
    lat_rad <- lat_deg * pi / 180
    lon_rad <- lon_deg * pi / 180
    x <- R_km * cos(lat_rad) * cos(lon_rad)
    y <- R_km * cos(lat_rad) * sin(lon_rad)
    list(x = x, y = y)
  }

  # Read all altitude slices and attach projected x/y (km)
  slices <- vector("list", length(altitudes_km))
  names(slices) <- paste0(altitudes_km, " km")
  for (k in seq_along(altitudes_km)) {
    alt_km <- altitudes_km[k]
    dat <- read_surface_file(alt_km)
    R_km <- 252 #+ alt_km
    xy <- latlon_to_xy(dat$lat_deg, dat$lon_deg, R_km)
    dat$x_plot <-  xy$y
    dat$y_plot <-  xy$x
    slices[[k]] <- dat
  }

  # Project jets and diffuse sources to surface plane (radius ~ 252 km)
  R_surf_km <- 252
  jets_xy <- latlon_to_xy(src$jets$lat_deg, src$jets$lon_deg, R_surf_km)
  diffuse_xy <- latlon_to_xy(src$diffuse$lat_deg, src$diffuse$lon_deg, R_surf_km)
  src$jets$x_plot    <-  jets_xy$y
  src$jets$y_plot    <-  jets_xy$x
  src$diffuse$x_plot <-  diffuse_xy$y
  src$diffuse$y_plot <-  diffuse_xy$x

  # Common plotting ranges in km (include sources and all altitudes)
  all_x <- range(c(
    unlist(lapply(slices, function(d) d$x_plot)),
    src$jets$x_plot,
    src$diffuse$x_plot
  ), na.rm = TRUE)
  all_y <- range(c(
    unlist(lapply(slices, function(d) d$y_plot)),
    src$jets$y_plot,
    src$diffuse$y_plot
  ), na.rm = TRUE)

  # ------------------------------------------------------------------
  # Colour palette and level selection
  # ------------------------------------------------------------------
  pal <- colorRampPalette(rev(brewer.pal(9, "BuGn")))  # dark -> bright green
  ncols <- 64

  # Output PNG
  png(file.path(results_dir, "dust_mass_density_at_altitudes.png"),
      width = 2200, height = 1400)
  on.exit(dev.off(), add = TRUE)

  # Layout: 2x2 panels
  par(bg = "black",
      mfrow = c(2, 2),
      mar   = c(2, 2, 3, 18),   # bottom, left, top, **right** (bigger right = narrower plot, wider legend)
      oma   = c(0, 2, 0, 2) ,   # outer margins; moderate space on right for all legends together
      cex.axis = 1.2 * 4,
      cex.lab = 1.4 * 4,
      cex.main = 1.6 * 4)

  # Titles per altitude (top of each panel)
  make_title <- function(alt_km) {
    sprintf("dust density at altitude %d km", alt_km)
  }

  # Choose z scale per panel (0 .. slightly above max)
  get_z_limits <- function(dens_g) {
    zmax <- max(dens_g, na.rm = TRUE)
    c(0, zmax)
  }

  # ------------------------------------------------------------------
  # Draw orbit and Saturn arrows (approximate, using x/y ranges)
  # ------------------------------------------------------------------
  draw_annotations <- function(xlim, ylim) {
    xmid <- mean(xlim)
    ymid <- mean(ylim)

    # Short Orbit arrow along right edge, pointing right
    y_orb <- ymid
    x0_orb <- xlim[2] - 0.1 * diff(xlim)
    x1_orb <- xlim[2] - 0.01 * diff(xlim)
    arrows(x0 = x0_orb, y0 = y_orb,
           x1 = x1_orb, y1 = y_orb,
           length = 0.3, col = "white", lwd = 2.2)
    text(x0_orb, y_orb,
         labels = "Orbit", col = "white", pos = 2, cex = 3.5)

    # Short Saturn arrow near top edge
    x_s <- xmid
    y_s0 <- ylim[2] - 0.1 * diff(ylim)
    y_s1 <- ylim[2] - 0.01 * diff(ylim)
    arrows(x0 = x_s, y0 = y_s0, x1 = x_s, y1 = y_s1,
           length = 0.3, col = "white", lwd = 2.2)
    text(x_s, y_s0,
         labels = "Saturn", col = "white", pos = 1, cex = 3.5)
  }

  # ------------------------------------------------------------------
  # Main loop over altitudes: filled-contour maps of mass density
  # in g/m^3, south pole at (0, 0) in projected x/y.
  # ------------------------------------------------------------------
  for (k in seq_along(altitudes_km)) {
    alt_km <- altitudes_km[k]
    dat <- slices[[k]]
    z <- dat$dens_g_m3

    # Interpolate scattered points to a regular grid in (x, y)
    nx <- 200
    ny <- 200
    x_seq <- seq(all_x[1], all_x[2], length.out = nx)
    y_seq <- seq(all_y[1], all_y[2], length.out = ny)

    itp <- akima::interp(
      x = dat$x_plot,
      y = dat$y_plot,
      z = z,
      xo = x_seq,
      yo = y_seq,
      duplicate = "mean",
      extrap = FALSE
    )

    zlim <- get_z_limits(itp$z)
    breaks <- seq(zlim[1], zlim[2], length.out = ncols + 1)
    cols <- pal(ncols)

    # Filled-contour style image with colour bar (with numeric labels)
    fields::image.plot(
      x = itp$x,
      y = itp$y,
      z = itp$z,
      xlim = all_x,
      ylim = all_y,
      zlim = zlim,
      breaks = breaks,
      col = cols,
      xlab = "",
      ylab = "",
      main = make_title(alt_km),
      axes = FALSE,
      asp = 1,
      legend.mar = 0.5,
      legend.args = list(
        text = expression(g/m^3),
        col  = "white",
        cex  = 2.7
      ),
      axis.args = list(
        col       = "white",
        col.axis  = "white",
        cex.axis  = 3.2
      )
    )

    box(col = "white")
    # No axis ticks or labels; keep only frame

    # Overlay jets (white diamonds) and diffuse sources (white open circles)
    points(src$jets$x_plot, src$jets$y_plot,
           pch = 17, col = "yellow", cex = 1.5, lwd = 2)  # filled triangles
    points(src$diffuse$x_plot, src$diffuse$y_plot,
           pch = 21, col = "yellow", bg = NA, cex = 1.7, lwd = 2)  # empty circles

    # Legend for jets / diffuse sources
    legend("bottomleft",
           legend = c("jets", "diffuse sources"),
           pch = c(17, 21),
           pt.bg = c(NA, NA),
           col = "yellow",
           bty = "n",
           text.col = "white",
           cex = 0.9 * 4)

    # Orbit/Saturn annotations
    draw_annotations(all_x, all_y)
  }

  invisible(NULL)
}
