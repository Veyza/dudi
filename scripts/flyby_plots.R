# Flyby plotting script for DUDI
# Flybys: 5, 17, 7 (7.1, 7.2, 7.3), 21 (21.1, 21.2, 21.3)
# Model data from ./results (flyby_profile output)
# CDA/HRD data from ./CDA_data
#
# Lines marked "# DATA: uncomment when CDA_data files are available" are commented out
# because the data files are not present. Uncomment them when you have the CDA/HRD files.

RESULTS_DIR <- "./results"
CDA_DIR <- "./CDA_data"


getprofileparams <- function(num) {
  if (num == 7.1) {
    fname1 <- file.path(CDA_DIR, "E7_M1.dat")
    fname2 <- file.path(CDA_DIR, "E7_mm1.dat")
    header1 <- "M1, R > 1.7"
    header2 <- "m1, R > 1.8"
    header3 <- "R > 1.75"
    timelinelabel <- "sec from 2 Nov 2009, 07:41:58"
    tmin <- -40
    tmax <- 40
    vel <- 7.7e3
    cospoint <- 0.96
  }
  if (num == 7.2) {
    fname1 <- file.path(CDA_DIR, "E7_M2.dat")
    fname2 <- file.path(CDA_DIR, "E7_mm2.dat")
    header1 <- "M2, R > 3.1"
    header2 <- "m2, R > 3.1"
    header3 <- "R > 3.1"
    timelinelabel <- "sec from 2 Nov 2009, 07:41:58"
    tmin <- -40
    tmax <- 40
    vel <- 7.7e3
    cospoint <- 0.96
  }
  if (num == 7.3) {
    fname1 <- file.path(CDA_DIR, "E7_M3.dat")
    fname2 <- NA
    header1 <- "M3, R > 6.7"
    header2 <- NA
    header3 <- "R > 6.7"
    timelinelabel <- "sec from 2 Nov 2009, 07:41:58"
    tmin <- -40
    tmax <- 40
    vel <- 7.7e3
    cospoint <- 0.96
  }
  if (num == 21.1) {
    fname1 <- file.path(CDA_DIR, "E21_M1.dat")
    fname2 <- file.path(CDA_DIR, "E21_mm1.dat")
    header1 <- "M1, R > 1.6"
    header2 <- "m1, R > 1.6"
    header3 <- "R > 1.6"
    timelinelabel <- "sec from 28 Oct 2015, 15:22:42"
    tmin <- -40
    tmax <- 40
    vel <- 8.5e3
    cospoint <- 0.975
  }
  if (num == 21.2) {
    fname1 <- file.path(CDA_DIR, "E21_M2.dat")
    fname2 <- file.path(CDA_DIR, "E21_mm2.dat")
    header1 <- "M2, R > 2.8"
    header2 <- "m2, R > 2.9"
    header3 <- "R > 2.85"
    timelinelabel <- "sec from 28 Oct 2015, 15:22:42"
    tmin <- -40
    tmax <- 40
    vel <- 8.5e3
    cospoint <- 0.975
  }
  if (num == 21.3) {
    fname1 <- file.path(CDA_DIR, "E21_M3.dat")
    fname2 <- NA
    header1 <- "M3, R > 6.2"
    header2 <- NA
    header3 <- "R > 6.2"
    timelinelabel <- "sec from 28 Oct 2015, 15:22:42"
    tmin <- -40
    tmax <- 40
    vel <- 8.5e3
    cospoint <- 0.975
  }

  data.frame(fn1 = fname1, fn2 = fname2, h1 = header1, h2 = header2, h3 = header3,
    camom = timelinelabel, tearly = tmin, tlate = tmax, speed = vel, pointing = cospoint)
}


densityplot <- function(num, data_only = FALSE) {
  if (!requireNamespace("latex2exp", quietly = TRUE)) {
    stop("latex2exp package required. Install with: install.packages('latex2exp')")
  }
  library(latex2exp)

  tmp <- num - round(num)
  profparams <- getprofileparams(num)
  # DATA: uncomment when CDA_data files are available
  # data1 <- read.table(profparams$fn1)
  # names(data1) <- c("time", "numdens")
  # errors1 <- sqrt(data1$numdens) / sqrt((data1$time[2] - data1$time[1]) *
  #   profparams$speed * 5e-3 * profparams$pointing)
  # data2 <- NULL
  # errors2 <- NULL
  # if (tmp < 0.29 && !is.na(profparams$fn2)) {
  #   data2 <- read.table(profparams$fn2)
  #   names(data2) <- c("time", "numdens")
  #   errors2 <- sqrt(data2$numdens) / sqrt((data2$time[2] - data2$time[1]) *
  #     profparams$speed * 1e-3 * profparams$pointing)
  # }

  # Model from flyby_profile output: ./results/E{num}_profile.dat
  modelfile <- file.path(RESULTS_DIR, paste0("E", as.character(num), "_profile.dat"))
  model <- NULL
  if (!data_only && file.exists(modelfile)) {
    model <- read.table(modelfile)
    names(model) <- c("time", "ndens", "frac1", "frac2", "frac3")
  }

  suff <- if (data_only) "dataonly" else ""
  fout <- file.path(RESULTS_DIR, paste0("E", as.character(round(num)),
    "M", as.character(round((num - round(num)) * 10)), "_comparison_3", suff, ".eps"))

  setEPS(width = 4.5, height = 4.5)
  postscript(fout)
  col1 <- "grey"
  par(mar = c(5, 5, 3, 3), lwd = 1)
  ymax <- if (!is.null(model)) max(model$ndens) * 1.2 else 1
  plot(1, type = "n",
    xlab = profparams$camom, ylab = latex2exp::TeX("number density, $m^{-3}$"),
    cex.axis = 1.3, yaxs = "i",
    xlim = c(profparams$tearly, profparams$tlate),
    ylim = c(0, ymax), cex.lab = 1.5, col = "white", pch = 18, cex = 1.5, cex.main = 1.5)
  # DATA: uncomment when CDA_data files are available
  # points(data1, cex = 1.5, col = col1, pch = 18)
  # arrows(data1$time, data1$numdens, data1$time, data1$numdens - errors1,
  #   length = 0.03, angle = 90, code = 3, col = col1)
  # arrows(data1$time, data1$numdens, data1$time, data1$numdens + errors1,
  #   length = 0.03, angle = 90, code = 3, col = col1)
  # if (!is.null(data2)) {
  #   col2 <- "cornsilk4"
  #   points(data2, cex = 1.5, col = col2, pch = 10)
  #   arrows(data2$time, data2$numdens, data2$time, data2$numdens - errors2,
  #     length = 0.03, angle = 90, code = 3, col = col2)
  #   arrows(data2$time, data2$numdens, data2$time, data2$numdens + errors2,
  #     length = 0.03, angle = 90, code = 3, col = col2)
  # }
  if (!data_only && !is.null(model)) {
    par(lwd = 3)
    lines(model$time, model$ndens, col = "purple", lwd = 3)
    legpos <- max(model$ndens) * 1.2
    legend(profparams$tearly - 3, legpos,
      legend = paste0("model ", profparams$h3),
      lwd = 3, col = "purple", lty = "solid", cex = 1.5, bty = "n")
  }
  # } else {
  #   legend(profparams$tearly - 3, max(data1$numdens) * 1.5,
  #     legend = paste("HRD, ", profparams$h1), pch = 18, col = col1, bty = "n", cex = 1.5)
  # }
  dev.off()
}


HRDplot3in1 <- function(num, data_only = FALSE) {
  if (!requireNamespace("latex2exp", quietly = TRUE)) {
    stop("latex2exp package required. Install with: install.packages('latex2exp')")
  }
  library(latex2exp)

  pchdata1 <- 10
  pchdata2 <- 18
  cexdata <- 2.1
  col1 <- "black"
  col2 <- "brown"
  tleg <- 3

  if (num == 7 || num == 21) {
    if (num == 7) {
      ticks1 <- c(0:5) * 3
      ticks2 <- c(0, 1, 2, 3)
      ticks3 <- c(0:3) * 0.1
    }
    if (num == 21) {
      ticks1 <- c(0:4) * 20
      ticks2 <- c(0:3) * 6
      ticks3 <- c(0:4) * 0.5
    }
    suff <- if (data_only) "dataonly" else ""
    fout <- file.path(RESULTS_DIR, paste0("HRD_E", as.character(num), "_data", suff, ".eps"))
    setEPS(width = 6, height = 9)
    postscript(fout)
    par(mfrow = c(3, 1), oma = c(10, 10, 2.0, 3.0))
    for (i in 1:3) {
      n <- num + i * 0.1
      profparams <- getprofileparams(n)
      # DATA: uncomment when CDA_data files are available
      # data1 <- read.table(profparams$fn1)
      # names(data1) <- c("time", "numdens")
      # errors1 <- sqrt(data1$numdens) / sqrt((data1$time[2] - data1$time[1]) *
      #   profparams$speed * 5e-3 * profparams$pointing)
      # data2 <- NULL
      # errors2 <- NULL
      # if (i < 3 && !is.na(profparams$fn2)) {
      #   data2 <- read.table(profparams$fn2)
      #   names(data2) <- c("time", "numdens")
      #   errors2 <- sqrt(data2$numdens) / sqrt((data2$time[2] - data2$time[1]) *
      #     profparams$speed * 1e-3 * profparams$pointing)
      # }
      modelfile <- file.path(RESULTS_DIR, paste0("E", as.character(n), "_profile.dat"))
      model <- NULL
      if (!data_only && file.exists(modelfile)) {
        model <- read.table(modelfile)
        names(model) <- c("time", "ndens", "frac1", "frac2", "frac3")
      }
      ymax <- if (!is.null(model)) max(model$ndens) * 1.2 else 1
      par(mar = c(0.0, 0.0, 2.6, 0.0))
      if (i < 3) {
        plot(1, type = "n",
          xlab = " ", ylab = " ", cex.axis = 1.8, yaxs = "i",
          xlim = c(profparams$tearly, profparams$tlate), xaxt = "n", yaxt = "n",
          ylim = c(0, ymax), cex.lab = 2, col = col1, pch = pchdata1,
          cex = cexdata, cex.main = 1.8)
        if (i < 2) {
          axis(side = 2, at = ticks1, pos = -43.1, lty = "solid", col = "black",
            las = 0, tck = -0.03, cex.axis = 2.8)
        } else {
          axis(side = 2, at = ticks2, pos = -43.1, lty = "solid", col = "black",
            las = 0, tck = -0.03, cex.axis = 2.8)
        }
      } else {
        plot(1, type = "n",
          xlab = " ", ylab = " ", cex.axis = 1.8, yaxs = "i",
          xlim = c(profparams$tearly, profparams$tlate), cex.axis = 2.5, xaxt = "n", yaxt = "n",
          ylim = c(0, ymax), cex.lab = 2, col = col1, pch = pchdata1,
          cex = cexdata, cex.main = 1.8)
        axis(side = 1, at = c(-40, -20, 0, 20, 40), pos = 0.0, lty = "solid",
          col = "black", las = 2, tck = -0.03, cex.axis = 2.8)
        axis(side = 2, at = ticks3, pos = -43.1, lty = "solid", col = "black",
          las = 0, tck = -0.03, cex.axis = 2.8)
      }
      # DATA: uncomment when CDA_data files are available
      # if (!is.null(data2)) {
      #   points(data2, cex = cexdata, col = col2, pch = pchdata2)
      #   arrows(data2$time, data2$numdens, data2$time, data2$numdens - errors2,
      #     length = 0.03, angle = 90, code = 3, col = col2)
      #   arrows(data2$time, data2$numdens, data2$time, data2$numdens + errors2,
      #     length = 0.03, angle = 90, code = 3, col = col2)
      # }
      # arrows(data1$time, data1$numdens, data1$time, data1$numdens - errors1,
      #   length = 0.03, angle = 90, code = 3, col = col1)
      # arrows(data1$time, data1$numdens, data1$time, data1$numdens + errors1,
      #   length = 0.03, angle = 90, code = 3, col = col1)
      if (!data_only && !is.null(model)) {
        lines(model$time, model$ndens, col = "purple", lwd = 3)
      }
      legend(tleg, ymax * 1.02,
        legend = latex2exp::TeX(paste(profparams$h1, "$\\mu m$")),
        pch = NA, col = NA, bty = "n", cex = 2.3)
    }
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "l", bty = "n", xaxt = "n", yaxt = "n")
    mtext(latex2exp::TeX("number density, $m^{-3} $"),
      side = 2, line = -5, outer = TRUE, cex = 2.5, las = 0)
    mtext(profparams$camom, side = 1, line = -3, outer = TRUE, cex = 2.5, las = 0)
    dev.off()
  }
}


readE5cdadata <- function() {
  # DATA: uncomment when CDA_data files are available
  # fname <- file.path(CDA_DIR, "E5_Datasheet_revision2021_V3_Anastasiia.txt")
  # d <- read.table(fname, header = TRUE, blank.lines.skip = TRUE, fill = TRUE)
  # d$Time <- as.character(d$Time)
  # hour <- numeric()
  # min <- numeric()
  # sec <- numeric()
  # for (i in seq_along(d$Time)) {
  #   hms <- unlist(strsplit(d$Time[i], ":"))
  #   hour[i] <- as.numeric(hms[1])
  #   min[i] <- as.numeric(hms[2])
  #   sec[i] <- as.numeric(hms[3])
  # }
  # d$hour <- hour
  # d$min <- min
  # d$sec <- sec + 3600 * hour + 60 * min
  # type4istype1 <- FALSE
  # type5istype3 <- FALSE
  # donotignoreS <- TRUE
  # camoment <- as.numeric(unlist(strsplit("19:06:40", ":")))
  # casec <- camoment[1] * 3600 + camoment[2] * 60 + camoment[3]
  # wind <- 9
  # tstep <- 3
  # tstart <- casec - 109
  # tfinish <- casec + 169
  # ttab <- numeric()
  # unknown <- numeric()
  # type1 <- numeric()
  # type2 <- numeric()
  # type3 <- numeric()
  # allspec <- numeric()
  # for (i in 1:((tfinish - tstart) / tstep)) {
  #   t <- tstart + (i - 1) * tstep
  #   ttab[i] <- t
  #   dtmp <- d[which(t - round(wind / 2 - 0.02) <= d$sec & d$sec <= t + round(wind / 2 - 0.02)), ]
  #   allspec[i] <- length(dtmp$Typ)
  #   if (type4istype1) {
  #     type1[i] <- length(which(dtmp$Typ == 1 | dtmp$Typ == 4 | dtmp$Typ == "4b"))
  #   } else {
  #     type1[i] <- length(which(dtmp$Typ == 1))
  #   }
  #   type2[i] <- length(which(dtmp$Typ == 2))
  #   if (type5istype3) {
  #     type3[i] <- length(which(dtmp$Typ == 3 | dtmp$Typ == 5))
  #   } else {
  #     type3[i] <- length(which(dtmp$Typ == 3))
  #   }
  #   unknown[i] <- length(which(dtmp$Typ == 0))
  #   if (donotignoreS) {
  #     ind <- which(dtmp$Typ == "S1" | dtmp$Typ == "S2" | dtmp$Typ == "S3" |
  #       dtmp$Typ == "s1" | dtmp$Typ == "s2" | dtmp$Typ == "s3")
  #     ind <- c(ind, which(dtmp$Typ == "S4" | dtmp$Typ == "S5" | dtmp$Typ == "s4" | dtmp$Typ == "s5"))
  #     if (length(ind) > 0) {
  #       s <- dtmp[ind, ]
  #       for (ii in 1:length(s$Typ)) {
  #         if (as.numeric(s$TypI[ii]) > 0) type1[i] <- type1[i] + as.numeric(s$TypI[ii]) / 100
  #         if ((!is.na(s$Typ4[ii]) & s$Typ4[ii] > 0) & type4istype1) {
  #           type1[i] <- type1[i] + as.numeric(s$Typ4[ii]) / 100
  #         }
  #         if (as.numeric(s$Typ2[ii]) > 0) type2[i] <- type2[i] + as.numeric(s$Typ2[ii]) / 100
  #         if (as.numeric(s$Typ3[ii]) > 0) type3[i] <- type3[i] + as.numeric(s$Typ3[ii]) / 100
  #         if (as.numeric(s$Typ3b[ii]) > 0 & type5istype3) {
  #           type3[i] <- type3[i] + as.numeric(s$Typ3b[ii]) / 100
  #         }
  #       }
  #     }
  #   } else {
  #     unknown[i] <- unknown[i] + length(which(dtmp$Typ == "S1" | dtmp$Typ == "S2" |
  #       dtmp$Typ == "S3" | dtmp$Typ == "s1" | dtmp$Typ == "s2" | dtmp$Typ == "s3"))
  #   }
  # }
  # data.frame(time = ttab - casec, t1 = type1, t2 = type2, t3 = type3,
  #   und = unknown, all = allspec)
  return(NULL)
}


readE17cdadata <- function(fname1) {
  # DATA: uncomment when CDA_data files are available
  # win <- 60
  # step <- 15
  # m <- 20
  # cdadata <- read.table(fname1, header = TRUE)
  # cdadata <- cdadata[9:length(cdadata$time), ]
  # t <- double()
  # type1 <- double()
  # type2 <- double()
  # type3 <- double()
  # ind <- integer()
  # unidentified <- double()
  # N <- double()
  # for (i in 1:m) {
  #   ind1 <- which(cdadata$time >= (cdadata$time[1] + (i - 1) * step))
  #   ind2 <- which(cdadata$time < (cdadata$time[1] + (i - 1) * step + win))
  #   if (min(ind1) < max(ind2)) {
  #     ind <- min(ind1):max(ind2)
  #     set <- data.frame(time = cdadata$time[ind], t1 = cdadata$t1[ind],
  #       t2 = cdadata$t2[ind], t3 = cdadata$t3[ind], und = cdadata$unidentified[ind])
  #     t[i] <- sum(set$time) / length(set$time)
  #     type1[i] <- sum(set$t1)
  #     type2[i] <- sum(set$t2)
  #     type3[i] <- sum(set$t3)
  #     unidentified[i] <- sum(set$und)
  #     N[i] <- length(set$time)
  #   }
  # }
  # data.frame(time = t, t1 = type1, t2 = type2, t3 = type3, all = N, und = unidentified)
  return(NULL)
}


flybyplot <- function(num, data_only = FALSE) {
  if (!requireNamespace("latex2exp", quietly = TRUE)) {
    stop("latex2exp package required. Install with: install.packages('latex2exp')")
  }
  library(latex2exp)

  # Model files: k=1 and k=2 (E5) or k=1.7 (E17)
  relation <- NULL
  relation2 <- NULL
  k1_label <- "k = 1"
  k2_label <- NULL
  if (num == 5) {
    modelfile1 <- file.path(RESULTS_DIR, "E5_profile.dat")
    modelfile2 <- file.path(RESULTS_DIR, "E5.2_profile.dat")
    k2_label <- "k = 2"
  } else if (num == 17) {
    modelfile1 <- file.path(RESULTS_DIR, "E17_profile.dat")
    modelfile2 <- file.path(RESULTS_DIR, "E17.17_profile.dat")
    k2_label <- "k = 1.7"
  } else {
    modelfile1 <- file.path(RESULTS_DIR, paste0("E", as.character(num), "_profile.dat"))
    modelfile2 <- NULL
  }

  if (!data_only && (num == 5 || num == 17)) {
    if (file.exists(modelfile1)) {
      tmptab <- read.table(modelfile1)
      relation <- data.frame(time = tmptab$V1, t1 = tmptab$V3, t2 = tmptab$V4, t3 = tmptab$V5)
    }
    if (!is.null(modelfile2) && file.exists(modelfile2)) {
      tmptab2 <- read.table(modelfile2)
      relation2 <- data.frame(time = tmptab2$V1, t1 = tmptab2$V3, t2 = tmptab2$V4, t3 = tmptab2$V5)
    }
  } else if (!data_only && (num > 7 && num < 8 || num > 21 && num < 22)) {
    modelfile1 <- file.path(RESULTS_DIR, paste0("E", as.character(num), "_profile.dat"))
    if (file.exists(modelfile1)) {
      tmptab <- read.table(modelfile1)
      relation <- data.frame(time = tmptab$V1, t1 = tmptab$V3, t2 = tmptab$V4, t3 = tmptab$V5)
    }
  }

  if (num == 5) {
    timelinelabel <- "sec from 9 Oct 2008, 19:06:40"
  }
  if (num == 17) {
    timelinelabel <- "sec from 27 Mar 2012, 18:30:09"
  }

  # Density plots for E7.x and E21.x
  if ((num > 7 && num < 8) || (num > 21 && num < 22)) {
    densityplot(num, data_only = data_only)
  }

  # Fraction plots for E5 and E17
  if (num == 17) {
    # DATA: uncomment when CDA_data files are available
    # cdadata <- readE17cdadata(file.path(CDA_DIR, "CDA_E17_data_New.dat"))
    # all_dust <- cdadata$t1 + cdadata$t2 + cdadata$t3
    # sterr1 <- sqrt(cdadata$t1 / all_dust / cdadata$all * (1.0 - cdadata$t1 / all_dust)) * 100
    # sterr2 <- sqrt(cdadata$t2 / all_dust / cdadata$all * (1.0 - cdadata$t2 / all_dust)) * 100
    # sterr3 <- sqrt(cdadata$t3 / all_dust / cdadata$all * (1.0 - cdadata$t3 / all_dust)) * 100
    # err1 <- sqrt((cdadata$und / (cdadata$all) * 100 / 3)^2 + sterr1^2)
    # err2 <- sqrt((cdadata$und / (cdadata$all) * 100 / 3)^2 + sterr2^2)
    # err3 <- sqrt((cdadata$und / (cdadata$all) * 100 / 3)^2 + sterr3^2)
    suff <- if (data_only) "dataonly" else ""
    fout <- file.path(RESULTS_DIR, paste0("E17_fit", suff, ".eps"))
    setEPS(width = 6, height = 6)
    postscript(fout)
    plot(1, type = "n", xlab = timelinelabel,
      ylab = "fraction of dust types, %", cex.axis = 1.3, cex.lab = 1.5,
      xlim = c(-61, 65), ylim = c(0, 100), xaxs = "i", yaxs = "i")
    if (!data_only && !is.null(relation)) {
      lines(relation$time, relation$t1 * 100, col = "black", lwd = 3, lty = "longdash")
      lines(relation$time, relation$t2 * 100, col = "blue", lwd = 3, lty = "longdash")
      lines(relation$time, relation$t3 * 100, col = "red", lwd = 3, lty = "longdash")
    }
    if (!data_only && !is.null(relation2)) {
      lines(relation2$time, relation2$t1 * 100, col = "black", lwd = 3, lty = "solid")
      lines(relation2$time, relation2$t2 * 100, col = "blue", lwd = 3, lty = "solid")
      lines(relation2$time, relation2$t3 * 100, col = "red", lwd = 3, lty = "solid")
    }
    # DATA: uncomment when CDA_data files are available
    # arrows(cdadata$time, cdadata$t1 / all_dust * 100 - err1, cdadata$time,
    #   cdadata$t1 / all_dust * 100 + err1, length = 0.03, angle = 90, code = 3)
    # arrows(cdadata$time, cdadata$t2 / all_dust * 100 - err2, cdadata$time,
    #   cdadata$t2 / all_dust * 100 + err2, length = 0.03, angle = 90, code = 3, col = "blue")
    # arrows(cdadata$time, cdadata$t3 / all_dust * 100 - err3, cdadata$time,
    #   cdadata$t3 / all_dust * 100 + err3, length = 0.03, angle = 90, code = 3, col = "red")
    # points(cdadata$time, cdadata$t1 / all_dust * 100, pch = 2, col = "black", cex = 1.6)
    # points(cdadata$time, cdadata$t2 / all_dust * 100, pch = 1, col = "blue", cex = 1.6)
    # points(cdadata$time, cdadata$t3 / all_dust * 100, pch = 0, col = "red", cex = 1.6)
    # lines(cdadata$time, cdadata$t1 / all_dust * 100, col = "black", lty = "dashed")
    # lines(cdadata$time, cdadata$t2 / all_dust * 100, col = "blue", lty = "dashed")
    # lines(cdadata$time, cdadata$t3 / all_dust * 100, col = "red", lty = "dashed")
    legend("topleft", legend = c("type I", "type II", "type III"),
      col = c("black", "blue", "red"), lty = c(5, 5, 5), lwd = 3,
      text.col = c("black", "blue", "red"), cex = 1.3, bty = "n")
    if (!is.null(relation2)) {
      legend("topleft", legend = c(k1_label, k2_label),
        col = c("gray40", "gray40"), lty = c(5, 1), lwd = 3,
        text.col = c("gray40", "gray40"), cex = 1.3, inset = c(0.28, 0), bty = "n")
    } else {
      legend("topleft", legend = k1_label,
        col = "gray40", lty = 5, lwd = 3, text.col = "gray40",
        cex = 1.3, inset = c(0.28, 0), bty = "n")
    }
    dev.off()
  }

  if (num == 5) {
    # DATA: uncomment when CDA_data files are available
    # cdadata <- readE5cdadata()
    # all_dust <- cdadata$t1 + cdadata$t2 + cdadata$t3
    # sterr1 <- sqrt(cdadata$t1 / all_dust / cdadata$all * (1.0 - cdadata$t1 / all_dust)) * 100
    # sterr2 <- sqrt(cdadata$t2 / all_dust / cdadata$all * (1.0 - cdadata$t2 / all_dust)) * 100
    # sterr3 <- sqrt(cdadata$t3 / all_dust / cdadata$all * (1.0 - cdadata$t3 / all_dust)) * 100
    # err1 <- sqrt((cdadata$und / (cdadata$all) * 100 / 3)^2 + sterr1^2)
    # err2 <- sqrt((cdadata$und / (cdadata$all) * 100 / 3)^2 + sterr2^2)
    # err3 <- sqrt((cdadata$und / (cdadata$all) * 100 / 3)^2 + sterr3^2)
    # ind1 <- which(cdadata$time > -65 & cdadata$time < 17)
    # ind2 <- which(cdadata$time > 35)
    # ind <- c(ind1, ind2)
    suff <- if (data_only) "dataonly" else ""
    fout <- file.path(RESULTS_DIR, paste0("E5_fit", suff, ".eps"))
    setEPS(width = 6, height = 6)
    postscript(fout)
    plot(1, type = "n", xlab = timelinelabel,
      ylab = "fraction of dust types, %", cex.axis = 1.3, cex.lab = 1.5,
      xlim = c(-60, 60), ylim = c(0, 100), xaxs = "i", yaxs = "i")
    rect(18, -10, 35, 110, density = 25, col = "darkgrey")
    if (!data_only && !is.null(relation)) {
      lines(relation$time, relation$t1 * 100, col = "black", lwd = 3, lty = "longdash")
      lines(relation$time, relation$t2 * 100, col = "blue", lwd = 3, lty = "longdash")
      lines(relation$time, relation$t3 * 100, col = "red", lwd = 3, lty = "longdash")
    }
    if (!data_only && !is.null(relation2)) {
      lines(relation2$time, relation2$t1 * 100, col = "black", lwd = 3, lty = "solid")
      lines(relation2$time, relation2$t2 * 100, col = "blue", lwd = 3, lty = "solid")
      lines(relation2$time, relation2$t3 * 100, col = "red", lwd = 3, lty = "solid")
    }
    # DATA: uncomment when CDA_data files are available
    # arrows(cdadata$time[ind], cdadata$t1[ind] / all_dust[ind] * 100 - err1[ind],
    #   cdadata$time[ind], cdadata$t1[ind] / all_dust[ind] * 100 + err1[ind],
    #   length = 0.03, angle = 90, code = 3)
    # arrows(cdadata$time[ind], cdadata$t2[ind] / all_dust[ind] * 100 - err2[ind],
    #   cdadata$time[ind], cdadata$t2[ind] / all_dust[ind] * 100 + err2[ind],
    #   length = 0.03, angle = 90, code = 3, col = "blue")
    # arrows(cdadata$time[ind], cdadata$t3[ind] / all_dust[ind] * 100 - err3[ind],
    #   cdadata$time[ind], cdadata$t3[ind] / all_dust[ind] * 100 + err3[ind],
    #   length = 0.03, angle = 90, code = 3, col = "red")
    # points(cdadata$time[ind], cdadata$t1[ind] / all_dust[ind] * 100,
    #   pch = 2, col = "black", cex = 1.6)
    # points(cdadata$time[ind], cdadata$t2[ind] / all_dust[ind] * 100,
    #   pch = 1, col = "blue", cex = 1.6)
    # points(cdadata$time[ind], cdadata$t3[ind] / all_dust[ind] * 100,
    #   pch = 0, col = "red", cex = 1.6)
    # lines(cdadata$time[ind1], cdadata$t1[ind1] / all_dust[ind1] * 100,
    #   col = "black", lty = "dashed")
    # lines(cdadata$time[ind1], cdadata$t2[ind1] / all_dust[ind1] * 100,
    #   col = "blue", lty = "dashed")
    # lines(cdadata$time[ind1], cdadata$t3[ind1] / all_dust[ind1] * 100,
    #   col = "red", lty = "dashed")
    # lines(cdadata$time[ind2], cdadata$t1[ind2] / all_dust[ind2] * 100,
    #   col = "black", lty = "dashed")
    # lines(cdadata$time[ind2], cdadata$t2[ind2] / all_dust[ind2] * 100,
    #   col = "blue", lty = "dashed")
    # lines(cdadata$time[ind2], cdadata$t3[ind2] / all_dust[ind2] * 100,
    #   col = "red", lty = "dashed")
    legend("topleft", legend = c("type I", "type II", "type III"),
      col = c("black", "blue", "red"), lty = c(5, 5, 5), lwd = 3,
      text.col = c("black", "blue", "red"), cex = 1.3, bty = "n")
    if (!is.null(relation2)) {
      legend("topleft", legend = c(k1_label, k2_label),
        col = c("gray40", "gray40"), lty = c(5, 1), lwd = 3,
        text.col = c("gray40", "gray40"), cex = 1.3, inset = c(0.28, 0), bty = "n")
    } else {
      legend("topleft", legend = k1_label,
        col = "gray40", lty = 5, lwd = 3, text.col = "gray40",
        cex = 1.3, inset = c(0.28, 0), bty = "n")
    }
    dev.off()
  }
}


numdensprofiles <- function(data_only = FALSE) {

  # HRD 3-in-1 for E7 and E21
  HRDplot3in1(7, data_only = data_only)
  HRDplot3in1(21, data_only = data_only)

}
