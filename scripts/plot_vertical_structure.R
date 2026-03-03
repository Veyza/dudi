composition_in_plane <- function(fnum, close = F){
	library(RColorBrewer)
	library(fields)
	library("viridis")
	library(latex2exp)
	library(plotrix)
	fnames <- character(2)
	if(close){
		fnames[1] <- "./results/salt_poor_0_plane_close.dat"
        fnames[2] <- "./results/salt_rich_0_plane_close.dat"
	}else{
		fnames[1] <- "./results/salt_poor_0_plane.dat"
        fnames[2] <- "./results/salt_rich_0_plane.dat"
	}

	t1 <- as.matrix(read.table(fnames[1]))
	t3 <- as.matrix(read.table(fnames[2]))
	r3 <- matrix(0, nrow = nrow(t1), ncol = ncol(t1))


    nn <- length(t1[1,])

	s <- t1 + t3

	for(i in seq_along(t1[1,])){
		for(ii in seq_along(t3[,1])){
			if(s[ii,i] > 1e-22){
				r3[ii,i] <- t3[ii,i] / s[ii,i] * 100
			}else{
				r3[ii,i] <- 0.0
			}
		}
	}
	print(max(r3))
	if(close){
		r3 <- r3[,15:150]
	}

	magma <- viridis::magma(256)
	mako  <- viridis::mako(256)

	#lev <- sort(c(2e-8, 1e-8, 8e-9, 6.5e-9, 5e-9, 3e-9, 2e-9, 1e-9, 5e-10, 1e-10, 1e-11, 1e-12, 1e-13)) * 1e-1
	deg = 1.0 / 3.0

	if(fnum == 0){
		if(close){
			lev2 <- c(0, 10, 20, 25, 27, 30, 33, 36, 40, 45, 50, 55, 60, 65, 75)
			lev3 <- c(5, 10, 17, c(3:9)*10)
			indlev <- c(1, 2, 3, 4, 5, 7)
		}else{
			lev2 <- c(0, 10, 20, 25, 27, 30, 33, 36, 40, 45, 50, 55, 60, 65, 75)
			lev3 <- c(1:14)*5
			indlev <- c(1, 2, 3, 5, 7, 9, 11)
		}

	}
	if(fnum == 0){
		if(close){
			x <- c(-150, 150)
			y <- c(237,372)
		}else{
			x <- c(-450, 450)
			y <- c(0,900)
		}


	#png(paste0("type3_proportion_E", paste0(as.character(fnum), ".png")), width = 900, height = 800)
	if(close){
		png("./results/type3_proportion_close.png", width = 900, height = 600)
	}else{
		png("./results/type3_proportion.png", width = 900, height = 800)
	}

		#mtitle <- paste0("type III percentage in E", paste0(as.character(fnum), " flyby plane"))
		mtitle = " "
		par(oma = c(2.0, 3.0, 2.0, 3), mar=c(4, 5, 6, 3))
		par(cex.axis = 2)
		graphics::filled.contour(x = seq(x[1], x[2], length.out = nrow(r3)), y = seq(from = y[1], to = y[2], length.out = ncol(r3)), r3, col = magma, levels = lev3,
	               plot.title = {par(cex.main=2.5, cex.lab = 2.5); title(main = mtitle,
	                                  xlab = "", ylab = " ")},
	               key.title = {par(cex.main=2)},
	               plot.axes = {axis(1); axis(2);
	               	#draw.circle(-1.5, y[2]*0.374, radius = 252, col = "lightgrey", border = "darkgrey");
	                if(fnum == 5 | fnum == 17){
	                	lines(xfl, yfl, lty = "longdash", lwd = 3, col = "grey");
	                	for(i in seq_along(tmarks)){
		            		points(tticks, pch = 3, cex = 3, col = "grey")
		            		if(fnum == 5){
		            			text(tticks[i,1]-40, tticks[i,2]-40, tmarks[i], cex = 2.5, col = "grey")
		            		}
		            		if(fnum == 17){
		            			text(tticks[i,1], tticks[i,2]-35, tmarks[i], cex = 2.5, col = "grey")
		            		}
	            	}
	                }
	            contour(x = seq(x[1], x[2], length.out = nrow(r3)), y = seq(from = y[1], to = y[2], length.out = ncol(r3)),
			               z = r3,
			               zlim = range(r3, finite = TRUE),
			               levels = lev3[indlev],
			               labels = paste(c(as.character(lev3[indlev])), "%"),
			               drawlabels = ! close,
			               method = "flattest",
			               lwd = 2, col = "white", labcex = 1.8,
							add = T);
	            if(close){
	            	draw.circle(1, 120, radius = 192.5, col = "azure4", border = "azure4");
	            	text(59, 362, paste(as.character(lev3[1]), "%"), col = "white", cex = 1.8)
					text(-117, 310, paste(as.character(lev3[2]), "%"), col = "white", cex = 1.8)
					text(125, 260, paste(as.character(lev3[3]), "%"), col = "white", cex = 1.8)
            	}else{
            		draw.circle(-1.5, -12, radius = 257, col = "azure4", border = "azure4");
            	}

	            }
	               #key.axes = axis(side = 4, cex.axis = 2, at = lev, labels = as.character(lev^(1.0/deg)))
	               #key.axes = axis(side = 4, cex.axis = 2, at = lev, labels = as.character(lev))

	         	)
		mtext("km from equatorial plane", side = 2, outer = TRUE, line = 0, cex = 2.5)
		mtext("km from polar axis projection         ", side = 1, outer = TRUE, line = 0, cex = 2.5)
		mtext("salt-rich dust proportion    ", side = 3, outer = TRUE, line = -1.5, cex = 3.5)
		mtext("Alexandria                       Damascus        ", side = 3, outer = TRUE, line = -4.7, cex = 2)
		#mtext("                                                                                                             %", side = 3, outer = T, line = -5.5, cex = 1.8)
		if(close){
			arrows(-45, 375, -150, 375, xpd = TRUE, angle = 15)
			arrows(5, 375, 100, 375, xpd = TRUE, angle = 15)
		}else{
			arrows(-150, 915, -370, 915, xpd = TRUE, angle = 15)
			arrows(0, 915, 220, 915, xpd = TRUE, angle = 15)
		}

	dev.off()

	}

	lev0 <- c(0.1, 0.3, 1, 10, 100, 200, 300, 500, 1000, 2000, 3000, 4000)
	png(paste0("number_density_E", paste0(as.character(fnum), ".png")), width = 900, height = 800)
		par(oma = c(2.0, 3.0, 0.0, 0.5), mar=c(4, 5, 4, 1.5))
		par(cex.axis = 2)
		graphics::filled.contour(x = seq(from = x[1], to = x[2], length.out = nrow(s)), y = seq(from = y[1], to = y[2], length.out = ncol(s)), s, col = mako, levels = lev0,
	               plot.title = {par(cex.main=2.5, cex.lab = 2.5);title(main = paste0("number density in E", paste0(as.character(fnum), " flyby plane")),
	                                  xlab = "km", ylab = " ")},
	               key.title = {par(cex.main=2)},
	               plot.axes = {axis(1); axis(2); cex.axes = 2;
	               	#draw.circle(-1.5, y[2]*0.374, radius = 252, col = "lightgrey", border = "grey");
	               	draw.circle(-1.5, -0, radius = 252, col = "lightgrey", border = "darkgrey");
	               	if(fnum == 5 | fnum == 17){
	                	lines(xfl, yfl, lty = "longdash", lwd = 3, col = "grey");
	                	for(i in seq_along(tmarks)){
		            		points(tticks, pch = 3, cex = 3, col = "grey")
		            		if(fnum == 5){
		            			text(tticks[i,1]-40, tticks[i,2]-40, tmarks[i], cex = 2.5, col = "grey")
		            		}
		            		if(fnum == 17){
		            			text(tticks[i,1], tticks[i,2]-35, tmarks[i], cex = 2.5, col = "grey")
		            		}
	            	}
	                }
	               }
	               #key.axes = axis(side = 4, cex.axis = 2, at = lev, labels = as.character(lev^(1.0/deg)))
	               #key.axes = axis(side = 4, cex.axis = 2, at = lev, labels = as.character(lev))

	         	)
		mtext("km", side = 2, outer = TRUE, line = 0, cex = 2.5)

	dev.off()
}
