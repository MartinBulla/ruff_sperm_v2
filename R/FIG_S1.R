# TOOLS
  require(here)
  source(here::here('R/tools.R'))
  require(GGally) 
  require(gtools)
  require(PerformanceAnalytics)
  require(psych)

  # function (from https://stackoverflow.com/questions/56728644/how-to-remove-significance-stars-from-chart-correlation-function)
  chart.Correlation.nostars <- function (R, histogram = TRUE, method = c("pearson", "kendall", "spearman"), ...) 
    {
      x = checkData(R, method = "matrix")
      if (missing(method)) 
        method = method[1]
      panel.cor <- function(x, y, digits = 2, prefix = "", 
                            use = "pairwise.complete.obs", method = "pearson", 
                            cex.cor, ...) {
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        r <- cor(x, y, use = use, method = method)
        txt <- format(c(r, 0.123456789), digits = digits)[1]
        txt <- paste(prefix, txt, sep = "")
        if (missing(cex.cor)) 
          cex <- 0.8/strwidth(txt)
        test <- cor.test(as.numeric(x), as.numeric(y), method = method)
        # Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
        #                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
        #                                                                           "**", "*", ".", " "))
        text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
        # text(0.8, 0.8, Signif, cex = cex, col = 2)
      }
      f <- function(t) {
        dnorm(t, mean = mean(x), sd = sd.xts(x))
      }
      dotargs <- list(...)
      dotargs$method <- NULL
      rm(method)
      hist.panel = function(x, ... = NULL) {
        par(new = TRUE)
        hist(x, col = "light gray", probability = TRUE, 
             axes = FALSE, main = "", breaks = "FD")
        lines(density(x, na.rm = TRUE), col = "red", lwd = 1)
        rug(x)
      }
      if (histogram) 
        pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
              diag.panel = hist.panel)
      else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor)
    }
  # modified pairs_panels - adjusting coefficient size, line colors and using loess also when no CIs are drawn
   pairs.panels_MB<- function (x, smooth = TRUE, scale = TRUE, density = TRUE, ellipses = FALSE, 
    digits = 2, method = "pearson", pch = 20, lm = FALSE, cor = TRUE, 
    jiggle = FALSE, factor = 2, hist.col = "light gray",dens.col = 'red',smooth.col ='red', show.points = TRUE, 
    rug = TRUE, breaks = "Sturges", cex.cor = 1, wt = NULL, smoother = FALSE, 
    stars = FALSE, ci = FALSE, alpha = 0.05, ...) 
    {
    "panel.hist.density" <- function(x, ...) {
        usr <- par("usr")
        on.exit(par("usr"))
        par(usr = c(usr[1], usr[2], 0, 1.5), tck = -0.03, mgp=c(2,0.4,0), las = 1)
        tax <- table(x)
        if (length(tax) < 11) {
            breaks <- as.numeric(names(tax))
            y <- tax/max(tax)
            interbreak <- min(diff(breaks)) * (length(tax) - 
                1)/41
            rect(breaks - interbreak, 0, breaks + interbreak, 
                y, col = hist.col)
        }
        else {
            h <- hist(x, breaks = breaks, plot = FALSE)
            breaks <- h$breaks
            nB <- length(breaks)
            y <- h$counts
            y <- y/max(y)
            rect(breaks[-nB], 0, breaks[-1], y, col = hist.col)
        }
        if (density) {
            tryd <- try(d <- density(x, na.rm = TRUE, bw = "nrd", 
                adjust = 1.2), silent = TRUE)
            if (!inherits(tryd, "try-error")) {
                d$y <- d$y/max(d$y)
                lines(d, col=dens.col)
            }
        }
        if (rug) 
            rug(x)
    }
    "panel.cor" <- function(x, y, prefix = "", ...) {
        usr <- par("usr")
        on.exit(par("usr"))
        par(usr = c(0, 1, 0, 1))
        if (is.null(wt)) {
            r <- cor(x, y, use = "pairwise", method = method)
        }
        else {
            r <- cor.wt(data.frame(x, y), w = wt[, c(1:2)])$r[1, 
                2]
        }
        txt <- format(c(round(r, digits), 0.123456789), digits = digits)[1]
        txt <- paste(prefix, txt, sep = "")
        if (stars) {
            pval <- r.test(sum(!is.na(x * y)), r)$p
            symp <- symnum(pval, corr = FALSE, cutpoints = c(0, 
                0.001, 0.01, 0.05, 1), symbols = c("***", "**", 
                "*", " "), legend = FALSE)
            txt <- paste0(txt, symp)
        }
        cex <- cex.cor * 0.8/(max(strwidth("0.12***"), strwidth(txt)))
        if (scale) {
            cex1 <- cex * (abs(r) + 0.6)/1.3# changed by MB from cex1 <- cex * abs(r)
            if (cex1 < 0.25) 
                cex1 <- 0.25
            text(0.5, 0.5, txt, cex = cex1)
        }
        else {
            text(0.5, 0.5, txt, cex = cex)
        }
    }
    "panel.smoother" <- function(x, y, pch = par("pch"), col.smooth = smooth.col, 
        span = 2/3, iter = 3, ...) {
        xm <- mean(x, na.rm = TRUE)
        ym <- mean(y, na.rm = TRUE)
        xs <- sd(x, na.rm = TRUE)
        ys <- sd(y, na.rm = TRUE)
        r = cor(x, y, use = "pairwise", method = method)
        if (jiggle) {
            x <- jitter(x, factor = factor)
            y <- jitter(y, factor = factor)
        }
        if (smoother) {
            smoothScatter(x, y, add = TRUE, nrpoints = 0)
        }
        else {
            if (show.points) 
                points(x, y, pch = pch, ...)
        }
        ok <- is.finite(x) & is.finite(y)
        if (any(ok)) {
            if (smooth & ci) {
                lml <- loess(y ~ x, degree = 1, family = "symmetric")
                tempx <- data.frame(x = seq(min(x, na.rm = TRUE), 
                  max(x, na.rm = TRUE), length.out = 47))
                pred <- predict(lml, newdata = tempx, se = TRUE)
                if (ci) {
                  upperci <- pred$fit + confid * pred$se.fit
                  lowerci <- pred$fit - confid * pred$se.fit
                  polygon(c(tempx$x, rev(tempx$x)), c(lowerci, 
                    rev(upperci)), col = adjustcolor(smooth.col, 
                    alpha.f = 0.5), border = NA)
                }
                lines(tempx$x, pred$fit, col = smooth.col, ...)
            }
            else {
                if (smooth) 
                lml <- loess(y ~ x, degree = 1, family = "symmetric")
                tempx <- data.frame(x = seq(min(x, na.rm = TRUE), 
                 max(x, na.rm = TRUE), length.out = 47))
                pred <- predict(lml, newdata = tempx, se = TRUE)
                lines(tempx$x, pred$fit, col = smooth.col, ...)
                  # original below
                  #lines(stats::lowess(x[ok], y[ok], f = span, 
                    #iter = iter), col = "pink")
            }
        }
        if (ellipses) 
            draw.ellipse(xm, ym, xs, ys, r, col.smooth = col.smooth, 
                ...)
    }
    "panel.lm" <- function(x, y, pch = par("pch"), col.lm = "red", 
        ...) {
        ymin <- min(y)
        ymax <- max(y)
        xmin <- min(x)
        xmax <- max(x)
        ylim <- c(min(ymin, xmin), max(ymax, xmax))
        xlim <- ylim
        if (jiggle) {
            x <- jitter(x, factor = factor)
            y <- jitter(y, factor = factor)
        }
        if (smoother) {
            smoothScatter(x, y, add = TRUE, nrpoints = 0)
        }
        else {
            if (show.points) {
                points(x, y, pch = pch, ylim = ylim, xlim = xlim, 
                  ...)
            }
        }
        ok <- is.finite(x) & is.finite(y)
        if (any(ok)) {
            lml <- lm(y ~ x)
            if (ci) {
                tempx <- data.frame(x = seq(min(x, na.rm = TRUE), 
                  max(x, na.rm = TRUE), length.out = 47))
                pred <- predict.lm(lml, newdata = tempx, se.fit = TRUE)
                upperci <- pred$fit + confid * pred$se.fit
                lowerci <- pred$fit - confid * pred$se.fit
                polygon(c(tempx$x, rev(tempx$x)), c(lowerci, 
                  rev(upperci)), col = adjustcolor("light grey", 
                  alpha.f = 0.8), border = NA)
            }
            if (ellipses) {
                xm <- mean(x, na.rm = TRUE)
                ym <- mean(y, na.rm = TRUE)
                xs <- sd(x, na.rm = TRUE)
                ys <- sd(y, na.rm = TRUE)
                r = cor(x, y, use = "pairwise", method = method)
                draw.ellipse(xm, ym, xs, ys, r, col.smooth = col.lm, 
                  ...)
            }
            abline(lml, col = col.lm, ...)
        }
    }
    "draw.ellipse" <- function(x = 0, y = 0, xs = 1, ys = 1, 
        r = 0, col.smooth, add = TRUE, segments = 51, ...) {
        angles <- (0:segments) * 2 * pi/segments
        unit.circle <- cbind(cos(angles), sin(angles))
        if (!is.na(r)) {
            if (abs(r) > 0) 
                theta <- sign(r)/sqrt(2)
            else theta = 1/sqrt(2)
            shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(c(theta, 
                theta, -theta, theta), ncol = 2, byrow = TRUE)
            ellipse <- unit.circle %*% shape
            ellipse[, 1] <- ellipse[, 1] * xs + x
            ellipse[, 2] <- ellipse[, 2] * ys + y
            if (show.points) 
                points(x, y, pch = 19, col = col.smooth, cex = 1.5)
            lines(ellipse, ...)
        }
    }
    "panel.ellipse" <- function(x, y, pch = par("pch"), col.smooth = "red", 
        ...) {
        segments = 51
        usr <- par("usr")
        on.exit(par("usr"))
        par(usr = c(usr[1] - abs(0.05 * usr[1]), usr[2] + abs(0.05 * 
            usr[2]), 0, 1.5))
        xm <- mean(x, na.rm = TRUE)
        ym <- mean(y, na.rm = TRUE)
        xs <- sd(x, na.rm = TRUE)
        ys <- sd(y, na.rm = TRUE)
        r = cor(x, y, use = "pairwise", method = method)
        if (jiggle) {
            x <- jitter(x, factor = factor)
            y <- jitter(y, factor = factor)
        }
        if (smoother) {
            smoothScatter(x, y, add = TRUE, nrpoints = 0)
        }
        else {
            if (show.points) {
                points(x, y, pch = pch, ...)
            }
        }
        angles <- (0:segments) * 2 * pi/segments
        unit.circle <- cbind(cos(angles), sin(angles))
        if (!is.na(r)) {
            if (abs(r) > 0) 
                theta <- sign(r)/sqrt(2)
            else theta = 1/sqrt(2)
            shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(c(theta, 
                theta, -theta, theta), ncol = 2, byrow = TRUE)
            ellipse <- unit.circle %*% shape
            ellipse[, 1] <- ellipse[, 1] * xs + xm
            ellipse[, 2] <- ellipse[, 2] * ys + ym
            points(xm, ym, pch = 19, col = col.smooth, cex = 1.5)
            if (ellipses) 
                lines(ellipse, ...)
        }
    }
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    if (missing(cex.cor)) 
        cex.cor <- 1
    for (i in 1:ncol(x)) {
        if (is.character(x[[i]])) {
            x[[i]] <- as.numeric(as.factor(x[[i]]))
            colnames(x)[i] <- paste(colnames(x)[i], "*", sep = "")
        }
    }
    n.obs <- nrow(x)
    confid <- qt(1 - alpha/2, n.obs - 2)
    if (!lm) {
        if (cor) {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, 
                lower.panel = panel.smoother, pch = pch, ...)
        }
        else {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.smoother, 
                lower.panel = panel.smoother, pch = pch, ...)
        }
    }
    else {
        if (!cor) {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm, 
                lower.panel = panel.lm, pch = pch, ...)
        }
        else {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, 
                lower.panel = panel.lm, pch = pch, ...)
        }
    }
    }
  # add image function (from https://stackoverflow.com/questions/27800307/adding-a-picture-to-plot-in-r)
    addImg <- function(
      obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
      x = NULL, # mid x coordinate for image
      y = NULL, # mid y coordinate for image
      width = NULL, # width of image (in x coordinate units)
      interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
    ){
      if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
      USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
      PIN <- par()$pin # The current plot dimensions, (width, height), in inches
      DIM <- dim(obj) # number of x-y pixels for the image
      ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
      WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
      HEIi <- WIDi * ARp # height in inches
      HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
      rasterImage(image = obj, 
        xleft = x-(width/2), xright = x+(width/2),
        ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
        interpolate = interpolate)
    }
# load data
  source(here::here('R/DAT_prepare.R'))      

  a[, animal := bird_ID]
  a[, id := bird_ID]
  b[, animal := bird_ID]
  b[, id := bird_ID]

  # prepare motility
    dd = d[!Morph%in%'Zebra finch']
    ddl = data.table(melt(dd[,.(bird_ID,month,Morph,age,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","motileCount"), variable.name = "Motility"))
    ddl[Motility == 'VAP' ,mot:='Average path']
    ddl[Motility == 'VCL' ,mot:='Curvilinear']
    ddl[Motility == 'VSL' ,mot:='Straight line']
    ddl[, animal := bird_ID]
    

    # use June values and for 4 males without June, May
      dd1 = dd[month == 'June']
      dd2 = dd[month == 'May']
      ddx = rbind(dd1,dd2[!bird_ID%in%dd1$bird_ID])
      ddxl = data.table(melt(ddx[,.(bird_ID,month,Morph,age,motileCount,VAP,VSL,VCL)], id.vars = c("bird_ID","month","Morph","age","motileCount"), variable.name = "Motility"))
      ddxl[Motility == 'VAP' ,mot:='Average path']
      ddxl[Motility == 'VCL' ,mot:='Curvilinear']
      ddxl[Motility == 'VSL' ,mot:='Straight line']
      ddxl[, animal := bird_ID]
  
  # pedigree
    p = fread('Data/Dat_parentage.txt',na.strings = c(NA_character_, ""))
    setnames(p, old = c('bird_ID', 'Mother', 'Father'), new = c('id', 'dam', 'sire'))
    p = p[!is.na(id),.(id,sire,dam)]
    pc = tidyped(p)
    #visped(pc) 
    p2 <- with(data.frame(pc), pedigreemm::pedigree(sire=Sire, dam=Dam, label=Ind))

# pairs.panels - single data
  bw_ = bw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')]

  png("Outputs/Fig_S1.png", width = 8*2, height = 8*2, units = "cm", res = 600)
  pairs.panels_MB(bw_,
             #oma=c(0,0,0,0),
             gap = 0,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = TRUE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             rug = FALSE,
             pch = 21,           # pch symbol
             bg = fills[bw$Morph],  # Background color of the symbol (pch 21 to 25)
             col = cols[bw$Morph],
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             #hist.col = 3,       # Histograms color
             stars = FALSE,       # If TRUE, adds significance level with stars
             ci = FALSE)          # If TRUE, adds confidence intervals   
  mtext("Single sperm", side=3, line=3)
  par(xpd=TRUE)
  y_ = -0.18
  points(c(0.17, 0.17+0.23, 0.17+0.45), c(y_,y_,y_), pch = 21, bg = fills, col = cols, cex = 1)
  addImg(img_i, x = 0.11, y = y_, width = 0.075)
  addImg(img_s, x = 0.11+0.23, y = y_, width = 0.075)
  addImg(img_fc, x = 0.11+0.45, y = y_, width = 0.075)
  dev.off()
# pairs.panels - averages data
  aw_ = aw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')]

  png("Outputs/Fig_S1b.png", width = 8*2, height = 8*2, units = "cm", res = 600)
  pairs.panels_MB(aw_,
             #oma=c(0,0,0,0),
             gap = 0,
             smooth = TRUE,      # If TRUE, draws loess smooths
             scale = TRUE,      # If TRUE, scales the correlation text font
             density = TRUE,     # If TRUE, adds density plots and histograms
             ellipses = FALSE,    # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             rug = FALSE,
             pch = 21,           # pch symbol
             bg = fills[aw$Morph],  # Background color of the symbol (pch 21 to 25)
             col = cols[aw$Morph],
             lm = FALSE,         # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = TRUE,         # If TRUE, reports correlations
             jiggle = FALSE,     # If TRUE, data points are jittered
             factor = 2,         # Jittering factor
             #hist.col = 3,       # Histograms color
             stars = FALSE,       # If TRUE, adds significance level with stars
             ci = FALSE)          # If TRUE, adds confidence intervals   
  mtext("Male means", side=3, line=3)
  par(xpd=TRUE)
  y_ = -0.18
  points(c(0.17, 0.17+0.23, 0.17+0.45), c(y_,y_,y_), pch = 21, bg = fills, col = cols, cex = 1)
  addImg(img_i, x = 0.11, y = y_, width = 0.075)
  addImg(img_s, x = 0.11+0.23, y = y_, width = 0.075)
  addImg(img_fc, x = 0.11+0.45, y = y_, width = 0.075)
  dev.off()

# not used - pair plot v2
 chart.Correlation.nostars(bw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')], histogram=TRUE, pch =20)    
  mtext("Single sperm", side=3, line=3)

# not used - ggpairs
  ggpairs(bw[, c('Acrosome', 'Nucleus','Head', 'Midpiece', 'Tail', 'Flagellum','Total', 'Midpiece_rel', 'Flagellum_rel')])

# END