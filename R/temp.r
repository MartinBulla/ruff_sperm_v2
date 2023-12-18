  
  fig <- function(p, fig, w, h, u = 'cm', ...) {
  png(file = glue("./Outputs/{fig}.png"), width = w, height = h, units = u, res = 300)
   par(mar = c(0.25, 0, 0, 1.2)
   , oma = c(2.1, 2.2, 0.2, 2.4), ps = 12, mgp = c(1.2, 0.35, 0), las = 1, cex = 1, col.axis = "black", font.main = 1, col.lab = "black", col.main = "black", fg = "black", cex.lab = 0.6, cex.main = 0.7, cex.axis = 0.5, tcl = -0.1, bty = "n", xpd = TRUE, lwd = 0.5)
    plot(p)
    dev.off()
  }
  
  # CV
  m = lm(scale(Length_avg) ~ Morph, a[part == i])

  a_cv = a[, raster::cv(Length_avg), by = list(Morph, part)]
  setnames(a_cv, 'V1', 'cv')

  ggplot(a_cv, aes(x = Morph, y = cv, color = part)) + geom_point()
  
  cv_ =  b[, cv(Length_µm), by = list(bird_ID, part, Morph, loc_avi)]


# aviaries

require(effects)
d[, motileCount_ln := scale(log(motileCount))]
dd <- d[!Morph %in% "Zebra finch"]

# use June values and for 4 males without June, May
dd1 <- dd[month == "June"]
dd2 <- dd[month == "May"]
ddx <- rbind(dd1, dd2[!bird_ID %in% dd1$bird_ID])

ddx[, access_to_females:=treat]
ddx[location == 9 & avi %in% c(1, 2), access_to_females := "partial"]
ddx[access_to_females == "prison", access_to_females := 'none']
ddx[access_to_females == "mixed", access_to_females := 'continuous']
ddx[, access_to_females := factor(access_to_females, levels = c("none", "partial", "continuous"))]

m <- lm(scale(VSL) ~ scale(log(motileCount)) + access_to_females + Morph, ddx)
png(file = here::here("Outputs/VSL_femaleAccess.png"), width = 25, height = 8, units = "cm", res = 300)
plot(allEffects(m), rows = 1, cols = 3)
dev.off()

m <- lm(scale(VAP) ~ scale(log(motileCount)) + access_to_females + Morph, ddx)
png(file = here::here("Outputs/VAP_femaleAccess.png"), width = 25, height = 8, units = "cm", res = 300)
plot(allEffects(m), rows = 1, cols = 3)
dev.off()

m <- lm(scale(VCL) ~ scale(log(motileCount)) + access_to_females
 + Morph, ddx)
png(file = here::here("Outputs/VCL_femaleAccess.png"), width = 25, height = 8, units = "cm", res = 300)
plot(allEffects(m), rows = 1, cols = 3)
dev.off()

ggplot()

tst = ddx[,.N, by = list(location, avi)]
setnames(tst, 'N', 'n_males')

ddx = merge(ddx,tst)

m <- lm(scale(VSL) ~ scale(log(motileCount)) + treat2 + Morph, ddx)

m <- lm(scale(VSL) ~ scale(log(motileCount)) , ddx)
plot(allEffects(m))
summary(m)


tsz <- ddx[, .N, by = list(location, treat2, avi, Morph)]
reshape(tsz, idvar = c("location", "treat2", "avi"), timevar = "Morph", direction = "wide")
