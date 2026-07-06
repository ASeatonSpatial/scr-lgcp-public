library(ascr)
library(RColorBrewer)
library(sp)
library(viridis)
library(Rcpp)
load("gibbon-data.RData")
sourceCpp("fitting-functions/fitting-functions.cpp")
source("fitting-functions/fitting-functions.R")

## A mask with only dense forest and very dense forest, and converting
## distances to km.
mask.forest.df <- vector(mode = "list", length = length(mask.full.df))
mask.forest <- vector(mode = "list", length = length(mask.full))
for (i in 1:length(mask.full.df)){
    mask.full.df[[i]]$VILLAGE <- mask.full.df[[i]]$VILLAGE/1000
    mask.forest.df[[i]] <- mask.full.df[[i]][mask.full.df[[i]]$FOREST_COVER != "NOT_DENSE", ]
    mask.forest[[i]] <- mask.full[[i]][mask.full.df[[i]]$FOREST_COVER != "NOT_DENSE", ]
}
fine.mask.df$VILLAGE <- fine.mask.df$VILLAGE/1000

## SCR model with an inhomogeneous Poisson point process, using both
## forest cover and distance to the nearest village as a covariate.
ihpp.fit <- fit.ascr(capt, traps, mask.forest, detfn = "hhn",
                       survey.length = rep(3, length(capt)),
                       ihd.opts = list(model = ~ FOREST_COVER + VILLAGE, covariates = mask.forest.df))
summary(ihpp.fit)

## Estimated density across PPWS.
is.forest <- fine.mask.df$FOREST_COVER != "NOT_DENSE"
ihpp.fit.pred <- numeric(nrow(fine.mask.df))
ihpp.fit.pred[is.forest] <- predict(ihpp.fit, newdata = fine.mask.df[is.forest, ])
## Converting to calling gibbon groups per 100 square km.
ihpp.fit.pred <- ihpp.fit.pred*100*100

## Setting things up for the model with penalised regression splines.
pred.df <- data.frame(x = 702553.4, y = 1404628.3,
                      FOREST_COVER = "VERY_DENSE",
                      VILLAGE = 5)
## Start values for optimisation.
start.par <- list(D_betas = log(0.001),
                  link_lambda0 = log(50),
                  link_sigma = log(500))

## Fitting model with smooth spatial and village effects.
fit.sm <- fit.scr.smooth(capt, traps, mask.full,
                         model = ~ FOREST_COVER +
                             s(x, y, k = 25) +
                             s(VILLAGE, k = 7),
                         mask.df = mask.full.df, pred.df = pred.df,
                         unsuitable = list(FOREST_COVER = "NOT_DENSE"),
                         n.occasions = 1, start.par = start.par, detfn = "HHN",
                         tmb.dir = "fitting-functions/tmb")
summary(fit.sm$sdrep)

## For comparison, another model with a log-linear effect of village.
fit.sm2 <- fit.scr.smooth(capt, traps, mask.full,
                              model = ~ FOREST_COVER +
                                s(x, y, k = 25) +
                                VILLAGE,
                              mask.df = mask.full.df, pred.df = pred.df,
                              unsuitable = list(FOREST_COVER = "NOT_DENSE"),
                              n.occasions = 1, start.par = start.par, detfn = "HHN",
                              tmb.dir = "fitting-functions/tmb")
summary(fit.sm2$sdrep)

## Estimated smooth village effect for 0-15 km for both models.
village.xx <- seq(0, 15, length.out = 100)
pred.village.sm <- predict.effect(fit.sm, data.frame(VILLAGE = village.xx),
                                  which.smooth = 2)
pred.village.sm <- pred.village.sm - mean(pred.village.sm)
pred.village.lin <- fit.sm2$fit$par[1] + fit.sm2$fit$par[3]*village.xx
pred.village.lin <- pred.village.lin - mean(pred.village.lin)


## Estimated smooth spatial effect over the whole region.
pred.spatial <- predict.effect(fit.sm, fine.mask.df, which.smooth = 1)
## Predicting density over the whole region.                               
pred.D <- predict.D(fit.sm, fine.mask.df)
## Converting to gibbon groups per 100 sq km.
pred.D <- pred.D*100*100
## Average estimated density.
mean(pred.D)

## New data frame for plotting.
plot.mask.df <- fine.mask.df
plot.mask.df$FOREST_COVER <- factor(plot.mask.df$FOREST_COVER,
                                    levels = c("NOT_DENSE", "DENSE", "VERY_DENSE"),
                                    labels = c("Not dense", "Dense", "Very dense"))
## Data frame of trap locations.
trap.df <- data.frame(do.call(rbind, traps))
colnames(trap.df) <- c("x", "y")
## Calculating trap centroids for each session.
midtraps <- t(sapply(traps, function(x) x[2, ]))
trap.dets <- sapply(capt, function(x) nrow(x$bincapt))

## Plotting comparisons between models.
par(mfrow = c(3, 2), mar = c(1, 0, 2, 7), oma = c(0, 0, 0, 0))
pal <- grDevices::colorRampPalette(c("grey92", "#3B77B6", "#08306B"))
dens.cols <- pal(100)
plot.surf(plot.mask.df$x, plot.mask.df$y, plot.mask.df$VILLAGE,
          traps = midtraps, trap.dets = trap.dets, scale = 10000,
          main = "(A) Distance to village (km)")
plot.surf(plot.mask.df$x, plot.mask.df$y, plot.mask.df$FOREST_COVER,
          traps = midtraps, trap.dets = trap.dets, main = "(B) Forest cover")
plot.surf(plot.mask.df$x, plot.mask.df$y, ihpp.fit.pred,
          traps = midtraps, trap.dets = trap.dets,
          cols = dens.cols,
          zlim = c(0, max(pred.D)),
          main = "(C) Estimated density from Model SCR-I:\n            density ~ forest + village",
          sub = bquote("Average estimated density:" ~ .(round(mean(ihpp.fit.pred), 1)) ~
                      "calling groups per 100" ~ km^2))
plot.surf(plot.mask.df$x, plot.mask.df$y, pred.D,
          traps = midtraps, trap.dets = trap.dets,
          cols = dens.cols,
          zlim = c(0, max(pred.D)),
          main = "(D) Estimated density from Model SCR-S:\n            density ~ forest + s(village) + s(x, y)",
          sub = bquote("Average estimated density:" ~ .(round(mean(pred.D), 1)) ~
                      "calling groups per 100" ~ km^2))
spar <- par(mar = c(4, 4, 2, 3))
plot.new()
plot.window(xlim = range(village.xx), ylim = range(c(pred.village.lin, pred.village.sm)))
box()
axis(1)
axis(2)
lines(village.xx, pred.village.sm)
lines(village.xx, pred.village.lin, lty = "dashed")
title(xlab = "Distance to nearest village (km)", ylab = "Effect on log-density")
mtext("(E) Estimated village effects", at = grconvertX(0, from = "ndc", to = "user"),
      side = 3, line = 0.5, adj = 0,
      cex = 0.66,
      font  = par("font.main"),
      col   = par("col.main"),
      family = par("family"))
legend("bottomright", lty = c("dashed", "solid"), legend = c("Model SCR-I: log-linear", "Model SCR-S: spline"))
par(spar)
plot.surf(plot.mask.df$x, plot.mask.df$y, pred.spatial,
          traps = midtraps, trap.dets = trap.dets,
          main = "")
mtext("(F) Estimated spatial effect from Model SCR-S", at = grconvertX(0.5, from = "ndc", to = "user"),
      side = 3, line = 0.5, adj = 0,
      cex = 0.66,
      font  = par("font.main"),
      col   = par("col.main"),
      family = par("family"))


