library(ascr)
library(geoR)
library(secr)
library(CircStats)
sourceCpp("fitting-functions/fitting-functions.cpp")
source("fitting-functions/fitting-functions.R")

## Creating a big grid of mask points.
x.range <- c(0, 100000)
y.range <- c(0, 100000)
xs <- seq(x.range[1], x.range[2], length.out = 100)
ys <- seq(y.range[1], y.range[2], length.out = 100)
xx <- rep(xs, length(ys))
yy <- rep(ys, each = length(xs))
mask.full <- cbind(xx, yy)
mask.area <- sqrt((max(x.range) - min(x.range))^2) *
  sqrt((max(y.range) - min(y.range))^2)

## Detection function parameter values. Order is g0, sigma, kappa.
pars <- c(0.95, 1500, 50)

## Generating session locations.
sess.xs <- seq(x.range[1] + 3*pars[2] , x.range[2] - 3*pars[2], length.out = 4)
sess.ys <- seq(y.range[1] + 3*pars[2], y.range[2] - 3*pars[2], length.out = 4)
sess.locs <- expand.grid(sess.xs, sess.ys)
n.sessions <- nrow(sess.locs)

## Generating trap locations.
traps <- vector(mode = "list", length = n.sessions)
for (i in 1:n.sessions) {
  traps[[i]] <- cbind(
    x = rep(sess.locs[i, 1], 3) + c(-1000, 0, 1000),
    y = rep(sess.locs[i, 2], 3)
  )
}
traps.all = read.traps(data = data.frame(do.call(rbind, traps)))

## Simulating a spatial covariate.
set.seed(625624)
cov <- grf(nrow(full.mask),
  grid = mask.full, xlims = x.range,
  ylims = y.range, cov.pars = c(2, 100000*1.2)
)$data

## A function relating the covariate to density via a saturating
## response.
D.calc <- function(df) {
  b0 <- -log(44426388/200000)
  alpha <- 1
  gamma <- 3
  tau <- -3.164784
  cov <- df$cov
  exp(b0 + alpha * (1 - exp(-gamma * (cov - tau))))
}

## Creating masks for each session for model fitting.
mask <- vector(mode = "list", length = n.sessions)
for (i in 1:n.sessions) {
  mask[[i]] <- make.mask(traps = data.frame(traps[[i]]),
                         buffer = 4*pars[2],
                         type = "trapbuffer")
}

## Creating covariates for model fitting by finding the nearest mask
## point.
mask.cov <- vector(mode = "list", length = n.sessions)
for (i in 1:n.sessions) {
  mask.cov.dists <- calc.dists(mask.full, mask[[i]])
  closest.cell <- apply(mask.cov.dists, 1, function(x) which(x == min(x))[1])
  mask.cov[[i]] <- data.frame(cov = cov[closest.cell])
}

## Simulating detection data.
sim.obj <- sim.capt.cov(
    traps = traps, mask = mask.full,
    df = data.frame(cov = cov), D.calc = D.calc,
    pars = pars,
    detfn = "HN"
)

## Fitting a model with a homogeneous Poisson process.
fit.a <- fit.ascr(capt = sim.obj$capt, traps = traps, mask = mask,  detfn = "hn")

## Fitting a mdoel with an inhomogeneous Poisson process, with a
## log-linear effect of the covariate.
fit.ih <- fit.ascr(capt = sim.obj$capt, traps = traps, mask = mask,
                   detfn = "hn", sv <- as.list(coef(fit.a)[-1]),
                   ihd.opts = list(model = ~ cov, covariates = mask.cov))

## Start value for models with smooth effects, based on the fitted
## models above.
start.par <- list(D_betas = log(coef(fit.a, "D")),
                  link_g0 = log(coef(fit.a, "g0") / (1 - coef(fit.a, "g0"))),
                  link_sigma = log(coef(fit.a, "sigma")),
                  link_kappa = log(coef(fit.a, "kappa")))

## Fitting a model with a smooth effect of the covariate.
fit.sm.cov <- fit.scr.smooth(sim.obj$capt, traps, mask,
                             model = ~ s(cov, k = 5), mask.df = mask.cov,
                             pred.df = data.frame(x = 50000, y = 50000, cov = 0),
                             start.par = start.par,
                             detfn = "HN", tmb.dir = "fitting-functions/tmb")

## Fitting a model with a smooth effect over space.
fit.sm.space <- fit.scr.smooth(sim.obj$capt, traps, mask,
                               model = ~ s(x, y, k = 25), mask.df = mask.cov,
                               pred.df = data.frame(x = 50000, y = 50000, cov = 0),
                               start.par = start.par,
                               detfn = "HN", tmb.dir = "fitting-functions/tmb")

## Predictions over space for the models.
pred.D.map.ih <- predict(fit.ih, newdata = data.frame(cov = cov))
pred.D.map.cov <- predict.D(fit.sm.cov, data.frame(x = mask.full[, 1],
                                                   y = mask.full[, 2],
                                                   cov = cov))
pred.D.map.space <- predict.D(fit.sm.space, data.frame(x = mask.full[, 1],
                                                       y = mask.full[, 2],
                                                       cov = cov))

## Prediction for the covariate for the two models that used it.
pred.D.cov <- predict.D(fit.sm.cov, data.frame(x = rep(50000, 1000),
                                               y = rep(50000, 1000),
                                               cov = cov.xx))
pred.D.cov.ih <- predict(fit.ih, newdata = data.frame(cov = cov.xx))

## These densities are arbitrary in terms of value; we could report in
## density per square metre, per hectare, per square km, or
## whatever). For nicer plotting and reporting in the paper, we scaled
## them so that true density varies from 0 to just over 1.

## Integrated squared error for each model.
ISE.cov <- sum(cell.area * (pred.D.map.cov - D)^2)
ISE.space <- sum(cell.area * (pred.D.map.space - D)^2)
ISE.a <- sum(cell.area * (fit.a$coeflist$D - D)^2)
ISE.ih <- sum(cell.area * (pred.D.map.ih - D)^2)

## Abundance estimation squared error.
## Covariate smooth model.
Nhat.cov <- sum(cell.area * pred.D.map.cov / 10000)
NSE.cov <- (Ntrue - Nhat.cov)^2

## Spatial smooth model
Nhat.space <- sum(cell.area * pred.D.map.space / 10000)
NSE.space <- (Ntrue - Nhat.space)^2

## Homogeneous model
Nhat.a <- sum(cell.area * fit.a$coeflist$D)
NSE.a <- (Ntrue - Nhat.a)^2

## Inhomogeneous model.
Nhat.ih <- sum(cell.area*pred.D.map.ih / 10000)
NSE.ih <- (Ntrue - Nhat.ih)^2
