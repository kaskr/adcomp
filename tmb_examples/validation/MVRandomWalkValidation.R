## MVRandomWalkValidation
##
## Estimate and validate a multivariate random walk model with correlated
## increments and correlated observations.
##
## Compare Thygesen et al (submitted, 2016): Validation of state space
## models fitted as mixed effects models
## Casper W. Berg and Kasper Kristensen, 2016

library(TMB)
compile("MVRandomWalkValidation.cpp")
dyn.load(dynlib("MVRandomWalkValidation"))

library(MASS)
simdata <- function(stateDim = 8, timeSteps = 100, rho = 0.9, rho2 = 0.9, sds = seq(0.5, 
    0.51, length = stateDim), sdObs = rep(2, stateDim), plot = FALSE, seed = 143) {
    set.seed(seed)
    
    corrMat <- matrix(0, stateDim, stateDim)
    corrMat2 <- matrix(0, stateDim, stateDim)
    for (i in 1:stateDim) {
        for (j in 1:stateDim) {
            corrMat[i, j] <- rho^abs(i - j)
            corrMat2[i, j] <- rho2^abs(i - j)
        }
    }
    Sigma <- corrMat * (sds %o% sds)
    Sigma2 <- corrMat2 * (sdObs %o% sdObs)
    d <- matrix(NA, timeSteps, stateDim)
    obs <- d
    ## init state
    d[1, ] <- rnorm(stateDim)
    i <- 1
    obs[i, ] <- d[i, ] + mvrnorm(1, rep(0, stateDim), Sigma = Sigma2)
    for (i in 2:timeSteps) {
        d[i, ] <- d[i - 1, ] + mvrnorm(1, rep(0, stateDim), Sigma = Sigma)
        obs[i, ] <- d[i, ] + mvrnorm(1, rep(0, stateDim), Sigma = Sigma2)
    }
    if (plot) {
        matplot(d, type = "l")
        matpoints(obs)
    }
    return(list(d = d, obs = obs, sds = sds, sdObs = sdObs))
}

stateDim <- 4
timeSteps <- 100

## Note: For reproducibility across machines ('mvrnorm' uses eigen
## decompostion)
if(TRUE) {
    ## Use cached data
    sim <- dget("MVRandomWalkValidation_data.R")
} else {
    ## Simulate data
    sim <- simdata(stateDim, timeSteps, seed = 12345)
}

d <- sim$d
obs <- sim$obs

data <- list(obs = t(obs))
parameters <- list(transf_rho = 0.1, transf_rhoObs = 0.1, logsds = sim$sds * 0, logsdObs = sim$sdObs * 
    0, u = data$obs * 0)
## Optimize using all data
obj <- MakeADFun(data, parameters, random = "u", DLL = "MVRandomWalkValidation")  ##,map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

pl <- obj$env$parList()

OSA <- oneStepPredict(obj, "obs", method = "fullGaussian")

res <- OSA$residual
dim(res) <- dim(data$obs)
res <- t(res)

## Bubble plot
bp <- function(x, y, v, scale = 3, ...) {
    plot(x, y, cex = sqrt(abs(v)) * scale, col = ifelse(v < 0, "tomato2", "blue"), 
        pch = ifelse(v < 0, 16, 1), ...)
    points(x[v > 0], y[v > 0], cex = sqrt(v[v > 0]) * scale, col = "blue", pch = 1, 
        ...)
}

## transform real axis to ]-1,1[
transf <- function(x) 2/(1 + exp(-2 * x)) - 1

## residuals for each of the 4 runs
OSA.resid <- array(NA, dim = c(4, timeSteps, stateDim))
OSA.resid[1, , ] <- res

## matrix with true/estimated parameters + AIC
param.tab <- matrix(NA, 5, length(opt$par) + 1)
param.tab[1, ] <- c(0.9, 0.9, seq(0.5, 1, length = stateDim), sdObs = rep(2, stateDim), 
    NA)
param.tab[2, ] <- unlist(c(transf(unlist(pl[1:2])), exp(unlist(pl[3])), exp(unlist(pl[4])), 
    2 * length(opt$par) + 2 * opt$objective))

########################################################################################
## Repeat with wrong model 1 (process AND observations assumed uncorrelated)
########################################################################################

map <- list(transf_rhoObs = factor(NA), transf_rho = factor(NA))
parameters$transf_rhoObs <- 1e-06
parameters$transf_rho <- 1e-06

obj <- MakeADFun(data, parameters, random = "u", DLL = "MVRandomWalkValidation", 
    map = map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

pl <- obj$env$parList()

OSA <- oneStepPredict(obj, "obs", method = "fullGaussian")

res <- OSA$residual
dim(res) <- dim(data$obs)
res <- t(res)

OSA.resid[2, , ] <- res

param.tab[3, ] <- unlist(c(transf(unlist(pl[1:2])), exp(unlist(pl[3])), exp(unlist(pl[4])), 
    2 * length(opt$par) + 2 * opt$objective))

########################################################################################
## Repeat with wrong model 2 (observations assumed uncorrelated)
########################################################################################

map <- list(transf_rhoObs = factor(NA))
parameters$transf_rhoObs <- 1e-06

obj <- MakeADFun(data, parameters, random = "u", DLL = "MVRandomWalkValidation", 
    map = map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

pl <- obj$env$parList()

OSA <- oneStepPredict(obj, "obs", method = "fullGaussian")

res <- OSA$residual
dim(res) <- dim(data$obs)
res <- t(res)
OSA.resid[3, , ] <- res

param.tab[4, ] <- unlist(c(transf(unlist(pl[1:2])), exp(unlist(pl[3])), exp(unlist(pl[4])), 
    2 * length(opt$par) + 2 * opt$objective))

########################################################################################
## Repeat with wrong model 3 (process assumed uncorrelated)
########################################################################################

map <- list(transf_rho = factor(NA))
parameters$transf_rho <- 1e-06

obj <- MakeADFun(data, parameters, random = "u", DLL = "MVRandomWalkValidation", 
    map = map)
opt <- nlminb(obj$par, obj$fn, obj$gr)

pl <- obj$env$parList()

OSA <- oneStepPredict(obj, "obs", method = "fullGaussian")

res <- OSA$residual
dim(res) <- dim(data$obs)
res <- t(res)

OSA.resid[4, , ] <- res

param.tab[5, ] <- unlist(c(transf(unlist(pl[1:2])), exp(unlist(pl[3])), exp(unlist(pl[4])), 
    2 * length(opt$par) + 2 * opt$objective))

############################# Plot results / make table

pdf("bubbles.pdf")
par(mfrow = c(4, 1), mar = c(4, 3, 1, 1))
mains <- c(paste("Model", 1:4))
for (i in 1:4) {
    bp(row(OSA.resid[i, , ]), col(OSA.resid[i, , ]), OSA.resid[i, , ], main = mains[i], 
        ylim = c(0, stateDim + 1), xlab = "t", scale = 2)
}
dev.off()

pdf("qqplots.pdf")
par(mfrow = c(4, 1), mar = c(4, 1, 1, 1))
for (i in 1:4) {
    qqnorm(OSA.resid[i, , ])
    abline(0, 1)
}
dev.off()


op <- par
## Try to squeeze it into one plot
mult <- 1.25
pdf("bubblesacf.pdf", width = 12 * mult, height = 8 * mult)
par(mar = c(4, 1, 4, 2))
mat <- matrix(c(1, 1, 1, 2, 3, 4, 4, 4, 5, 6, 7, 7, 7, 8, 9, 10, 10, 10, 11, 12), 
    4, 5, byrow = TRUE)
nf <- layout(mat)

acf.pval <- function(x) {
    acfC <- x$acf[-1]
    acfC <- na.omit(acfC)
    2 - 2 * pnorm(abs(acfC) * sqrt(x$n.used))
}
calc.limval <- function(p, acf) qnorm((2 - p)/2)/sqrt(acf$n.used)

KSvals <- numeric(4)

for (i in 1:4) {
    
    tmp <- OSA.resid[i, , ]
    tmp <- rbind(tmp, matrix(NA, 11, 4))
    
    myacf <- acf(as.vector(tmp), lag.max = 10, na.action = na.pass, plot = FALSE)
    xl <- ifelse(i == 4, "t", "")
    
    bp(row(OSA.resid[i, , ]), col(OSA.resid[i, , ]), OSA.resid[i, , ], main = mains[i], 
        ylim = c(0, stateDim + 1), scale = 2, ylab = "", xlab = xl)
    
    kspval <- round(ks.test(OSA.resid[i, , ], pnorm)$p.value, 3)
    KSvals[i] <- kspval
    text(95, 4.75, labels = paste0("KS p-value: ", kspval), cex = 1)
    xl <- ifelse(i == 4, "Lag", "")
    
    if (i == 1) 
        acf(as.vector(tmp), lag.max = 10, na.action = na.pass, main = "ACF time direction", 
            xlab = xl) else acf(as.vector(tmp), lag.max = 10, na.action = na.pass, main = "", xlab = xl)
    
    tmp <- OSA.resid[i, , ]
    tmp <- cbind(tmp, matrix(NA, nrow(tmp), 4))
    
    myacf <- acf(as.vector(t(tmp)), lag.max = 3, na.action = na.pass, plot = FALSE)
    
    if (i == 1) 
        acf(as.vector(t(tmp)), lag.max = 3, na.action = na.pass, main = "ACF state direction", 
            xlab = xl) else 
    acf(as.vector(t(tmp)), lag.max = 3, na.action = na.pass, main = "", xlab = xl)
    
}
dev.off()

if(FALSE) {
    library(xtable)
    param.tab <- cbind(param.tab, c(NA, KSvals))
    ## Table
    rownames(param.tab) <- c("True", paste("Model", 1:4))
    colnames(param.tab) <- c("$\\rho_X$",
                             "$\\rho_Y$",
                             paste("$\\sigma_{X",
                                   paste0(1:stateDim, "}$"), sep = ","),
                             paste("$\\sigma_{Y",
                                   paste0(1:stateDim, "}$"), sep = ","),
                             "AIC", "KS test")
    cat(print(xtable(t(param.tab),
                     caption = paste(
                         "True/estimated parameters, AIC, and",
                         "p-values from the Kolomogorov-Smirnov",
                         "tests for normality of the residuals"),
                     align = rep("c", nrow(param.tab) + 1),
                     digits = 3),
              type = "latex",
              sanitize.text.function = function(x)x),
        file = "table.tex")
}
