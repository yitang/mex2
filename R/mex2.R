## main funciton 1
## estimate extreme dependence parameter alpha and beta
## calcualte the residual matrix Z
mex2 <- function(dt, dqu = 0.995){
    dt.laplace2 <- dt[, lapply(.SD, toLaplace)]
    y.cond <- dt.laplace2[[1]]
    y.cond.th <- quantile(y.cond, dqu)
    ind <- y.cond >= y.cond.th
    y.cond <- y.cond[ind]
    dt.laplace2 <- dt.laplace2[ind, ]
    dt.dep2 <- dt.laplace2[, lapply(.SD, function(x)
                                    my.qfun(y.cond, x)), .SDcols = -1]
    tmp <- rbind(dt.dep2[1:4, ], dt.laplace2[, -1, with = F])
    z.mat <- tmp[, lapply(.SD, function(x)
                          getZMatrix(y.cond = y.cond, x[-c(1:4)], para = x[1:4]))]
    ## format output
    coef <- data.table(id = names(dt.dep2), t(dt.dep2[1:2, ]))
    setnames(coef, names(coef)[2:3], c("a", "b"))
    HT.model <- list(coef = coef,
                     z = z.mat, y = dt.laplace2, y.threshold = y.cond.th)
    return(HT.model)
}

## main function 2
## diagnosis plots
diag.mex <- function(HT.model){
    ## cast
    y.big <- HT.model$y
    coef <- HT.model$coef
    z <- HT.model$z
    ## residual pltos
    y0 <- y.big[[1]]
    zz <- melt(z, variable.name = "id", value.name = "z", variable.factor = F)
    zz[, y0 := y0]
    p0 <- ggplot() + geom_point(data = zz, aes(y0, z)) + facet_wrap(~id)
    ## condifdence interval plot
    yy <- melt(y.big[, -1, with=F], variable.name = "id", value.name = "y", variable.factor = F)
    yy[, y0 := y0]
    band <- copy(zz)
    band[, lower := quantile(z, 0.025), by = id]
    band[, upper := quantile(z, 0.975), by = id]
    setkey(band, id)
    band <- unique(band)
    band <- merge(coef, band, by = "id")[, list(id, a, b, lower, upper)]
    y.seq <- seq(min(y0), max(y0), len = 1e3)
    band2 <- band[, list(y.seq,
                         lower = lower * y.seq^b + y.seq * a,
                         upper = upper * y.seq^b + y.seq * a),
                  by = id]
    p1 <- ggplot() + geom_point(data = yy, aes(y0, y)) + geom_ribbon(data = band2, aes(y.seq, ymin = lower, ymax = upper), alpha = 0.2) + facet_wrap(~id)
    return(list(p0, p1))
}


## main funciton 3
## simulate
samp.cond <- function(HT.model){
    ## cast
    z <- HT.model$z
    alpha <- HT.model$coef$a
    beta <- HT.model$coef$b
    y0 <- HT.model$y.threshold
    ## simulate
    alpha <- y0*alpha
    beta <- y0^beta
    tmp <- as.data.table(rbind(rbind(alpha,beta), z))
    sim <- tmp[, lapply(.SD, function(x)
                        x[-c(1:2)] * x[2] + x[1])]
    sim <- data.table(y0 = y0, sim)
    return(sim)
}



## tranfer x to laplapce distribution, using expiral distirbution
toLaplace <- function(x, mu = 0, b = 1){
    ## u <- rank(x, ties = "random") / (1+length(x))
    u <- rank(x) / (1+length(x))
    y <-     mu - b * sign(u - 0.5) * log(1 - 2 * abs(u - 0.5))
    return(y)
}

## laplace distirbution functions
qLaplace <- function(p, mu = 0, b = 1){
     ## mu is locaton param, b is scale para
    mu - b * sign(p - 0.5) * log(1 - 2 * abs(p - 0.5))
}
pLaplace <- function(x, mu = 0, b = 1){
    ## mu is location para, b is scale para
    1/2 + 1/2 * sign(x - mu) * (1  - exp(- abs(x - mu) / b))
}

## tranfrom laplace (y) to x distributiuon, using GPD extrpolate
y2x <- function(y, x, qu, coef){
    u <- pLaplace(y)
    xx <- rep(NA, len = length(u))
    threshold <- quantile(x, qu)
    ind <- u <= qu
    xx[ind] <- quantile(x, u[ind])
    if (any(!ind )){
        sig <- coef[2]
        xi <- coef[3]
        xx[!ind] <- texmex::qgpd(1- (1 - u[!ind]) / (1 - qu), sigma = sig, xi = xi, u = 0) + threshold
    }
    return(xx)
}

## calcualte Z matrix
getZMatrix <- function(y.cond, x, para){
    a <- para[1]
    b <- para[2]
    cee <- para[3]
    d <- para[4]
    if (is.na(a))
        rep(NA, length(x)) else {
            if (a < 10^(-5) & b < 0)
                a <- cee - d * log(y.cond) else a <- a * y.cond
            (x - a)/(y.cond^b)
        }
}



########### the folloiwng fnctions are copied from texemx pacgage ###
## a wrapper for NLH optimisation
## returns dependnecne parameters(alpha, beta) and likelihood
my.qfun <- function(y.cond, y.dep, aLow = -1 + 10^(-10), margins = "laplace", constrain = FALSE, v = 10, maxit = 10000, start = c(0.01, 0.01) , nOptim = 1) {
    o <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit), yex = y.cond, ydep = y.dep, constrain = constrain, v = v, aLow = aLow), silent = FALSE)
    ## optim catch block.
    if (class(o) == "try-error") {
        warning("Error in optim call from mexDependence")
        o <- as.list(o)
        o$par <- rep(NA, 6)
        o$value <- NA
    } else if (o$convergence != 0) {
        warning("Non-convergence in mexDependence")
        o <- as.list(o)
        o$par <- rep(NA, 6)
    } else if (nOptim > 1) { ## run optim k time.s
        for (i in 2:nOptim) {
            o <- try(optim(par = o$par, fn = Qpos, control = list(maxit = maxit), yex = y.cond, ydep = y.dep, constrain = constrain, v = v,
                           aLow = aLow), silent = TRUE)
            if (class(o) == "try-error") {
                warning("Error in optim call from mexDependence")
                o <- as.list(o)
                o$par <- rep(NA, 6)
                o$value <- NA
                (break)()
            } else if (o$convergence != 0) {
                warning("Non-convergence in mexDependence")
                o <- as.list(o)
                o$par <- rep(NA, 6)
                (break)()
            }
        }
    }
    ## gumbel margins and negative dependence
    if (!is.na(o$par[1])) {
                                        # gumbel margins and negative dependence
        if (margins == "gumbel" & o$par[1] <= 10^(-5) & o$par[2] < 0) {
            lo <- c(10^(-10), -Inf, -Inf, 10^(-10), -Inf, 10^(-10))
            Qneg <- function(yex, ydep, param) {
                param <- param[-1]
                b <- param[1]
                cee <- param[2]
                d <- param[3]
                m <- param[4]
                s <- param[5]
                obj <- function(yex, ydep, b, cee, d, m, s) {
                    mu <- cee - d * log(yex) + m * yex^b
                    sig <- s * yex^b
                    log(sig) + 0.5 * ((ydep - mu)/sig)^2
                }
                res <- sum(obj(yex, ydep, b, cee, d, m, s))
                res
            }
            o <- try(optim(c(0, 0, 0, 0, 0, 1), Qneg, method = "L-BFGS-B", lower = lo, upper = c(1, 1 - 10^(-10), Inf, 1 - 10^(-10), Inf, Inf),
                           yex = y.cond, ydep = y.dep), silent = TRUE)
            if (class(o) == "try-error" || o$convergence != 0) {
                warning("Non-convergence in mexDependence")
                o <- as.list(o)
                o$par <- rep(NA, 6)
            }
        } else {
                                        # end if gumbel margins and neg dependence
            Z <- (y.dep - y.cond * o$par[1])/(y.cond^o$par[2])
            o$par <- c(o$par[1:2], 0, 0, mean(Z), sd(Z))
        }
    }
    return(c(o$par[1:6], o$value))  # Parameters and negative loglik
}  # Close qfun <- function(

## a wrapper for lower level functions.
## return negative log likeliihohod
Qpos <- function(param, yex, ydep, constrain, v, aLow) {
    a <- param[1]
    b <- param[2]
    res <- PosGumb.Laplace.negProfileLogLik(yex, ydep, a, b, constrain, v, aLow)  # defined in file mexDependenceLowLevelFunctions
    res$profLik
}



                                        # check constraints on parameters under constrained Laplace estimation
ConstraintsAreSatisfied <- function(a,b,z,zpos,zneg,v){
    C1e <- a <= min(1, 1 - b*min(z)*v^(b-1), 1 - v^(b-1)*min(z) + min(zpos)/v) &
        a <= min(1, 1 - b*max(z)*v^(b-1), 1 - v^(b-1)*max(z) + max(zpos)/v)

    C1o <- a <= 1 &
        a > 1 - b * min(z) * v^(b-1) &
            a > 1 - b * max(z) * v^(b-1) &
                (1 - 1/b)*(b*min(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + min(zpos) > 0 &
                    (1 - 1/b)*(b*max(z))^(1/(1-b)) * (1-a)^(-b/(1 - b)) + max(zpos) > 0

    C2e <- -a <= min(1, 1 + b*v^(b-1)*min(z), 1 + v^(b-1)*min(z) - min(zneg)/v) &
        -a <= min(1, 1 + b*v^(b-1)*max(z), 1 + v^(b-1)*max(z) - max(zneg)/v)

    C2o <- -a <= 1 &
        -a > 1 + b*v^(b-1)*min(z) &
            -a > 1 + b*v^(b-1)*max(z) &
                (1-1/b)*(-b*min(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - min(zneg) > 0 &
                    (1-1/b)*(-b*max(z))^(1/(1-b))*(1+a)^(-b/(1-b)) - max(zneg) > 0

    if (any(is.na(c(C1e, C1o, C2e, C2o)))) {
        warning("Strayed into impossible area of parameter space")
        C1e <- C1o <- C2e <- C2o <- FALSE
    }

    (C1e | C1o) && (C2e | C2o)
}

                                        # positive dependence Gumbel and pos or neg dependence Laplace neg likelihood function
PosGumb.Laplace.negloglik <- function(yex, ydep, a, b, m, s, constrain, v, aLow) {
    BigNumber <- 10^40
    WeeNumber <- 10^(-10)

    if(a < aLow[1] | s < WeeNumber | a > 1-WeeNumber  | b > 1-WeeNumber) {
        res <- BigNumber
    } else {
        mu <- a * yex + m * yex^b
        sig <- s * yex^b

        res <- sum(0.5 * log(2*pi) + log(sig) + 0.5 * ((ydep - mu)/sig)^2)

        if (is.infinite(res)){
            if (res < 0){
                res <- -BigNumber
            } else {
                res <- BigNumber
            }
            warning("Infinite value of Q in mexDependence")
        } else if (constrain){
                                        #v <- v * max(yex)
            zpos <- range(ydep - yex) # q0 & q1
            z <- range((ydep - yex * a) / (yex^b)) # q0 & q1
            zneg <- range(ydep + yex) # q0 & q1

            if (!ConstraintsAreSatisfied(a,b,z,zpos,zneg,v)){
                res <- BigNumber
            }
        }
    }
    res
}

PosGumb.Laplace.negProfileLogLik <- function(yex, ydep, a, b, constrain, v, aLow) {
    Z <- (ydep - yex * a) / (yex^b)

    m <- mean(Z)
    s <- sd(Z)

    res <- PosGumb.Laplace.negloglik(yex,ydep,a,b,m=m,s=s,constrain,v,aLow=aLow)
    res <- list(profLik=res,m=m, s=s)
    res
}
