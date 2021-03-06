### automated approchh in slelecting GPD fitting thedholds.

library(texmex)
ds <- as.data.frame(nhi.x[, 2, with=F])
n <- 100
th <- seq(0.8, 0.99, len = n)
ma <- vector("list", len = n)
for (i in seq_len(n)) {
    cat("\t", i)
    ma[[i]] <- migpd(ds, mqu = th[i])
}

ma.coef <- sapply(ma, coef)
ma.coef <- as.data.table(t(ma.coef))
setnames(ma.coef, names(ma.coef), c("th", "qu", "sigma", "xi", "ub"))

par(mfrow=c(2,2))
ma.coef[,
        {
            plot(th, sigma)
            plot(th, xi)
        }]



x11()
par(mfrow=c(1,2))
plot(res <- gpdRangeFit(ds[,1], umin = quantile(ds[,1], 0.8), umax = quantile(ds[, 1], 0.999)))

par(mfrow=c(1, 2))
fi <- gpd.fitrange(ds[, 1], umin = quantile(ds[, 1], 0.8), umax =quantile(ds[,1], 0.999), nint = 100)

x11()
par(mfrow=c(1,2))
plot(fi$thresholds, log(fi$mle[, 1]))
plot(fi$thresholds, fi$mle[, 2])

x <- nhi.x[[2]]
## phi and xi v.s. th
res <- gpdRangeFit(x, umin = quantile(x, 0.8), umax = quantile(x, 0.999), nint = 30)
plot(res)
phi <- data.table(th = res$th,
                  mle = res$par[, 1],
                  lb = res$lo[, 1],
                  ub = res$hi[, 1])
xi <-  data.table(th = res$th,
                  mle = res$par[, 2],
                  lb = res$lo[, 2],
                  ub = res$hi[, 2])
par <- data.table(rbind(phi, xi), par = rep(c("phi", "xi"), c(nrow(phi), nrow(xi))))
## ggplot(par, aes(x=th, y= mle)) + geom_point() + geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.5) + facet_wrap(~par, scale = "free")
x11()
ggplot(par, aes(x=th, y= mle)) + geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub), alpha = 0.5) + facet_wrap(~par, scale = "free")


### mean excendence life (mrl) plot,
## Threshold choice for the fitting of the GPD is guided by the shape of the Mean Residual Life plot.
## A threshold which is suitably high will have a corresponding mrl plot which is approximately linear
## in shape above the threshold (up to sampling variation). 

par(mfrow=c(1, 2))
plot(texmex::mrl(x), umin = quantile(x, 0.8), umax = quantile(x, 0.999))
ismev::mrl.plot(x, umin = quantile(x, 0.8), umax = quantile(x, 0.999))

pdf("1000select.gpd.th.pdf", width = 12, height = 7)
n <- ncol(nhi.x)
for (i in seq_len(n - 1)){
            cat(i/(n-1))
##    for (i in sample(1:(n-1), 20)) {

    x <- nhi.x[[i + 1 ]]
mrl <- data.table(texmex::mrl(x, umin = quantile(x, 0.8), umax = quantile(x, 0.999))$mrl)
## ggplot(mrl, aes(threshold, MRL)) + geom_point() + geom_errorbar(aes(ymin = lo, ymax = hi)) + geom_smooth()
p1 <- ggplot(mrl[complete.cases(mrl), ], aes(threshold, MRL)) + geom_point() + geom_errorbar(aes(ymin = lo, ymax = hi)) + geom_smooth()

l <- loess(MRL ~ threshold, mrl)
e <- c(NA, diff(l$fitted))
mrl[, e := e]
p2 <- ggplot(mrl, aes(threshold, e)) + geom_point() + geom_abline(intercept = 0, slope = 0)
p <- arrangeGrob(p1, p2, ncol = 2 )
    print(p)
}
dev.off()





### choose the optimal threarhold for gpd fitting
### compare the results with benchmark (mqu = 0.995)

## n <- ncol(nhi.x)
## for (i in seq_len(n - 1)){
samp.cells <- sample(names(nhi.x)[-1], 20)
col.ind <- match(samp.cells, names(nhi.x))
pdf("../figure/select.gpd.th.pdf", width = 12, height = 7)
for (i in seq_along(col.ind)) {
    cat("\t", i/20)
    x <- nhi.x[[col.ind[i] ]]

    res <- gpdRangeFit(x, umin = quantile(x, 0.8), umax = quantile(x, 0.995), nint = 100)
    ## plot(res)
    phi <- data.table(th = res$th,
                      mle = res$par[, 1],
                      lb = res$lo[, 1],
                      ub = res$hi[, 1])
    xi <-  data.table(th = res$th,
                      mle = res$par[, 2],
                      lb = res$lo[, 2],
                      ub = res$hi[, 2])
    par <- data.table(rbind(phi, xi), par = rep(c("phi", "xi"), c(nrow(phi), nrow(xi))))
##     ggplot(par, aes(x=th, y= mle)) + geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub), alpha = 0.5) + facet_wrap(~par, scale = "free")

    xi[, n.id := 1:nrow(xi)]
    xi[, n.cover := {
        sum(lb <= xi$mle & ub >= xi$mle)
    }, by = n.id]
    opt.th <- xi[which.max(n.cover), th]

    p <- ggplot(par, aes(x=th, y= mle)) + geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub), alpha = 0.5) + facet_wrap(~par, scale = "free") + geom_errorbar(data = par[th == opt.th], aes(ymin = lb, ymax = ub), col = 2) + labs(title = samp.cells[i])
    print(p)
    ma <- migpd(as.data.frame(x), mqu = 0.995)
    ma.opt <- migpd(as.data.frame(x), mth = opt.th)

    par(mfcol = c(2, 4))
    plot(ma)
    plot(ma.opt)

}
dev.off()



## n <- ncol(nhi.x)
## for (i in seq_len(n - 1)){
## pdf("../figure/select.gpd.th.pdf", width = 12, height = 7)
n <- ncol(nhi.x) - 1
n <- nrow(LT.model)
## gev.para <- vector("list", len = n)
gpd.para <- vector("list", len = n)
for (i in seq_len(n)) {
##     for (i in 1:10){
    cat("\t", i/n)
    i.id <- LT.model[i, id]
    x <- nhi.x[[as.character(i.id)]]
    res <- gpdRangeFit(x, umin = quantile(x, 0.8), umax = quantile(x, 0.995), nint = 10)
    ## plot(res)
    phi <- data.table(th = res$th,
                      mle = res$par[, 1],
                      lb = res$lo[, 1],
                      ub = res$hi[, 1])
    xi <-  data.table(th = res$th,
                      mle = res$par[, 2],
                      lb = res$lo[, 2],
                      ub = res$hi[, 2])
    par <- data.table(rbind(phi, xi), par = rep(c("phi", "xi"), c(nrow(phi), nrow(xi))))
    xi[, n.id := 1:nrow(xi)]
    xi[, n.cover := {
        sum(lb <= xi$mle & ub >= xi$mle)
    }, by = n.id]
    opt.th <- xi[which.max(n.cover), th]
    opt.sigma <- exp(phi[th == opt.th, mle])
    opt.xi <- xi[opt.th, mle]
    ei <- extremalIndex(x, threshold = opt.th)
    dc <- declust(ei, r = 4) ## trigger.window
    lambda <- dc$nClusters / 34

    (p <- ggplot(par, aes(x=th, y= mle)) + geom_point() + geom_errorbar(aes(ymin = lb, ymax = ub), alpha = 0.5) + facet_wrap(~par, scale = "free") + geom_errorbar(data = par[th == opt.th], aes(ymin = lb, ymax = ub), col = 2) + labs(title = samp.cells[i]))

    gpd.para <- data.table(id = i.id, th = opt.th, sigma = opt.sigma, xi = opt.xi, lambda = lambda)
   ##  para <- GPD2GEV(mu = opt.th, sigma = opt.sigma, xi = opt.xi, lambda = lambda)
##     names(para) <- c("mu", "sigma", "xi")
##     gev.para[[i]] <- data.table(id = i.id, t(para))
}

gev.para <- gpd.para[, GPD2GEV(mu = th, sigma, xi, lambda), by = id]

g <- copy(gev.para)
sapply(c(10, 34, 100, 1000, 10000), function(i)
       {
           g[, paste0("rl", i) := qgev(1 - 1/i, mu, sigma, xi), by = id]
       })

keep.cols <- names(g)[grep("rl", names(g), fixed = TRUE)]
keep.cols <- c("id", keep.cols)
gg <- g[, keep.cols, with=F]
gg <- melt(gg, id.vars = "id", variable.name = "rl", value.name = "precip", variable.factor = FALSE)

gg <- merge(gg, nhi.info, by = "id")
p <- ggplot(gg, aes(lon, lat)) + geom_tile(aes(fill = round(precip))) + facet_wrap(~rl)


ggplot(gg[rl=="rl100" & precip <= 800], aes(lon, lat)) + geom_tile(aes(fill = round(precip))) + facet_wrap(~rl) + scale_fill_gradientn(colours = rev(rainbow(7, start = 0, end = 0.65)))









x <- rgpd(1e6, sigma = 10, xi = 0.2)
res <- gpdRangeFit(x, umin = quantile(x, 0.8), umax = quantile(x, 0.9995))
    phi <- data.table(th = res$th,
                      mle = res$par[, 1],
                      lb = res$lo[, 1],
                      ub = res$hi[, 1])
    xi <-  data.table(th = res$th,
                      mle = res$par[, 2],
                      lb = res$lo[, 2],
                      ub = res$hi[, 2])


i.th <- res$th[9]
coef(migpd(as.data.frame(x), mth = i.th))






th <- seq(quantile(x, 0.8), quantile(x, 0.99), len = 10)
coef <- list()
for (i in 1:10){
    i.th <- th[i]
    coef[[i]] <- coef(migpd(as.data.frame(x), mth = i.th))
}


s <- sapply(coef, function(x) x[3])
