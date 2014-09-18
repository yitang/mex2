source(".RProfile")
this.id <- "44405"
load(link$nhi,v=T)
load(link$nhi.info)
load("info/EachCell/44405.RData")
nearby <- merge(nearby[, list(id = nearby)], cfsr.cell, by = "id")
nearby <- nearby[lon >= -118&lat == 59.5]
const_qu = 0.995
this.ids <- as.character(nearby$id)

df <- as.data.frame(nhi.x[, this.ids, with=F])
require(texmex)
ht <- mex(df, mqu = const_qu, which = 1, dqu = const_qu, penalty = "none")
z <- as.data.table(ht$dependence$Z)
y <- as.data.table(ht$margins$transformed)
y.big <- y[ y[[1]] > quantile(y[[1]], const_qu), ]
coef <- t(coef(ht$dependence)[1:2, ])
coef <- data.table( id = rownames(coef), coef)
y0 <- y.big[[1]]

## residual plots
zz <- melt(z, variable.name = "id", value.name = "z", variable.factor = F)
zz[, y0 := y0]
p0 <- ggplot() + geom_point(data = zz, aes(y0, z)) + facet_wrap(~id)


## tail distribtion plots
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


## my hacks

###################### formulate data with lag #######################

dt <- as.data.table(df)
dt.laplace2 <- dt[, lapply(.SD, toLaplace)]

aLow <- -1 + 10^(-10)
marTransform <- "mixture"
margins <- "laplace"
start = c(0.01, 0.01)
constrain = TRUE
v= 10
nOptim = 1
maxit <- 1e4


y.cond <- dt.laplace2[[1]]
y.cond.th <- quantile(y.cond, 0.995)
ind <- y.cond >= y.cond.th
y.cond <- y.cond[ind]
dt.laplace2 <- dt.laplace2[ind, ]
dt.dep2 <- dt.laplace2[, lapply(.SD, function(x)
  my.qfun(y.cond, x,
          aLow = aLow, margins = margins, constrain = TRUE, maxit = maxit, start = start, v= v,  nOptim = 1)),
  .SDcols = -1]

tmp <- rbind(dt.dep2[1:4, ], dt.laplace2[, -1, with = F])
z.mat2 <- tmp[, lapply(.SD, function(x)
  getZMatrix(y.cond = y.cond, x[-c(1:4)], para = x[1:4]))]

# cast and compare
y.big <- dt.laplace2
y0 <- y.big[[1]]
coef <- data.table(t(dt.dep2[1:2, ]), id = names(dt.dep2))
setnames(coef, names(coef)[1:2], c("a", "b"))
zz <- melt(z.mat2, variable.name = "id", value.name = "z", variable.factor = F)
zz[, y0 := y0]
p0.hack <- ggplot() + geom_point(data = zz, aes(y0, z)) + facet_wrap(~id)

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
p1.hack <- ggplot() + geom_point(data = yy, aes(y0, y)) + geom_ribbon(data = band2, aes(y.seq, ymin = lower, ymax = upper), alpha = 0.2) + facet_wrap(~id)




###########################################################################
## check my hacks                                                        ##
###########################################################################




## mexDependence2(ma, which = 1, dqu = 0.995)
setwd("texmex_mod")
f <- list.files()
sapply(f, source)
setwd("..")

nearby <- c("34319", "34320", "34321", "34322", "34323", "34324", "34325",
            "34326", "34327", "34328", "34329", "34330", "34331", "35026"
)

###################### formulate data with lag #######################
dt <- as.list(nhi.x[,
                    rep(nearby[-1], 1 + const_R), with=F])
lag.ind <- equal.range( 1:length(dt), const_R + 1)
lag.ind[, lag := block - 1]
lag.ind <- lag.ind[lag != 0, ]
for (i in seq_len(nrow(lag.ind))){
  col.ind <- lag.ind[i, ind]
  lag <- lag.ind[i, lag]
  dt[[col.ind]] <- dt[[col.ind]][-c(1:lag)]
}
dt <- as.data.table(dt)
dt <- cbind(nhi.x[, nearby[1], with = F], dt)## add conditioning varibale
dt <- dt[-c((nrow(dt) - const_R + 1 ) : (nrow(dt))), ]

###############################################
## para <-gpd.fit.coefs[, names(dt), with=F] ##
## dt <- rbind(para, dt)                     ##
###############################################
setnames(dt, 1:ncol(dt), paste0("x", 1:ncol(dt) - 1))

dt.laplace2 <- dt[, lapply(.SD, toLaplace)]

## mexDependence2(migpd.obj, 1, dqu = 0.995)
aLow <- -1 + 10^(-10)
marTransform <- "mixture"
margins <- "laplace"
start = c(0.01, 0.01)
constrain = TRUE
v= 10
nOptim = 1
maxit <- 1e6


y.cond <- dt.laplace2[[1]]
y.cond.th <- quantile(y.cond, 0.995)
ind <- y.cond >= y.cond.th
y.cond <- y.cond[ind]
dt.laplace2 <- dt.laplace2[ind, ]
dt.dep2 <- dt.laplace2[, lapply(.SD, function(x)
  my.qfun(y.cond, x,
          aLow = aLow, margins = margins, constrain = TRUE, maxit = 1e6, start = start, v= v, nOptim = nOptim)),
  .SDcols = -1]






dt.laplace <- dt[, lapply(.SD, function(x)
  mexTransform2(x[-c(1:2)],
                qu = 0.995,
                gpd.coef=x[1:2]))]

y.cond <- dt.laplace[[1]]
y.cond.th <- quantile(y.cond, 0.995)
ind <- y.cond >= y.cond.th
y.cond <- y.cond[ind]
dt.laplace <- dt.laplace[ind, ]




dt.dep <- dt.laplace[, lapply(.SD, function(x)
  my.qfun(y.cond, x,
          aLow = aLow, margins = margins, constrain = TRUE, maxit = 1e6, start = start, v= v)),
  .SDcols = -1]


tmp <- rbind(dt.dep2[1:4, ], dt.laplace2[, -1, with = F])
z.mat2 <- tmp[, lapply(.SD, function(x)
  getZMatrix(y.cond = y.cond, x[-c(1:4)], para = x[1:4]))]
tmp <- rbind(dt.dep2[1:2, ], z.mat2)
samp2 <- tmp[, lapply(.SD, function(x)
  x[-c(1:2)] * (y.cond.th^x[2]) + (y.cond.th*x[1]))]




test.z <- cbind(z.mat2, res2$dependence$Z)

ma2 <- migpd(as.data.frame(dt), mqu = 0.995, penalty = "none")
ma2$models <- as.data.table(lapply(ma2$models, function(i) i$coefficients))
res2 <- mexDependence2(ma2, which = 1, dqu = 0.995)

## dependence parametes are
cbind(dt.dep2, res2$dependence$coefficients, res.texmex$dependence$coefficients)[1:6, ]
x1      V1      x1
1: 0.20474 0.20474 0.23646
2: 0.87566 0.87566 0.87121
3: 0.00000 0.00000 0.00000
4: 0.00000 0.00000 0.00000
5: 0.96615 0.96615 0.93400
6: 0.09681 0.09681 0.09751



setwd(proj.folder)
source("refile.R")
load("info/GPD_Paras_AllCell_20140812_1224.RData", v=T)
load("info/EachCell/100449.RData", v=T )
nearby <- as.character(nearby$nearby)
dt <- as.list(nhi.x[,
                    rep(nearby[-1], 1 + const_R), with=F])
lag.ind <- equal.range( 1:length(dt), const_R + 1)
lag.ind[, lag := block - 1]
lag.ind <- lag.ind[lag != 0, ]
for (i in seq_len(nrow(lag.ind))){
  col.ind <- lag.ind[i, ind]
  lag <- lag.ind[i, lag]
  dt[[col.ind]] <- dt[[col.ind]][-c(1:lag)]
}
dt <- as.data.table(dt)
dt <- cbind(nhi.x[, nearby[1], with = F], dt)## add conditioning varibale
dt <- dt[-c((nrow(dt) - const_R + 1 ) : (nrow(dt))), ]
dt.org <- copy(dt)

## my hacks
dt <- copy(dt.org)
para <-gpd.fit.coefs[, names(dt), with=F]
dt <- rbind(para, dt)
setnames(dt, 1:ncol(dt), paste0("x", 1:ncol(dt) - 1))

dt.laplace2 <- dt[, lapply(.SD, function(x)
  toUnifExtro_laplace(x[-c(1:2)],
                      qu = 0.995,
                      coef=x[1:2]))]
dt.laplace2[dt.laplace2[[1]] >= quantile(dt.laplace2[[1]], 0.995), ]
aLow <- -1 + 10^(-10)
marTransform <- "mixture"
margins <- "laplace"
start = c(0.01, 0.01)
constrain = TRUE
v= 10
nOptim = 1
maxit <- 1e4


y.cond <- dt.laplace2[[1]]
y.cond.th <- quantile(y.cond, 0.995)
ind <- y.cond >= y.cond.th
y.cond <- y.cond[ind]
dt.laplace2 <- dt.laplace2[ind, ]
dt.dep2 <- dt.laplace2[, lapply(.SD, function(x)
  my.qfun(y.cond, x,
          aLow = aLow, margins = margins, constrain = TRUE, maxit = maxit, start = start, v= v,  nOptim = 1)),
  .SDcols = -1]

tmp <- rbind(dt.dep2[1:4, ], dt.laplace2[, -1, with = F])
z.mat2 <- tmp[, lapply(.SD, function(x)
  getZMatrix(y.cond = y.cond, x[-c(1:4)], para = x[1:4]))]
tmp <- rbind(dt.dep2[1:2, ], z.mat2)
samp2 <- tmp[, lapply(.SD, function(x)
  x[-c(1:2)] * (y.cond.th^x[2]) + (y.cond.th*x[1]))]


## use texmex package
dt <- copy(dt.org)
res.texmex <- texmex::mex(as.data.frame(dt), which = 1, mqu = 0.995, dqu = 0.995, maxit = maxit, penalty = "none")
## dt.laplace.texmex <- as.data.table(res.texmex$margins$transformed)
## dt.laplace.texmex[dt.laplace.texmex[[1]] >= quantile(dt.laplace.texmex[[1]], 0.995), ]
tmp <- res.texmex$dependence$coefficients[1:2, ]
tmp <- as.data.table(rbind(tmp, res.texmex$dependence$Z))
samp.texmex <- tmp[, lapply(.SD, function(x)
  x[-c(1:2)] * (y.cond.th^x[2]) + (y.cond.th*x[1]))]

## check the difference bewenee texemx package and my hacks
diff <- unlist(samp2 - samp.texmex)
hist(diff)
