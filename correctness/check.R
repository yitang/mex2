
{
this.thread <- 1
setwd("\\\\SKI-CQV5LZ1\\vmshare")
source(".RProfile")
if (!exists("nearby.all"))
    load(link$cell.info.nearby, v=T)
i <- this.thread
utils::setWindowTitle(paste0("evo.",i))
source("0_methods.R")
source(link$my.texmex)

## file <- paste0("keydata/table", this.thread, ".triggers.RData") ## all.triggers
## load(file, v=TRUE)
load("keydata/block1.triggers.RData", v = TRUE)
setkey(all.triggers, id, s.ymd)

## distribute works to this thread
tmp <- all.triggers[, length(s.ymd), by = id][order(V1, decreasing = TRUE),] ## distrbute ids evenly to threads according number of simulation per id.
this.thread.ids <- tmp[seq(this.thread, nrow(tmp), by = 10), ]$id
all.triggers <- all.triggers[J(this.thread.ids), ]
setkey(all.triggers, id, s.ymd)

## simulation
unique.triggers <- unique(all.triggers$id)
all.evo.samp <- vector("list", len = length(unique.triggers))
names(all.evo.samp) <- unique.triggers
count <- 0
t0 <- proc.time()
for (i in unique.triggers){
    ##for (i in sample(unique.triggers, 100)){
    count <-  count + 1
    cat("\t", count / length(unique.triggers))
    load(paste0("info/EachCell_old/", i, ".RData"))# LT.model, HT.model
    this.trigger.info <- all.triggers[J(i), ]
    this.evo <- vector("list", len = nrow(this.trigger.info))
    for (j in seq_len(nrow(this.trigger.info))){
        this.day <- this.trigger.info[j, s.ymd]
        this.day.rain <- this.trigger.info[j, rain]
        HT.model$y.threshold <- this.day.rain
        ## rejection sampling: y0 is the largest in day 0.
        day0 <- grep(".0", names(this.evo.samp), fixed = TRUE)
        apply(evo.samp, 1, function(x)
              qLaplace(x[1]) > max(x[day0]))
        evo.samp <- samp.cond(HT.model)
        this.evo.samp <- evo.samp[sample(nrow(evo.samp), 1), ]

        ## this.evo.samp <- cbind(this.trigger.infos[i, s.ymd], this.trigger.info, this.evo.samp)
        this.evo[[j]] <- cbind(s.ymd = this.day, id = i, this.evo.samp)
    }
    this.evo <- rbindlist(this.evo) ## cast to data.table
    all.evo.samp[[paste(i)]] <- this.evo ## match by name
}
t1 <- proc.time()
cat("\n Takes about ", (t1-t0)[3], " seconds to run evo simulation for all triggers")
cat("\n memoery usage: \t")
print(pryr::mem_used())
save(all.evo.samp, file = paste0("keydata/block1.Evo.samp_", this.thread, ".Rdata"))
}




library(texmex)
load("info/NearbyCells_All.RData",v=T)
load("nhi.x.only.RData")
source("texmex_hack/mex2.R")

i = 37960
nearby <- nearby.all[J(i),]$nearby
nearby <- setdiff(nearby, i)
dt <- nhi.x[, as.character(c(i, nearby)), with=F]
r2 <- mex(as.data.frame(dt), mqu = 0.995, dqu = 0.995, which = 1, constrain = TRUE, penalty = "none")
coef2 <- coef(r2$dependence)[c(1, 2), ]
coef2 <- data.table(id = colnames(coef2), a = coef2[1, ], b = coef2[2, ])

r2.v5 <- mex(as.data.frame(dt), mqu = 0.995, dqu = 0.995, which = 1, constrain = TRUE, penalty = "none", v=5)
coef2.v5 <- coef(r2.v5$dependence)[c(1, 2), ]
coef2.v5 <- data.table(id = colnames(coef2.v5), a = coef2.v5[1, ], b = coef2.v5[2, ])

r2.qu0.99 <- mex(as.data.frame(dt), mqu = 0.99, dqu = 0.995, which = 1, constrain = TRUE, penalty = "none")
coef2.qu0.99 <- coef(r2.qu0.99$dependence)[c(1, 2), ]
coef2.qu0.99 <- data.table(id = colnames(coef2.qu0.99), a = coef2.qu0.99[1, ], b = coef2.qu0.99[2, ])

r1 <- mex2(dt, dqu = 0.995)
coef1 <- r1$coef
r1.qu0.99 <- mex2(dt, dqu = 0.995)
coef1.qu0.99 <- r1.qu0.99$coef

par(mfrow=c(2,3))
coef2[, plot(a, b, main =  "texmex")]
coef2.v5[, plot(a, b, main = "texmex, v5")]
coef2.qu0.99[, plot(a, b, main = "texmex, qu 0.99")]
coef1[, plot(a, b, main = "hack")]

coef1.qu0.99[, plot(a, b, main = "hack, qu 0.99")]


load(paste0("info/EachCell_old/", i, ".RData"))# LT.model, HT.model
