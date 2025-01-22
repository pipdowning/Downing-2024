# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

# R code for:

				### MICHENER'S GROUP-SIZE PARADOX IN COOPERATIVELY BREEDING BIRDS ###

# Author: Philip A. Downing
# Contact: philip.downing@oulu.fi

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

## Packages ##

library(ape)
library(MCMCglmm)
library(metafor)
library(doBy)
library(lattice)
library(nlme)

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

### DATA MANIPULATION ###

# reproductive success data
gsDataLong <- read.csv(".../all_data.csv")
head(gsDataLong)
tail(gsDataLong)
str(gsDataLong)						# 136 fledlging success estimates in different group sizes
range(gsDataLong$groupSize)			# 2 to 17
length(unique(gsDataLong$animal))	# 26 species

# group size frequency data
groupSizes <- read.csv(".../group_frequencies.csv")
groupSizes <- groupSizes[-which(is.na(groupSizes$percGrps)),]
head(groupSizes)
tail(groupSizes)
str(groupSizes)						# 163 estimates of different sized groups
range(groupSizes$groupSize)			# 2 to 17
length(unique(groupSizes$animal))	# 23 species

# zscore total reproductive success (t_rs) for each species separately
t_rsValues <- tapply(gsDataLong$mean, list(gsDataLong$common), scale)
str(t_rsValues)
gsDataLong$t_rs <- 0
speciesNames <- unique(gsDataLong$common)
for(i in 1:length(speciesNames)){
gsDataLong$t_rs[which(gsDataLong$common ==speciesNames[i])] <- t_rsValues[[which(names(t_rsValues) == speciesNames[[i]])]]}
# check it worked
scale(gsDataLong$mean[which(gsDataLong$common == "azure-winged magpie")])
gsDataLong$t_rs[which(gsDataLong$common == "azure-winged magpie")]

# calculate and zscore per capita reproductive success (pc_rs) for each species separately
gsDataLong$pc_rs <- gsDataLong$mean / gsDataLong$groupSize
pc_rsValues <- tapply(gsDataLong$pc_rs, list(gsDataLong$common), scale)
str(pc_rsValues)
gsDataLong$pc_rs <- 0
for(i in 1:length(speciesNames)){
gsDataLong$pc_rs[which(gsDataLong$common ==speciesNames[i])] <- pc_rsValues[[which(names(pc_rsValues) == speciesNames[[i]])]]}
# check it worked
scale(gsDataLong$pc_rs[which(gsDataLong$common == "azure-winged magpie")])
gsDataLong$pc_rs[which(gsDataLong$common == "azure-winged magpie")]

# add log transformed group size (ln_gs) and percent groups (ln_pg) to the data
gsDataLong$ln_gs <- log(gsDataLong$groupSize)
gsDataLong$ln_pg <- log(gsDataLong$perG)
gsDataLong$ln_N <- log(gsDataLong$N)

## phylogenetic trees ##
trees <- read.nexus(".../BirdZilla1300.nex")
is.ultrametric(trees[[1]])

## replace Placid Greenbul with a close relative to match with the tree ##
gsDataLong$animal[which(gsDataLong$animal == "Phyllastrephus_placidus")] <- "Phyllastrephus_cabanisi"

## trim the trees ##
dropTip <- trees[[1]]$tip.label[which(is.na(match(trees[[1]]$tip.label, gsDataLong$animal)))]
gsTrees <- lapply(trees, drop.tip, dropTip, trim.internal=T)
gsTrees <- lapply(gsTrees, makeNodeLabel, method = "number")
plot(gsTrees[[1]])
length(gsTrees[[1]]$tip.label)		# reduced to 25 species
gsDataLong$animal[which((gsDataLong$animal %in% gsTrees[[1]]$tip.label) == FALSE)]
gsTrees[[1]]$tip.label[which((gsTrees[[1]]$tip.label %in% gsDataLong$animal) == FALSE)]

## variance-covariance matrices ##
sample(1:1300, 1)
birdCorA <- vcv.phylo(gsTrees[[604]], cor=TRUE)
# delete the four species missing data on the frequency of different group sizes
drop <- unique(gsDataLong$animal[which(is.na(gsDataLong$ln_pg))])
gsTreesB <- lapply(gsTrees, drop.tip, drop, trim.internal=T)
gsTreesB <- lapply(gsTreesB, makeNodeLabel, method = "number")
birdCorB <- vcv.phylo(gsTreesB[[604]], cor=TRUE)
# note that because sample at random for the 'birdCor' VCV
# metafor models might give slightly different parameter estimates to those reported below

## priors ##

priorA <- list(R=list(V=1,nu=0.002), G=list(G1=list(V=diag(2),nu=0.002), G2=list(V=1,nu=0.002)))


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

				### Part A: per capita and total reproductive success across group sizes ###

## Exploratory Plots ##

# log group size
hist(gsDataLong$ln_gs)			# normal(ish)

# standardised per capita reproductive success #
hist(gsDataLong$pc_rs)			# normal(ish)
plot(pc_rs ~ ln_gs, gsDataLong, pch = 21, bg = factor(common))
xyplot(pc_rs ~ ln_gs | common, data=gsDataLong)

# standardised total reproductive success #
hist(gsDataLong$t_rs)			# normal(ish)
plot(t_rs ~ ln_gs, gsDataLong, pch = 21, bg = factor(common))
xyplot(t_rs ~ ln_gs | common, data=gsDataLong)


## Random Slope Models (MCMCglmm) ##
# see Hadfield, J. (2021) MCMCmmm Course Notes

# standardised per capita reproductive success #

INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
mod.start <- MCMCglmm(pc_rs ~ ln_gs, random= ~ us(1 + ln_gs):common + animal, data = gsDataLong, nitt=1000, thin=1, burnin=0, pr=TRUE, pl=TRUE, prior=priorA, ginverse=list(animal=INtree$Ainv), family ="gaussian", verbose=FALSE)
mod <- mod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=mod$Liab[1,], R=mod$VCV[1,6], G=list(G1=matrix(ncol=2, nrow=2, mod$VCV[1,1:4]), G2=mod$VCV[1,5]))
  mod <- MCMCglmm(pc_rs ~ ln_gs, random= ~ us(1 + ln_gs):common + animal, data = gsDataLong, ginverse=list(animal=INtree$Ainv), family ="gaussian", prior=priorA, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=999, start=start, verbose=FALSE)
  if(i > 300){
    mod.start$VCV[i-300,] <- mod$VCV[1,]
    mod.start$Sol[i-300,] <- mod$Sol[1,]
    mod.start$Liab[i-300,] <- mod$Liab[1,]
  }
}
rs_pc_rs_modA.1 <- mod.start
rs_pc_rs_modA.2 <- mod.start
rs_pc_rs_modA.3 <- mod.start
save(rs_pc_rs_modA.1, file=".../rs_pc_rs_modA.1")
save(rs_pc_rs_modA.2, file=".../rs_pc_rs_modA.2")
save(rs_pc_rs_modA.3, file=".../rs_pc_rs_modA.3")

# chain convergence
hist(rs_pc_rs_modA.1$Liab)
plot(rs_pc_rs_modA.1$VCV)     		# parameters well mixed
plot(rs_pc_rs_modA.1$Sol[,1:2])	  	# intercept and slope well mixed
autocorr(rs_pc_rs_modA.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(rs_pc_rs_modA.1$Sol[,1:2]) 	# correlation between successive samples < 0.1 for all components
rs_pc_rs_modA.Sols <- mcmc.list(list(rs_pc_rs_modA.1$Sol[,1:2], rs_pc_rs_modA.2$Sol[,1:2], rs_pc_rs_modA.3$Sol[,1:2]))
plot(rs_pc_rs_modA.Sols)
gelman.diag(rs_pc_rs_modA.Sols)    	 	# upper CI = 1.00 suggesting convergence
heidel.diag(rs_pc_rs_modA.1$VCV)    		# all parameters passed halfwidth
heidel.diag(rs_pc_rs_modA.1$Sol[,1:2])  # intercept and slope passed halfwidth

# model parameters
summary(rs_pc_rs_modA.1)
posterior.mode(rs_pc_rs_modA.1$Sol[,1:2])		# intercept = 1.16, slope = -0.65
HPDinterval(rs_pc_rs_modA.1$Sol[,1:2])			# intercept = 0.29 to 1.72, slope = -1.01 to -0.13

# intercept slope correlation
int.slope.cor <- rs_pc_rs_modA.1$VCV[,2]/sqrt(rs_pc_rs_modA.1$VCV[,1] * rs_pc_rs_modA.1$VCV[,4])
posterior.mode(int.slope.cor)
autocorr(int.slope.cor)


# standardised total reproductive success #

INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
mod.start <- MCMCglmm(t_rs ~ ln_gs, random= ~ us(1 + ln_gs):common + animal, data = gsDataLong, nitt=1000, thin=1, burnin=0, pr=TRUE, pl=TRUE, prior=priorA, ginverse=list(animal=INtree$Ainv), family ="gaussian", verbose=FALSE)
mod <- mod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=mod$Liab[1,], R=mod$VCV[1,6], G=list(G1=matrix(ncol=2, nrow=2, mod$VCV[1,1:4]), G2=mod$VCV[1,5]))
  mod <- MCMCglmm(t_rs ~ ln_gs, random= ~ us(1 + ln_gs):common + animal, data = gsDataLong, ginverse=list(animal=INtree$Ainv), family ="gaussian", prior=priorA, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=999, start=start, verbose=FALSE)
  if(i > 300){
    mod.start$VCV[i-300,] <- mod$VCV[1,]
    mod.start$Sol[i-300,] <- mod$Sol[1,]
    mod.start$Liab[i-300,] <- mod$Liab[1,]
  }
}
rs_t_rs_modA.1 <- mod.start
rs_t_rs_modA.2 <- mod.start
rs_t_rs_modA.3 <- mod.start
save(rs_t_rs_modA.1, file=".../rs_t_rs_modA.1")
save(rs_t_rs_modA.2, file=".../rs_t_rs_modA.2")
save(rs_t_rs_modA.3, file=".../rs_t_rs_modA.3")

# chain convergence
hist(rs_t_rs_modA.1$Liab)
plot(rs_t_rs_modA.1$VCV)     		# parameters well mixed
plot(rs_t_rs_modA.1$Sol[,1:2])	  	# intercept and slope well mixed
autocorr(rs_t_rs_modA.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(rs_t_rs_modA.1$Sol[,1:2]) 	# correlation between successive samples < 0.1 for all components
rs_t_rs_modA.Sols <- mcmc.list(list(rs_t_rs_modA.1$Sol[,1:2], rs_t_rs_modA.2$Sol[,1:2], rs_t_rs_modA.3$Sol[,1:2]))
plot(rs_t_rs_modA.Sols)
gelman.diag(rs_t_rs_modA.Sols)     		# upper CI = 1.00 suggesting convergence
heidel.diag(rs_t_rs_modA.1$VCV)    		# all parameters passed halfwidth
heidel.diag(rs_t_rs_modA.1$Sol[,1:2])   # intercept and slope passed halfwidth

# model parameters
summary(rs_t_rs_modA.1)
posterior.mode(rs_t_rs_modA.1$Sol[,1:2])		# intercept = -1.64, slope = 1.10
HPDinterval(rs_t_rs_modA.1$Sol[,1:2])		# intercept = -2.11 to -1.12, slope = 0.80 to 1.53

# intercept slope correlation
int.slope.cor <- rs_t_rs_modA.1$VCV[,2]/sqrt(rs_t_rs_modA.1$VCV[,1] * rs_t_rs_modA.1$VCV[,4])
posterior.mode(int.slope.cor)
autocorr(int.slope.cor)


## Two-Stage Models (metafor) ##

# see Viechtbauer, W. (2010) https://www.metafor-project.org/doku.php/tips:two_stage_analysis

# standardised per capita reproductive success #

# get slope and sampling variance estimates
res.List <- lmList(groupedData(pc_rs ~ ln_gs | common, data = gsDataLong[,c("pc_rs", "ln_gs", "common")]))
summary(res.List)		# 17/26 slopes negative
plot(augPred(res.List), grid=TRUE)
b <- lapply(res.List, coef)
V <- lapply(res.List, vcov)

# combine in a dataframe
metData <- data.frame(slope = unlist(lapply(b, "[[", 2)),  var = unlist(lapply(V, "[[", 4)))
metData$common <- rownames(metData)
metData$animal <- gsDataLong$animal[match(metData$common, gsDataLong$common)]

# run an intercept only model
ts_pc_rs_modA <- rma.mv(slope ~ 1, V=var, random = list(~ 1 | common, ~ 1 | animal), R = list(animal=birdCorA), data=metData)
summary(ts_pc_rs_modA)   	# B = -0.88, lwr = -1.77, upr = -0.01	, Q = 184, p < 0.001


# standardised total reproductive success #

# get slope and sampling variance estimates
res.List <- lmList(groupedData(t_rs ~ ln_gs | common, data = gsDataLong[,c("t_rs", "ln_gs", "common")]))
summary(res.List)		# 24/26 slopes positive
plot(augPred(res.List), grid=TRUE)
b <- lapply(res.List, coef)
V <- lapply(res.List, vcov)

# combine in a dataframe
metData <- data.frame(slope = unlist(lapply(b, "[[", 2)),  var = unlist(lapply(V, "[[", 4)))
metData$common <- rownames(metData)
metData$animal <- gsDataLong$animal[match(metData$common, gsDataLong$common)]

# run an intercept only model
ts_t_rs_modA <- rma.mv(slope ~ 1, V=var, random = list(~ 1 | common, ~ 1 | animal), R = list(animal=birdCorA), data=metData)
summary(ts_t_rs_modA)   	# B = 2.09, lwr = 1.79, upr = 2.39, Q = 110, p < 0.0001


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

	### Part B: per capita and total reproductive success and the frequency of group sizes ###

## Group Size Frequencies ##

# total number of groups across species
sum(groupSizes$N, na.rm=TRUE)		# 16101

# summarise data by group size across species
counts <- data.frame(groupSize = 2:17, count = tapply(groupSizes$N, list(groupSizes$groupSize), sum, na.rm=TRUE))
counts
sum(counts$count)		# 16101

# work out percentage of each group size
counts$perc <- (counts$count / sum(counts$count)) * 100
sum(counts$perc)		# 100
sum(counts$perc[1:2])	# pairs and trios comprise 79% of groups
sum(counts$perc[3:4])	# threes and fours comprise 18% of groups
sum(counts$perc[5:16])	# groups with 5 + members comprise 3% of groups


## Exploratory Plots ##

# histogram of group size distributions
hist(rep(counts$groupSize, counts$count), breaks=1:18, ylab="Number of groups", xlab="Group size", axes=F)
axis(side=1, pos=0, las=1, at=1:16+0.5, labels=2:17, cex.axis=1, col="black")
axis(side=2)

# histogram of group sizes within species
xyplot(N ~ groupSize | common, data=groupSizes)
xyplot(percGrps ~ groupSize | common, data=groupSizes)

# log % groups
hist(gsDataLong$ln_pg)			# normal(ish)

## standardised per capita reproductive success ##
plot(pc_rs ~ ln_pg, gsDataLong, pch = 21, bg = factor(common))
xyplot(pc_rs ~ ln_pg | common, data=gsDataLong)

## standardised total reproductive success ##
plot(t_rs ~ ln_pg, gsDataLong, pch = 21, bg = factor(common))
xyplot(t_rs ~ ln_pg | common, data=gsDataLong)


## Random Slope Models (MCMCglmm) ##

# standardised per capita reproductive success #

INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
mod.start <- MCMCglmm(pc_rs ~ ln_pg, random= ~ us(1 + ln_pg):common + animal, data = gsDataLong[-which(is.na(gsDataLong$ln_pg)),], nitt=1000, thin=1, burnin=0, pr=TRUE, pl=TRUE, prior=priorA, ginverse=list(animal=INtree$Ainv), family ="gaussian", verbose=FALSE)
mod <- mod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=mod$Liab[1,], R=mod$VCV[1,6], G=list(G1=matrix(ncol=2, nrow=2, mod$VCV[1,1:4]), G2=mod$VCV[1,5]))
  mod <- MCMCglmm(pc_rs ~ ln_pg, random= ~ us(1 + ln_pg):common + animal, data = gsDataLong[-which(is.na(gsDataLong$ln_pg)),], ginverse=list(animal=INtree$Ainv), family ="gaussian", prior=priorA, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=999, start=start, verbose=FALSE)
  if(i > 300){
    mod.start$VCV[i-300,] <- mod$VCV[1,]
    mod.start$Sol[i-300,] <- mod$Sol[1,]
    mod.start$Liab[i-300,] <- mod$Liab[1,]
  }
}
rs_pc_rs_modB.1 <- mod.start
rs_pc_rs_modB.2 <- mod.start
rs_pc_rs_modB.3 <- mod.start
save(rs_pc_rs_modB.1, file=".../rs_pc_rs_modB.1")
save(rs_pc_rs_modB.2, file=".../rs_pc_rs_modB.2")
save(rs_pc_rs_modB.3, file=".../rs_pc_rs_modB.3")

# chain convergence
hist(rs_pc_rs_modB.1$Liab)
plot(rs_pc_rs_modB.1$VCV)     		# all parameters
plot(rs_pc_rs_modB.1$Sol[,1:2])	  	# intercept and slope well mixed
autocorr(rs_pc_rs_modB.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(rs_pc_rs_modB.1$Sol[,1:2]) 	# correlation between successive samples < 0.1 for all components
rs_pc_rs_modB.Sols <- mcmc.list(list(rs_pc_rs_modB.1$Sol[,1:2], rs_pc_rs_modB.2$Sol[,1:2], rs_pc_rs_modB.3$Sol[,1:2]))
plot(rs_pc_rs_modB.Sols)
gelman.diag(rs_pc_rs_modB.Sols)    	 	# upper CI = 1.00 suggesting convergence
heidel.diag(rs_pc_rs_modB.1$VCV)    		# all parameters passed halfwidth
heidel.diag(rs_pc_rs_modB.1$Sol[,1:2])  # intercept and slope passed halfwidth

# model parameters
summary(rs_pc_rs_modB.1)
posterior.mode(rs_pc_rs_modB.1$Sol[,1:2])		# intercept = -0.63, slope = 0.31
HPDinterval(rs_pc_rs_modB.1$Sol[,1:2])			# intercept = -1.16 to -0.19, slope = 0.09 to 0.46

# intercept slope correlation
int.slope.cor <- rs_pc_rs_modB.1$VCV[,2]/sqrt(rs_pc_rs_modB.1$VCV[,1] * rs_pc_rs_modB.1$VCV[,4])
posterior.mode(int.slope.cor)
autocorr(int.slope.cor)


# standardised total reproductive success #

INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
mod.start <- MCMCglmm(t_rs ~ ln_pg, random= ~ us(1 + ln_pg):common + animal, data = gsDataLong[-which(is.na(gsDataLong$ln_pg)),], nitt=1000, thin=1, burnin=0, pr=TRUE, pl=TRUE, prior=priorA, ginverse=list(animal=INtree$Ainv), family ="gaussian", verbose=FALSE)
mod <- mod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=mod$Liab[1,], R=mod$VCV[1,6], G=list(G1=matrix(ncol=2, nrow=2, mod$VCV[1,1:4]), G2=mod$VCV[1,5]))
  mod <- MCMCglmm(t_rs ~ ln_pg, random= ~ us(1 + ln_pg):common + animal, data = gsDataLong[-which(is.na(gsDataLong$ln_pg)),], ginverse=list(animal=INtree$Ainv), family ="gaussian", prior=priorA, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=999, start=start, verbose=FALSE)
  if(i > 300){
    mod.start$VCV[i-300,] <- mod$VCV[1,]
    mod.start$Sol[i-300,] <- mod$Sol[1,]
    mod.start$Liab[i-300,] <- mod$Liab[1,]
  }
}
rs_t_rs_modB.1 <- mod.start
rs_t_rs_modB.2 <- mod.start
rs_t_rs_modB.3 <- mod.start
save(rs_t_rs_modB.1, file=".../rs_t_rs_modB.1")
save(rs_t_rs_modB.2, file=".../rs_t_rs_modB.2")
save(rs_t_rs_modB.3, file=".../rs_t_rs_modB.3")

# chain convergence
hist(rs_t_rs_modB.1$Liab)
plot(rs_t_rs_modB.1$VCV)     		# all parameters well mixed
plot(rs_t_rs_modB.1$Sol[,1:2])	  	# intercept and slope well mixed
autocorr(rs_t_rs_modB.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(rs_t_rs_modB.1$Sol[,1:2]) 	# correlation between successive samples < 0.1 for all components
rs_t_rs_modB.Sols <- mcmc.list(list(rs_t_rs_modB.1$Sol[,1:2], rs_t_rs_modB.2$Sol[,1:2], rs_t_rs_modB.3$Sol[,1:2]))
plot(rs_t_rs_modB.Sols)
gelman.diag(rs_t_rs_modB.Sols)     		# upper CI = 1.00 suggesting convergence
heidel.diag(rs_t_rs_modB.1$VCV)    		# 
heidel.diag(rs_t_rs_modB.1$Sol[,1:2]) 	# intercept and slope passed halfwidth

# model parameters
summary(rs_t_rs_modB.1)
posterior.mode(rs_t_rs_modB.1$Sol[,1:2])		# intercept = 0.97, slope = -0.36
HPDinterval(rs_t_rs_modB.1$Sol[,1:2])		# intercept = 0.56 to 1.41, slope = -0.52 to -0.22

# intercept slope correlation
int.slope.cor <- rs_t_rs_modB.1$VCV[,2]/sqrt(rs_t_rs_modB.1$VCV[,1] * rs_t_rs_modB.1$VCV[,4])
posterior.mode(int.slope.cor)
autocorr(int.slope.cor)


## Two-Stage Models (metafor) ##

# standardised per capita reproductive success #

# get slope and sampling variance estimates
res.List <- lmList(groupedData(pc_rs ~ ln_pg | common, data = na.omit(gsDataLong[,c("pc_rs", "ln_pg", "common")])))
summary(res.List)		# 14/22 slopes positive
plot(augPred(res.List), grid=TRUE)
b <- lapply(res.List, coef)
V <- lapply(res.List, vcov)

# combine in a dataframe
metData <- data.frame(slope = unlist(lapply(b, "[[", 2)),  var = unlist(lapply(V, "[[", 4)))
metData$common <- rownames(metData)
metData$animal <- gsDataLong$animal[match(metData$common, gsDataLong$common)]

# run an intercept only model
ts_pc_rs_modB <- rma.mv(slope ~ 1, V=var, random = list(~ 1 | common, ~ 1 | animal), R = list(animal=birdCorB), data=metData)
summary(ts_pc_rs_modB)   	# B = 0.42, lwr = 0.20, upr = 0.64	, Q = 49, p = 0.0005


# standardised total reproductive success #

# get slope and sampling variance estimates
res.List <- lmList(groupedData(t_rs ~ ln_pg | common, data = na.omit(gsDataLong[,c("t_rs", "ln_pg", "common")])))
summary(res.List)		# 18/22 slopes negative
plot(augPred(res.List), grid=TRUE)
b <- lapply(res.List, coef)
V <- lapply(res.List, vcov)

# combine in a dataframe
metData <- data.frame(slope = unlist(lapply(b, "[[", 2)),  var = unlist(lapply(V, "[[", 4)))
metData$common <- rownames(metData)
metData$animal <- gsDataLong$animal[match(metData$common, gsDataLong$common)]

# run an intercept only model
ts_t_rs_modB <- rma.mv(slope ~ 1, V=var, random = list(~ 1 | common, ~ 1 | animal), R = list(animal=birdCorB), data=metData)
summary(ts_t_rs_modB)   	# B = -0.58, lwr = -0.73, upr = -0.44, Q = 27, p = 0.159


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

										### Part C: study effort ###

## include weighting term to account for different sample sizes between studies ##
## note however that the fact there are few large groups is biological and not due to sampling bias ##


# standardised per capita reproductive success ~ group size #

MEV <- 1/gsDataLong$N
INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
mod.start <- MCMCglmm(pc_rs ~ ln_gs, random= ~ us(1 + ln_gs):common + animal, data = gsDataLong, nitt=1000, thin=1, burnin=0, pr=TRUE, pl=TRUE, prior=priorA, ginverse=list(animal=INtree$Ainv), family ="gaussian", verbose=FALSE, mev=MEV)
mod <- mod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=mod$Liab[1,], R=mod$VCV[1,7], G=list(G1=matrix(ncol=2, nrow=2, mod$VCV[1,1:4]), G2=mod$VCV[1,5]))
  mod <- MCMCglmm(pc_rs ~ ln_gs, random= ~ us(1 + ln_gs):common + animal, data = gsDataLong, ginverse=list(animal=INtree$Ainv), family ="gaussian", prior=priorA, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=999, start=start, verbose=FALSE, mev=MEV)
  if(i > 300){
    mod.start$VCV[i-300,] <- mod$VCV[1,]
    mod.start$Sol[i-300,] <- mod$Sol[1,]
    mod.start$Liab[i-300,] <- mod$Liab[1,]
  }
}
rs_pc_rs_modCA.1 <- mod.start
rs_pc_rs_modCA.2 <- mod.start
rs_pc_rs_modCA.3 <- mod.start
save(rs_pc_rs_modCA.1, file=".../rs_pc_rs_modCA.1")
save(rs_pc_rs_modCA.2, file=".../rs_pc_rs_modCA.2")
save(rs_pc_rs_modCA.3, file=".../rs_pc_rs_modCA.3")

# chain convergence
hist(rs_pc_rs_modCA.1$Liab)
plot(rs_pc_rs_modCA.1$VCV)     		# all parameters well mixed
plot(rs_pc_rs_modCA.1$Sol[,1:2])	  	# intercept and slope estimates well mixed
autocorr(rs_pc_rs_modCA.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(rs_pc_rs_modCA.1$Sol[,1:2]) 	# correlation between successive samples < 0.1 for all components
rs_pc_rs_modCA.Sols <- mcmc.list(list(rs_pc_rs_modCA.1$Sol[,1:2], rs_pc_rs_modCA.2$Sol[,1:2], rs_pc_rs_modCA.3$Sol[,1:2]))
plot(rs_pc_rs_modCA.Sols)
gelman.diag(rs_pc_rs_modCA.Sols)     		# upper CI = 1.00 suggesting convergence
heidel.diag(rs_pc_rs_modCA.1$VCV)    		# all parameters passed halfwidth
heidel.diag(rs_pc_rs_modCA.1$Sol[,1:2]) 		# intercept and slope passed halfwidth

# model parameters
summary(rs_pc_rs_modCA.1)
posterior.mode(rs_pc_rs_modCA.1$Sol[,1:2])	 # int = 0.79, gs = -0.70
HPDinterval(rs_pc_rs_modCA.1$Sol[,1:2])		 # int = 0.24 to 1.80, gs = -1.30 to -0.07

# intercept slope correlation
int.slope.cor <- rs_pc_rs_modCA.1$VCV[,2]/sqrt(rs_pc_rs_modCA.1$VCV[,1] * rs_pc_rs_modCA.1$VCV[,4])
posterior.mode(int.slope.cor)
autocorr(int.slope.cor)


# standardised total reproductive success ~ group size #

MEV <- 1/gsDataLong$N
INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
mod.start <- MCMCglmm(t_rs ~ ln_gs, random= ~ us(1 + ln_gs):common + animal, data = gsDataLong, nitt=1000, thin=1, burnin=0, pr=TRUE, pl=TRUE, prior=priorA, ginverse=list(animal=INtree$Ainv), family ="gaussian", verbose=FALSE, mev=MEV)
mod <- mod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=mod$Liab[1,], R=mod$VCV[1,7], G=list(G1=matrix(ncol=2, nrow=2, mod$VCV[1,1:4]), G2=mod$VCV[1,5]))
  mod <- MCMCglmm(t_rs ~ ln_gs, random= ~ us(1 + ln_gs):common + animal, data = gsDataLong, ginverse=list(animal=INtree$Ainv), family ="gaussian", prior=priorA, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=999, start=start, verbose=FALSE, mev=MEV)
  if(i > 300){
    mod.start$VCV[i-300,] <- mod$VCV[1,]
    mod.start$Sol[i-300,] <- mod$Sol[1,]
    mod.start$Liab[i-300,] <- mod$Liab[1,]
  }
}
rs_t_rs_modCA.1 <- mod.start
rs_t_rs_modCA.2 <- mod.start
rs_t_rs_modCA.3 <- mod.start
save(rs_t_rs_modCA.1, file=".../rs_t_rs_modCA.1")
save(rs_t_rs_modCA.2, file=".../rs_t_rs_modCA.2")
save(rs_t_rs_modCA.3, file=".../rs_t_rs_modCA.3")

# chain convergence
hist(rs_t_rs_modCA.1$Liab)
plot(rs_t_rs_modCA.1$VCV)     		# all parameters well mixed
plot(rs_t_rs_modCA.1$Sol[,1:2])	  	# intercept and slope estimates well mixed
autocorr(rs_t_rs_modCA.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(rs_t_rs_modCA.1$Sol[,1:2]) 	# correlation between successive samples < 0.1 for all components
rs_t_rs_modCA.Sols <- mcmc.list(list(rs_t_rs_modCA.1$Sol[,1:2], rs_t_rs_modCA.2$Sol[,1:2], rs_t_rs_modCA.3$Sol[,1:2]))
plot(rs_t_rs_modCA.Sols)
gelman.diag(rs_t_rs_modCA.Sols)     		# upper CI = 1.00 suggesting convergence
heidel.diag(rs_t_rs_modCA.1$VCV)    		# all parameters passed halfwidth
heidel.diag(rs_t_rs_modCA.1$Sol[,1:2])  # intercept and slope passed halfwidth

# model parameters
summary(rs_t_rs_modCA.1)
posterior.mode(rs_t_rs_modCA.1$Sol[,1:2])	# int = -1.64, gs = 1.10
HPDinterval(rs_t_rs_modCA.1$Sol[,1:2])		# int = -2.11 to -1.12, gs = 0.78 to 1.57

# intercept slope correlation
int.slope.cor <- rs_t_rs_modCA.1$VCV[,2]/sqrt(rs_t_rs_modCA.1$VCV[,1] * rs_t_rs_modCA.1$VCV[,4])
posterior.mode(int.slope.cor)
autocorr(int.slope.cor)


# standardised per capita reproductive success ~ the frequency of groups #

MEV <- 1/gsDataLong$N[-which(is.na(gsDataLong$ln_pg))]
INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
mod.start <- MCMCglmm(pc_rs ~ ln_pg, random= ~ us(1 + ln_pg):common + animal, data = gsDataLong[-which(is.na(gsDataLong$ln_pg)),], nitt=1000, thin=1, burnin=0, pr=TRUE, pl=TRUE, prior=priorA, ginverse=list(animal=INtree$Ainv), family ="gaussian", verbose=FALSE, mev=MEV)
mod <- mod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=mod$Liab[1,], R=mod$VCV[1,7], G=list(G1=matrix(ncol=2, nrow=2, mod$VCV[1,1:4]), G2=mod$VCV[1,5]))
  mod <- MCMCglmm(pc_rs ~ ln_pg, random= ~ us(1 + ln_pg):common + animal, data = gsDataLong[-which(is.na(gsDataLong$ln_pg)),], ginverse=list(animal=INtree$Ainv), family ="gaussian", prior=priorA, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=999, start=start, verbose=FALSE, mev=MEV)
  if(i > 300){
    mod.start$VCV[i-300,] <- mod$VCV[1,]
    mod.start$Sol[i-300,] <- mod$Sol[1,]
    mod.start$Liab[i-300,] <- mod$Liab[1,]
  }
}
rs_pc_rs_modCB.1 <- mod.start
rs_pc_rs_modCB.2 <- mod.start
rs_pc_rs_modCB.3 <- mod.start
save(rs_pc_rs_modCB.1, file=".../rs_pc_rs_modCB.1")
save(rs_pc_rs_modCB.2, file=".../rs_pc_rs_modCB.2")
save(rs_pc_rs_modCB.3, file=".../rs_pc_rs_modCB.3")

# chain convergence
hist(rs_pc_rs_modCB.1$Liab)
plot(rs_pc_rs_modCB.1$VCV)     		# all parameters well mixed
plot(rs_pc_rs_modCB.1$Sol[,1:2])	  	# intercept and slope estimates well mixed
autocorr(rs_pc_rs_modCB.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(rs_pc_rs_modCB.1$Sol[,1:2]) 	# correlation between successive samples < 0.1 for all components
rs_pc_rs_modCB.Sols <- mcmc.list(list(rs_pc_rs_modCB.1$Sol[,1:2], rs_pc_rs_modCB.2$Sol[,1:2], rs_pc_rs_modCB.3$Sol[,1:2]))
plot(rs_pc_rs_modCB.Sols)
gelman.diag(rs_pc_rs_modCB.Sols)     		# upper CI = 1.00 suggesting convergence
heidel.diag(rs_pc_rs_modCB.1$VCV)    		# all parameters passed halfwidth
heidel.diag(rs_pc_rs_modCB.1$Sol[,1:2]) 		# intercept and slope passed halfwidth

# model parameters
summary(rs_pc_rs_modCB.1)
posterior.mode(rs_pc_rs_modCB.1$Sol[,1:2])		# int = -0.87, slope = 0.32
HPDinterval(rs_pc_rs_modCB.1$Sol[,1:2])			# int = -1.20 to -0.11, slope = 0.08 to 0.47

# intercept slope correlation
int.slope.cor <- rs_pc_rs_modCB.1$VCV[,2]/sqrt(rs_pc_rs_modCB.1$VCV[,1] * rs_pc_rs_modCB.1$VCV[,4])
posterior.mode(int.slope.cor)
autocorr(int.slope.cor)


# standardised total reproductive success ~ the frequency of groups #

MEV <- 1/gsDataLong$N[-which(is.na(gsDataLong$ln_pg))]
INtree <- inverseA(gsTrees[[1]], nodes="TIPS")
mod.start <- MCMCglmm(t_rs ~ ln_pg, random= ~ us(1 + ln_pg):common + animal, data = gsDataLong[-which(is.na(gsDataLong$ln_pg)),], nitt=1000, thin=1, burnin=0, pr=TRUE, pl=TRUE, prior=priorA, ginverse=list(animal=INtree$Ainv), family ="gaussian", verbose=FALSE, mev=MEV)
mod <- mod.start
for(i in 1:1300){
  INtree <- inverseA(gsTrees[[i]], nodes="TIPS")
  start <- list(Liab=mod$Liab[1,], R=mod$VCV[1,7], G=list(G1=matrix(ncol=2, nrow=2, mod$VCV[1,1:4]), G2=mod$VCV[1,5]))
  mod <- MCMCglmm(t_rs ~ ln_pg, random= ~ us(1 + ln_pg):common + animal, data = gsDataLong[-which(is.na(gsDataLong$ln_pg)),], ginverse=list(animal=INtree$Ainv), family ="gaussian", prior=priorA, pr=TRUE, pl=TRUE, nitt=1000, thin=1, burnin=999, start=start, verbose=FALSE, mev=MEV)
  if(i > 300){
    mod.start$VCV[i-300,] <- mod$VCV[1,]
    mod.start$Sol[i-300,] <- mod$Sol[1,]
    mod.start$Liab[i-300,] <- mod$Liab[1,]
  }
}
rs_t_rs_modCB.1 <- mod.start
rs_t_rs_modCB.2 <- mod.start
rs_t_rs_modCB.3 <- mod.start
save(rs_t_rs_modCB.1, file=".../rs_t_rs_modCB.1")
save(rs_t_rs_modCB.2, file=".../rs_t_rs_modCB.2")
save(rs_t_rs_modCB.3, file=".../rs_t_rs_modCB.3")

# chain convergence
hist(rs_t_rs_modCB.1$Liab)
plot(rs_t_rs_modCB.1$VCV)     		# all parameters well mixed
plot(rs_t_rs_modCB.1$Sol[,1:2])	  	# intercept and slope well mixed
autocorr(rs_t_rs_modCB.1$VCV)   		# correlation between successive samples < 0.1 for all components
autocorr(rs_t_rs_modCB.1$Sol[,1:2]) 	# correlation between successive samples < 0.1 for all components
rs_t_rs_modCB.Sols <- mcmc.list(list(rs_t_rs_modCB.1$Sol[,1:2], rs_t_rs_modCB.2$Sol[,1:2], rs_t_rs_modCB.3$Sol[,1:2]))
plot(rs_t_rs_modCB.Sols)
gelman.diag(rs_t_rs_modCB.Sols)     		# upper CI = 1.00 suggesting convergence
heidel.diag(rs_t_rs_modCB.1$VCV)    		# units passed halfwidth
heidel.diag(rs_t_rs_modCB.1$Sol[,1:2])	# intercept and slope passed halfwidth

# model parameters
summary(rs_t_rs_modCB.1)
posterior.mode(rs_t_rs_modCB.1$Sol[,1:2])	# int = 1.01, pg = -0.36
HPDinterval(rs_t_rs_modCB.1$Sol[,1:2])		# int = 0.55 to 1.38, pg = -0.54 to -0.23

# intercept slope correlation
int.slope.cor <- rs_t_rs_modCB.1$VCV[,2]/sqrt(rs_t_rs_modCB.1$VCV[,1] * rs_t_rs_modCB.1$VCV[,4])
posterior.mode(int.slope.cor)
autocorr(int.slope.cor)


# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #

											### THE END ###

# ------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------------------------------------------------------ #