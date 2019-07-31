###############################################################################
##  STORM PETREL ABUNDANCE ESTIMATION ON FILFLA ISLAND  ##
###############################################################################
## written by steffen.oppel@rspb.org.uk on 4 January 2014

library(RODBC)
library(R2jags)


################################################################################################################
###################### LOAD AND MANIPULATE DATA 	  ######################################################
################################################################################################################

####################### LOADING RINGING DATA           #########################################################
setwd("A:\\RSPB\\Malta\\Raw_data")
setwd("C:\\STEFFEN\\RSPB\\Malta\\Raw_data")

SP<-odbcConnectAccess2007('ESP_Malta.accdb')
ESP <- sqlQuery(SP, "SELECT * FROM capthist")  ## changed query to only have data from 2000
odbcClose(SP)




### FORMAT DATA ###
ESP[is.na(ESP)]<-0
yobs<-as.matrix(ESP[,2:dim(ESP)[2]])
yobs<-ifelse(yobs>0,1,0)



#######################################################################
#
# 6. Estimation of the size of a closed population
# 
#######################################################################


# Augment data set by 5000 potential individuals
nz <- 20000
yaug <- rbind(yobs, array(0, dim = c(nz, dim(yobs)[2])))


# Specify model in BUGS language
sink("model.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)			# Inclusion indicators
   for (j in 1:T){
      yaug[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p		# Can only be detected if z=1
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(yaug = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), p = runif(1, 0, 1))

# Parameters monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 25000
nt <- 5
nb <- 15000
nc <- 3

# Call WinBUGS from R (BRT 2.5 hrs)
#out <- bugs(win.data, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())
out <- jags(win.data, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
out2<-out$BUGSoutput


# Summarize posteriors
print(out2, dig = 3)
hist(out2$sims.list$N, nclass = 150, col = "gray", main = "European Storm Petrel population Filfla Island", ylab = "likelihood to be correct size", xlab = "Population size", las = 1, xlim = c(0, 30000))
abline(v = dim(yobs)[1], lwd = 3, col='red')



# 6.2.4. Individual (random) effects: the heterogeneity model Mh
# Define function to simulate data under Mh
data.fn <- function(N = 100, mean.p = 0.4, T = 5, sd = 1){
   yfull <- yobs <- array(NA, dim = c(N, T) )
   mean.lp <- log(mean.p / (1-mean.p))
   p.vec <- plogis(mean.lp+ rnorm(N, 0, sd))

   for (i in 1:N){
      yfull[i,] <- rbinom(n = T, size = 1, prob = p.vec[i])
      }

   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   hist(p.vec, xlim = c(0,1), nclass = 20, col = "gray", main = "", xlab = "Detection probability", las = 1)
   return(list(N = N, p.vec = p.vec, mean.lp = mean.lp, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data <- data.fn()

# Aggregate capture-histories and augment data set
y <- sort(apply(data$yobs, 1, sum), decreasing = TRUE)
nz <- 300
yaug <- c(y, rep(0, nz))
yaug

# Specify model in BUGS language
sink("model.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
mean.lp <- logit(mean.p)
mean.p ~ dunif(0, 1)
tau <- 1 / (sd * sd)
sd ~ dunif(0, 5)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   logit(p[i]) <- eps[i]
   eps[i] ~ dnorm(mean.lp, tau)I(-16, 16)	# See web appendix A in Royle (2009)
   p.eff[i] <- z[i] * p[i]
   y[i] ~ dbin(p.eff[i], T)
   }

# Derived quantities
N <- sum(z[])
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = yaug, M = length(yaug), T = ncol(data$yobs))

# Initial values
inits <- function() list(z = rep(1, length(yaug)), sd = runif(1, 0.1, 0.9))

# Parameters monitored
params <- c("N", "mean.p", "sd", "omega")

# MCMC settings
ni <- 25000
nt <- 2
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 6 min)
out <- bugs(win.data, inits, params, "model.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)
hist(out$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 250))
abline(v = data$C, col = "black", lwd = 3)


# 6.2.5. Combined effects: model Mth
# Define function to simulate data under Mth
data.fn <- function(N = 100, T = 5, mean.p = 0.4, time.effects = runif(T, -1, 1), sd = 1){
   yfull <- yobs <- p <- array(NA, dim = c(N, T) )
   mean.lp <- log(mean.p / (1-mean.p))         # mean p on logit scale
   eps <- rnorm(N, 0, sd)                      # Individual effects

   for (j in 1:T){
      pp <- p[,j] <- plogis(mean.lp + time.effects[j] + eps)
      yfull[,j] <- rbinom(n = N, size = 1, prob = pp)
      }

   ever.detected <- apply(yfull, 1, max)
   C <- sum(ever.detected)
   yobs <- yfull[ever.detected == 1,]
   cat(C, "out of", N, "animals present were detected.\n")
   cat("Mean p per occasion:", round(apply(p, 2, mean), 2), "\n")
   par(mfrow = c(2,1))
   plot(plogis(mean.lp + time.effects), xlab = "Occasion", type = "b", main = "Approx. mean p at each occasion", ylim = c(0, 1))
   hist(plogis(mean.lp + eps), xlim = c(0, 1), col = "gray", main = "Approx. distribution of p at average occasion")
   return(list(N = N, mean.lp = mean.lp, time.effects = time.effects, sd = sd, eps = eps, C = C, T = T, yfull = yfull, yobs = yobs))
   }

data <- data.fn()

# data<-data.fn(T = 10, mean.p = 0.2, time.effects = runif(10, 0, 0), sd = 0)	# M0
# data<-data.fn(T = 10, mean.p = 0.5, time.effects = runif(10, 0, 0), sd = 1)	# Mh
# data <- data.fn(T = 10, sd = 0)	# Mt

# Augment data set
nz <- 300
yaug <- rbind(data$yobs, array(0, dim=c(nz, data$T)))

# Specify model in BUGS language
sink("model.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for (j in 1:T){
   mean.lp[j] <- log(mean.p[j] / (1 - mean.p[j]) ) # Define logit 
   mean.p[j] ~ dunif(0, 1)
   }
tau <- 1 / (sd * sd)
sd ~ dunif(0,5)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   eps[i] ~ dnorm(0, tau)I(-16, 16)	# See web appendix A in Royle (2009)
   for (j in 1:T){
      lp[i,j] <- mean.lp[j] + eps[i]
      p[i,j] <- 1 / (1 + exp(-lp[i,j])) 			# Define logit
      p.eff[i,j] <- z[i] * p[i,j]
      y[i,j] ~ dbern(p.eff[i,j])
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} 
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = yaug, M = nrow(yaug), T = ncol(yaug))

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), sd = runif(1, 0.1, 0.9))

# Parameters monitored
params <- c("N", "mean.p", "mean.lp", "sd", "omega")

# MCMC settings
ni <- 25000
nt <- 2
nb <- 5000
nc <- 3

# Call WinBUGS from R (BRT 47 min)
out <- bugs(win.data, inits, params, "model.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(out, dig = 3)
hist(out$sims.list$N, nclass = 50, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(80, 200))
abline(v = data$C, col = "black", lwd = 3)


# 6.3. Analysis of a real data set: model Mtbh for species richness estimation
# Read in data and look at them
p610 <- read.table("p610.txt", header = TRUE)
y <- p610[,5:9]                           # Grab counts
y[y > 1] <- 1                             # Counts to det-nondetections
C <- sum(apply(y, 1, max)) ; print(C)     # Number of observed species
table(apply(y, 1, sum))                   # Capture-frequencies

# Specify model in BUGS language
sink("M_tbh.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for  (j in 1:T){
   alpha[j] <- log(mean.p[j] / (1-mean.p[j])) # Define logit 
   mean.p[j] ~ dunif(0, 1) 	# Detection intercepts
   }
gamma ~ dnorm(0, 0.01)
tau <- 1 / (sd * sd)
sd ~ dunif(0, 3)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   eps[i] ~ dnorm(0, tau)I(-16, 16)

   # First occasion: no term for recapture (gamma)
   y[i,1] ~ dbern(p.eff[i,1])
   p.eff[i,1] <- z[i] * p[i,1]
   p[i,1] <- 1 / (1 + exp(-lp[i,1]))
   lp[i,1] <- alpha[1] + eps[i]

   # All subsequent occasions: includes recapture term (gamma)
   for (j in 2:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))   
      lp[i,j] <- alpha[j] + eps[i] + gamma * y[i,(j-1)]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} 
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = as.matrix(y), M = nrow(y), T = ncol(y))

# Initial values
inits <- function() list(z = rep(1, nrow(y)), sd = runif(1, 0.1, 0.9))

# Parameters monitored
params <- c("N", "mean.p", "gamma", "sd", "omega")

# MCMC settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 24 min)
out <- bugs(win.data, inits, params, "M_tbh.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors and plot posterior for N
print(out, dig = 3)
par(mfrow = c(1,2))
hist(out$sims.list$N, breaks = 35, col = "gray", main = "", xlab = "Community size", las = 1, xlim = c(30, 100), freq = FALSE)
abline(v = C, col = "black", lwd = 3)

# Define model
sink("M0.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
p ~ dunif(0, 1)

# Likelihood
for (i in 1:M){
   z[i] ~ dbern(omega)
   for (j in 1:T){
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p
      } #j
   } #i

# Derived quantities
N <- sum(z[])
} # end model
",fill = TRUE)
sink()

# Initial values
inits <- function() list(z = rep(1, nrow(y)))

# Define parameters to be monitored
params <- c("N", "p", "omega")

# MCMC settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 1 min)
out0 <- bugs(win.data, inits, params, "M0.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = bugs.dir, working.directory = getwd())

# Inspect output
print(out0, dig = 3)


# 6.4. Capture-recapture models with individual covariates: model Mt+X
# 6.4.1. Individual covariate model for species richness estimation
p610 <- read.table("p610.txt", header = TRUE)
y <- p610[,5:9]                         # Grab counts
y[y > 1] <- 1                           # Convert to det-nondetections
ever.observed <- apply(y, 1, max)
wt <- p610$bm[ever.observed == 1]       # Body mass
yy <- as.matrix(y[ever.observed == 1,]) # Detection histories
dimnames(yy) <- NULL

mlog <- mean(log(p610$bm^(1/3)))
sdlog <- sd(log(p610$bm^(1/3)))
hist(p610$bm^(1/3), xlim = c(0, 30), nclass = 25, freq = FALSE, col = "gray")
lines(density(rlnorm(n = 10^6, meanlog = mlog, sdlog = sdlog)), col = "blue", lwd = 3)

# Augment both data sets
nz = 150
yaug <- rbind(yy, array(0, dim = c(nz, ncol(yy))))
logwt3 <- c(log(wt^(1/3)), rep(NA, nz))

# Specify model in BUGS language
sink("M_t+X.txt")
cat("
model {

# Priors
omega ~ dunif(0, 1)
for (j in 1:T){
   alpha[j] <- log(mean.p[j] / (1-mean.p[j]))
   mean.p[j] ~ dunif(0, 1)
   }
beta ~ dnorm(0, 0.01)
mu.size ~ dnorm(0, 0.01)
tau.size <- 1 / pow(sd.size, 2)
sd.size ~ dunif(0, prior.sd.upper)   # Provide upper bound as data

# Likelihood
for (i in 1:M){  # Loop over individuals
   z[i] ~ dbern(omega)
   size[i] ~ dnorm(mu.size, tau.size)I(-6, 6)
   for (j in 1:T){  # Loop over occasions
      y[i,j] ~ dbern(p.eff[i,j])
      p.eff[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp[i,j]))   
      lp[i,j] <- alpha[j] + beta * size[i]
      } #j
   } #i

# Derived quantities
N <- sum(z[])
}
",fill = TRUE)
sink()

# Bundle data
win.data <- list(y = yaug, size = logwt3 - mean(logwt3, na.rm = TRUE), M = nrow(yaug), T = ncol(yaug), prior.sd.upper = 3)

# Initial values
inits <- function() list(z = rep(1, nrow(yaug)), beta = runif(1, 0, 1), mu.size = rnorm(1, 0, 1))

# Parameters monitored
params <- c("N", "mean.p", "beta", "omega", "mu.size", "sd.size")

# MCMC settings
ni <- 50000
nt <- 4
nb <- 10000
nc <- 3

# Call WinBUGS from R (BRT 19 min)
outX <- bugs(win.data, inits, params, "M_t+X.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors and plot posterior for N
print(outX, dig = 3)

hist(outX$sims.list$N, breaks = 100, col = "gray", main = "", xlab = "Community size", las = 1, xlim = c(30, 100), freq = FALSE)
abline(v = 31, col = "black", lwd = 3)

pred.wt <- seq(5, 2000, length.out = 100)	# Cov. vals for prediction
pred.wt.st <- log(pred.wt^(1/3))- mean(logwt3, na.rm = TRUE) # Transform them in the same was as in the analysis
pred.p<- plogis(log(mean(outX$mean$mean.p)/(1- mean(outX$mean$mean.p))) + outX$mean$beta * pred.wt.st) # Compute predicted response
plot(pred.wt, pred.p, type = "l", lwd = 3, col = "blue", las = 1, frame.plot = FALSE, ylim = c(0, 0.5))


# 6.4.2. Individual covariate model for population size estimation
# Read in data and look at shell width distribution
pinna <- read.table("pinna.txt", header = TRUE)
y <- cbind(pinna$d1, pinna$d2)
size <- pinna$width
hist(size, col = "gray", nclass = 50, xlim = c(0, 30), freq = FALSE)
lines(density(rnorm(10^6, mean = mean(size), sd = sd(size))), col = "blue", lwd = 3)

# Augment both data sets
nz = 150
yaug <- rbind(y, array(0, dim = c(nz, ncol(y))))
size <- c(size, rep(NA, nz))

# Bundle data
win.data <- list(y = yaug, size = size - mean(size, na.rm = TRUE), M = nrow(yaug), T = ncol(yaug), prior.sd.upper = 5)

# MCMC settings
ni <- 2500
nt <- 2
nb <- 500
nc <- 3

# Call WinBUGS from R (BRT 1 min)
outXX <- bugs(win.data, inits, params, "M_t+X.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(outXX, dig = 2)

Plot posterior for N and prediction of p
par(mfrow = c(1,2), mar = c(4.5, 4, 2, 1))
hist(outXX$sims.list$N, breaks = 30, col = "gray", main = "", xlab = "Population size", las = 1, xlim = c(143, 220), freq = FALSE)
abline(v = 143, col = "black", lwd = 3)

pred.size <- seq(0, 30, length.out = 1000)	# Cov. vals for prediction
pred.size.st <- pred.size - mean(size, na.rm = TRUE) # Transform them
pred.p<- plogis(log(mean(outXX$mean$mean.p)/(1- mean(outXX$mean$mean.p))) + outXX$mean$beta * pred.size.st) # Compute predicted detection prob.
plot(pred.size, pred.p, type = "l", lwd = 3, col = "blue", las = 1, frame.plot = FALSE, ylim = c(0, 1), xlab = "Shell width (cm)", ylab = "Predicted detection probability")










