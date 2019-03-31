#-----------------------------------------------------------------------------------------
# November 9 2015
# MK Lab Meeting from Steph Peacock <stephanie.peacock@ualberta.ca>
# Using JAGS for Bayesian analyses
# with some material stolen from 
# http://www.johnmyleswhite.com/notebook/2010/08/20/using-jags-in-r-with-the-rjags-package/
#-----------------------------------------------------------------------------------------

library(dclone) #automatically loads rjags
library(lattice)
library(lme4)

#setwd("~/Google Drive/Data cloning")

#######################################################################################
## A. Simple Ricker model
#######################################################################################

#-----------------------------------------------------------------------------------------
# 1) Simulate data

a<-1.2			# growth rate
b<-1.2/1200	# density dependence
sig<-0.8		# process error
N<-100			# timesteps

y<-numeric(N+1)	# vector of population sizes
y1<-600			# initial condition

# Simulate data:
set.seed(9823)
y[1]<-y1
for(i in 2:(N+1)) y[i]<-y[i-1]*exp(rnorm(n=1, a-b*y[i-1], sig))

S<-y[1:N];		# spawners
Y<-log(y[2:(N+1)]/S)	# productivity = log(recruits/spawner)

plot(y, type="l", xlab="time", ylab="spawners"); abline(h=a/b, lty=3)
plot(S, Y, xlab="spawners", ylab="productivity"); abline(a=a, b=-b, col=2); abline(h=0)

#-----------------------------------------------------------------------------------------
# 2) Fit using simple lm() in R:
fit1<-lm(Y~S)
summary(fit1)

#-----------------------------------------------------------------------------------------
# 3) Fit using optim()
Ricker.logLik<-function(params){
  sig<-exp(params[3])
  Y.hat<-params[1]-params[2]*S
  LL<-dnorm(Y, mean=Y.hat, sd=sig, log=TRUE)
  return(-sum(LL))
}

# What is the logLikelihood of our parameters from lm()
Ricker.logLik(c(summary(fit1)$coefficients[1,1], -summary(fit1)$coefficients[2,1], log(summary(fit1)$sigma)))

fit2<-optim(par=c(0,0,0), fn=Ricker.logLik)
fit2

#-----------------------------------------------------------------------------------------
# 4) Fit using JAGS

# Need to write model function IN JAGS LANGUAGE (see user manual)
ricker.model<-function(){	
  # define priors on parameter (mean and PRECISION)
  a ~ dnorm(0, 0.0001)
  b ~ dnorm(0, 0.0001)
  sigma ~ dunif(0, 100)
  tau <- pow(sigma, -2)	
  
  # simulate model and likelihood
  for(i in 1:N){ 					# for each data point
    Y.hat[i] <- a-b*S[i] 		# expected value based on parameters
    Y[i] ~ dnorm(Y.hat[i], tau) # probability of observation Y given expectation Y.hat - this is where you're FITTING based on predictions
  }
}
#for the dc.fit here, because it's a for loop, the list objects are your data, make sure the data are in the same order (i.e. daphnia 5 needs to be at index 5 in the different dfs)
fit3<-jags.fit(data=list('Y'=Y, 'S'=S, 'N'=N), params=c("a", "b", "sigma"), model=ricker.model, n.chains=3, n.adapt=100)

plot(fit3)
summary(fit3) 

#-----------------------------------------------------------------------------------------
# 5) Data cloning: Are parameters estimable?

#very little is different here, but tell it what to multiply by ('N' in this case)
fit4<-dc.fit(data=list('Y'=Y, 'S'=S, 'N'=N), params=c("a", "b", "sigma"), model=ricker.model, n.clones=c(1:5), multiply=c("N"), n.chains=3, n.adapt=100)

dct<-dctable(fit4)
plot(dct, type="var")
# Looks good!

#######################################################################################
## B. Hierarchical Ricker model
#######################################################################################

#-----------------------------------------------------------------------------------------
# 1) Simulate data

a<-1.2								# growth rate
b<-1.2/1200							# density dependence
sig<-0.3							# process error
sigP<-0.8							# sd in growth rate among populations

Nt<-50								# number of timesteps
Np<-20								# number of populations
thetaP<-rnorm(Np, mean=0, sd=sigP)	# population random-effect on growth rate

t<-rep(1:(Nt+1), Np)				# vector of timesteps
P<-rep(1:Np, each=(Nt+1))			# vector of population numbers
y<-numeric((Nt+1)*Np)				# vector of population sizes (to be filled)
y[which(t==1)]<-runif(Np, 100, 1000)# initialize with random number between 100 and 1000

# Simulate data:
for(i in 1:Np){ #for each population
  for(j in 2:(Nt+1)){ #for each timestep
    y[which(P==i&t==j)]<-y[which(P==i&t==(j-1))]*exp(rnorm(n=1, (a+thetaP[i])-b*y[which(P==i&t==(j-1))], sig))
  }
}

# Make dataframe of spawners and productivity
data<-data.frame(S=numeric(Nt*Np), Y=numeric(Nt*Np), Pop=numeric(Nt*Np))
N<-Nt*Np

k<-0
for(i in 1:Np){
  for(j in 1:Nt){
    k<-k+1
    data$S[k]<- y[which(P==i&t==j)]
    data$Y[k]<- log(y[which(P==i&t==(j+1))]/y[which(P==i&t==j)])
    data$Pop[k]<-i
  }
}

xyplot(Y~S|Pop, data=data)

#-----------------------------------------------------------------------------------------
# 2) Fit using simple lmer() in R:

hfit1<-lmer(Y~S+(1|Pop), data=data)
summary(hfit1)

hfit1.par<-c(summary(hfit1)$coefficients[1,1], -summary(hfit1)$coefficients[2,1], summary(hfit1)$sigma, attr(summary(hfit1)$varcor[[1]], "stddev"))
#-----------------------------------------------------------------------------------------
# 4) Fit using JAGS

# Transform data to matrices to make model syntax easier
S<-matrix(data$S, nrow=Np, ncol=Nt, byrow=TRUE)
Y<-matrix(data$Y, nrow=Np, ncol=Nt, byrow=TRUE)

H.ricker.model<-function(){
  # Priors on parameters
  a ~ dnorm(0, 0.0001)
  b ~ dnorm(0, 0.0001)
  sig ~ dunif(0, 100)
  sigP ~ dunif(0, 100)
  tauP <- pow(sigP, -2)
  tau <- pow(sig, -2)
  
  # Hierarchical part - a different thetaP for each population:
  for(i in 1:Np){
    thetaP[i] ~ dnorm(0, tauP)
  }
  
  # Model simulation and likelihood
  for(j in 1:Np){ #for each population
    for(i in 1:Nt){ #for each time step
      Y.hat[j,i] <- a + thetaP[j] - b*S[j,i] # calculate model prediction
      Y[j,i] ~ dnorm(Y.hat[j,i], tau) #likelihood
    } #end i
  } # end j	
}


hfit2<-jags.fit(data=list('S'=S, 'Y'=Y, 'Np'=Np, 'Nt'=Nt), params=c("a", "b", "sig", "sigP"), model=H.ricker.model)

S<-summary(hfit2)
par(mfrow=c(2,2))
for(i in 1:4){
  hist(hfit2[[1]][,i], main=colnames(hfit2[[1]])[i], col="#00000030", border=NA, xlab="", freq=FALSE)
  hist(hfit2[[2]][,i], main=colnames(hfit2[[2]])[i], col="#FF000030", add=TRUE, freq=FALSE, border=NA)
  hist(hfit2[[3]][,i], main=colnames(hfit2[[3]])[i], col="#00FF0030", add=TRUE, freq=FALSE, border=NA)
  abline(v=S[[1]][i,1], lty=2)
  abline(v=c(a,b,sig,sigP)[i])
  abline(v=hfit1.par[i], lty=3)
}
legend("topright", lty=c(1,2,3), c("Real", "JAGS", "lmer"))

#######################################################################################
## C. Hierarchical Ricker model with observation and process error
#######################################################################################

#-----------------------------------------------------------------------------------------
# 1) Simulate data

a<-1.2			# growth rate
b<-1.2/1200		# density dependence
sig.proc<-0.2		# process error
sig.obs<-0.15		# observation error
N<-100			# timesteps

y<-numeric(N)	# vector of population sizes
y1<-600			# initial condition

# Simulate data:
set.seed(9823)
# Process (same as in A)
y[1]<-y1
for(i in 2:N) y[i]<-y[i-1]*exp(rnorm(n=1, a-b*y[i-1], sig.proc))
# Observation
Y<-y*exp(rnorm(N, 0, sig.obs))

plot(y, type="l", xlab="time", ylab="spawners"); abline(h=a/b, lty=3)
lines(Y, col=2)

#-----------------------------------------------------------------------------------------
# 2) Fit using JAGS

# Need to write model function IN JAGS LANGUAGE (see user manual)
SS.ricker.model<-function(){		
  # define priors on parameter (normal distn: mean and PRECISION)
  a ~ dnorm(0, 0.0001)
  b ~ dnorm(0, 0.0001)
  sig.proc ~ dunif(0, 100)
  sig.obs ~ dunif(0, 100)
  tau.proc <- pow(sig.proc, -2)	
  tau.obs <- pow(sig.obs, -2)	
  
  # simulate timeseries of TRUE values
  y.true[1]<-y1 # initialize true data
  for(i in 2:N){ 					
    y.true[i] ~ dlnorm(log(y.true[i-1]*exp(a-b*y.true[i-1])), tau.proc)
  }
  
  # interface with data: observation error
  for(i in 1:N){
    Y[i] ~ dlnorm(log(y.true[i]), tau.obs)
  }
  
}

SSfit<-jags.fit(data=list('Y'=Y, 'N'=N, 'y1'=y1), params=c("a", "b", "sig.proc", "sig.obs"), model=SS.ricker.model, n.adapt=10000, n.iter=20000, thin=2)



plot(SSfit)
summary(SSfit)
gelman.diag(SSfit)

