library(tibble)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(mvtnorm)
library(coda)
library(rjags)
library("Rlab")   
set.seed(16)

                                                           ### ESERCIZIO 1 ###

# def functions
media <- function(parameter, posterior) {return(sum(parameter * posterior)*3/10**4)}
var <- function(parameter, posterior, m) {return(sum((parameter-m)**2*posterior)*3/10**4)}
percentile <- function(posterior, parameter, perc) {
  cumulative <- cumsum(posterior*3/10**4)
  index <- which.min(abs(cumulative - perc))
  return(parameter[index])
}

# data
death <- c(0:5)
corp1 <- c(109, 65, 22, 3, 1, 0)
corp2 <- c(144, 91, 32, 11, 2, 0)

lambda <- seq(0.0001, 3, length=10**4)


# UNIFORM PRIOR (with 'brutal force')

#Corp1

# def useful functions:
# - single likelihood (with only one data)
single_log_like <- function(n, lambda) {
  P <- dpois(n, lambda)
  return(log(P))
}
# - total likelihood (with all data)
total_log_like <- function(n_death, count, fun, lambda) {
  # n_death = number of die soldiers in one interval
  # count = number of intervals
  log_like <- 0
  for (i in c(1:length(count)) ) {
    log_like <- log_like + count[i]*fun(n_death[i], lambda)
  }
  return(log_like)
}

# def the prior
unif <- rep(1, 10**4)

# compute the logLikelihood
Log_like1 <- rep(0,10**4)
for (i in c(1:10**4)) {Log_like1[i] <- total_log_like(death, corp1, single_log_like, lambda[i])}

# in general log(like*prior) = log(like)+log(prior) but in this case log(prior)=log(1)=0
Log_post1 <- Log_like1                                   

#normalization
Z <- sum(exp(Log_post1)*3/10**4)
Post_corp1 <- exp(Log_post1)/Z
#plot(lambda, Post_corp1, type='l')


#Corp2

# compute the logLikelihood
Log_like2 <- rep(0,10**4)
for (i in c(1:10**4)) {Log_like2[i] <- total_log_like(death, corp2, single_log_like, lambda[i])}

# in general log(like*prior) = log(like)+log(prior) but in this case log(prior)=log(1)=0
Log_post2 <- Log_like2                               

#normalization
Z <- sum(exp(Log_post2)*3/10**4)
Post_corp2 <- exp(Log_post2)/Z
#plot(lambda, Post_corp2, type='l')



# JEFFREYS' PRIOR (using the fact that the Gamma is the conjugate prior of the Poisson distribution,
# and that the Jeffreys' prior can be written as a gamma function

#Corp1

# Parameters for my Gamma distribution
alpha1 <- sum(death*corp1) + 0.5
lambda1 <- sum(corp1)
Log_post1_J <- dgamma(lambda, shape=alpha1, scale = 1/lambda1, log = FALSE)

# normalization
Z <- sum(Log_post1_J*3/10**4)
Post_corp1_J <- Log_post1_J/Z
#plot(lambda, Post_corp1_J, type='l')


#Corp2

# Parameters for my Gamma distribution
alpha2 <- sum(death*corp2) + 0.5
lambda2 <- sum(corp2)
Log_post2_J <- dgamma(lambda, shape=alpha2, scale = 1/lambda2, log = FALSE)

# normalization
Z <- sum(Log_post2_J*3/10**4)
Post_corp2_J <- Log_post2_J/Z
#plot(lambda, Post_corp2_J, type='l')



mean <- c()
median <- c()
variance <- c()
lower_limit <- c()
upper_limit <- c()

# all the posteriors inside a matrix
P <- matrix(NA, nrow=10**4, ncol=4)
P[,1] <- Post_corp1
P[,2] <- Post_corp2
P[,3] <- Post_corp1_J
P[,4] <- Post_corp2_J

for (i in c(1:4)) {
  mu <- media(lambda, P[,i])
  mean <- c(mean, mu)
  median <- c(median, percentile(P[,i], lambda, 0.5))
  variance <- c(variance, var(lambda, P[,i], mu))
  lower_limit <- c(lower_limit, percentile(P[,i], lambda, 0.025))
  upper_limit <- c(upper_limit, percentile(P[,i], lambda, 0.975))
}

title <- c('Corp 1 with uniform prior', 'Corp 2 with uniform prior', 'Corp 1 with Jeffreys prior', 'Corp 2 with Jeffreys prior')

pdf('Es1_5.pdf')
par(mfrow = c(2,2))
for (i in c(1:4)) {
  plot(lambda, P[,i], type='l', main=title[i], xlab=' lambda = death rate', ylab='p(lambda | data)', xlim=c(0,1.5))
  polygon (c(lower_limit[i], lambda[lambda>lower_limit[i] & lambda<upper_limit[i]], upper_limit[i]),
           c(0, P[,i][lambda>lower_limit[i] & lambda<upper_limit[i]], 0), col=rgb(1, 0, 0,0.5))
  abline(v = mean[i], col='green', lty=1)
  abline(v = median[i], col='blue', lty=2)
  legend(x = "topright", legend = c("Posterior for alpha", '95% cred. int.', "Mean", "Median"),
         lty = c(1, 1, 2, 2), col = c(1, 'red','green', 'blue'), lwd = 2, cex=0.6) 
  minor.tick(nx=5, ny=4, tick.ratio=0.5) 
}
out <- dev.off() 

summary1 <- tibble(name_prior=c('Corp1_Unif', 'Corp2_Unif', 'Corp1_Jeff', 'Corp2_Jeff'),
                   mean=mean, median=median, variance=variance, lower_limit=lower_limit, upper_limit=upper_limit) 


cat(' #########  ESERCIZIO 1   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('See Es1_5.pdf to visualize the result ', '\n', '\n')
cat('The following table summarises the obtained results:', '\n')
print(summary1)
cat('----------------------------------------------------------', '\n')	



                                                           ### ESERCIZIO 2 ###

# def metrop function
metrop <- function(func, thetaInit, Nburnin, Nsamp, sampleCov, verbose, 
                   demo=FALSE, ...) {
  
  Ntheta   <- length(thetaInit)
  thetaCur <- thetaInit
  funcCur  <- func(thetaInit, ...) # log10
  funcSamp <- matrix(data=NA, nrow=Nsamp, ncol=2+Ntheta) 
  # funcSamp will be filled and returned
  nAccept  <- 0
  acceptRate <- 0
  if(demo) {
    thetaPropAll <- matrix(data=NA, nrow=Nsamp, ncol=Ntheta)
  }
  
  for(n in 1:(Nburnin+Nsamp)) {
    
    # Metropolis algorithm. No Hastings factor for symmetric proposal
    if(is.null(dim(sampleCov))) { # theta and sampleCov are scalars
      thetaProp <- rnorm(n=1, mean=thetaCur, sd=sqrt(sampleCov))
    } else {
      thetaProp <- rmvnorm(n=1, mean=thetaCur, sigma=sampleCov, 
                           method="eigen")
    }
    funcProp  <- func(thetaProp, ...) 
    logMR <- sum(funcProp) - sum(funcCur) # log10 of the Metropolis ratio
    #cat(n, thetaCur, funcCur, ":", thetaProp, funcProp, "\n")
    if(logMR>=0 || logMR>log10(runif(1, min=0, max=1))) {
      thetaCur   <- thetaProp
      funcCur    <- funcProp
      nAccept    <- nAccept + 1
      acceptRate <- nAccept/n
    }
    if(n>Nburnin) {
      funcSamp[n-Nburnin,1:2] <- funcCur
      funcSamp[n-Nburnin,3:(2+Ntheta)] <- thetaCur
      if(demo) {
        thetaPropAll[n-Nburnin,1:Ntheta] <- thetaProp
      }
    }
    
    # Diagnostics
    if( is.finite(verbose) && (n%%verbose==0 || n==Nburnin+Nsamp) ) {
      s1 <- noquote(formatC(n,          format="d", digits=5, flag=""))
      s2 <- noquote(formatC(Nburnin,    format="g", digits=5, flag=""))
      s3 <- noquote(formatC(Nsamp,      format="g", digits=5, flag=""))
      s4 <- noquote(formatC(acceptRate, format="f", digits=4, width=7, 
                            flag=""))
      # cat(s1, "of", s2, "+", s3, s4, "\n")
    }
    
  }
  
  if(demo) {
    return(list(funcSamp=funcSamp, thetaPropAll=thetaPropAll))
  } else {
    return(funcSamp)
  }
}

# data
death <- c(0:5)
corp1 <- c(109, 65, 22, 3, 1, 0)
corp2 <- c(144, 91, 32, 11, 2, 0)

obsdata1 <- c(sum(death*corp1), sum(corp1))
obsdata2 <- c(sum(death*corp2), sum(corp2))

# def log10 unif prior function
logprior_unif <- function(lambda) {
  logPrior <- 0
  return(logPrior)
}

# def log10 Jeff prior function
logprior_Jeff <- function(lambda) {
  logPrior <- log10(dgamma(lambda, shape = 0.5, scale = 10**4))
  return(logPrior)
}

# def log10 likelihood function
loglike <- function(lambda , obsdata) {
  logLike <- log10(dgamma(lambda, shape=obsdata[1]+1, rate=obsdata[2]))
  return(logLike)
}

# def log10 posterior as input for metrop function with unif prior
logpost_unif <- function(lambda, obsdata) {
  logprior <- logprior_unif(lambda)
  if(is.finite(logprior )) { 
    return( c(logprior, loglike(lambda, obsdata)) )
  } else {
    return( c(-Inf , -Inf) )
  }
}

# def log10 posterior as input for metrop function with Jeff prior
logpost_Jeff <- function(lambda, obsdata) {
  logprior <- logprior_Jeff(lambda)
  if(is.finite(logprior )) { 
    return( c(logprior, loglike(lambda, obsdata)) )
  } else {
    return( c(-Inf , -Inf) )
  }
}

# compute the fours mcmc
lambdaInit <- 0.6
sampleCov <- 0.05
allSamp_1_unif <- metrop(func=logpost_unif , thetaInit=lambdaInit ,
                         Nburnin=200, Nsamp=10**5,
                         sampleCov=sampleCov , verbose=1e3, obsdata=obsdata1)

allSamp_1_Jeff <- metrop(func=logpost_Jeff , thetaInit=lambdaInit ,
                         Nburnin=200, Nsamp=10**5,
                         sampleCov=sampleCov , verbose=1e3, obsdata=obsdata1)

allSamp_2_unif <- metrop(func=logpost_unif , thetaInit=lambdaInit ,
                         Nburnin=200, Nsamp=10**5,
                         sampleCov=sampleCov , verbose=1e3, obsdata=obsdata2)

allSamp_2_Jeff <- metrop(func=logpost_Jeff , thetaInit=lambdaInit ,
                         Nburnin=200, Nsamp=10**5,
                         sampleCov=sampleCov , verbose=1e3, obsdata=obsdata2)

# all results in a single matrix
Samples <- matrix(c(allSamp_1_unif[,3], allSamp_1_Jeff[,3], allSamp_2_unif[,3], allSamp_2_Jeff[,3]), ncol=4)

mean <- c()
median <- c()
var <- c()
l_limit <- c()
h_limit <- c()

title <- c('Corp1 with uniform prior', 'Corp1 with Jeffrey prior',
           'Corp2 with uniform prior', 'Corp2 with Jeffrey prior')

pdf('Es2_MCMC_corp1.pdf')
par(mfrow = c(2,2))

# the following two for cycles compute the same operations, at the beginning was just one cycle but I have to replicate 
# it otherwise the plots, all in the same page, were not readable.

# creating plot and computing results for mcmc about corp1
for (i in c(1:2)) {
  jump <- seq(from=1, to=nrow(Samples), by=25) # thin by factor 25
  filter_allSamp <- Samples[jump,i]
  
  nr <- length(filter_allSamp)
  
  l <- quantile(filter_allSamp, probs = 0.025)
  h <- quantile(filter_allSamp, probs = 0.975)
  media <- mean(filter_allSamp)
  med <- quantile(filter_allSamp, probs = 0.5)
  v <- var(filter_allSamp)
  
  l_limit <- c(l_limit, l)
  h_limit <- c(h_limit, h)
  median <- c(median, med)
  mean <- c(mean, media)
  var <- c(var, v)
  
  plot(c(1:nr), filter_allSamp,
       type="l",
       xlab="iteration",
       ylab='lambda', main=title[i])
  abline(h=media, lwd=1, lty=2, col='red')
  legend(x = "topright", legend = c("Mean"),
         lty = c(2), col = c('red'), lwd = 1, box.lty=0) 
  minor.tick(nx=4, ny=4, tick.ratio=0.5) 
  
  postDen <- density(filter_allSamp, n=10**3)
  
  plot(postDen$x, postDen$y, type='l', 
       col='navy', lwd = 2,
       xlab='lambda',
       ylab="p(lambda | data)", main=title[i] )
  polygon (c(l, postDen$x[postDen$x>l & postDen$x<h], h),
           c(0, postDen$y[postDen$x>l & postDen$x<h], 0), col=rgb(1, 0, 0,0.5))
  abline(v=media, lwd=1, lty=2)
  legend(x = "topright", legend = c("Mean"),
         lty = c(2), col = c('black'), lwd = 1, box.lty=0)
  minor.tick(nx=4, ny=4, tick.ratio=0.5) 
}
out <- dev.off()

pdf('Es2_MCMC_corp2.pdf')
par(mfrow = c(2,2))

# creating plot and computing results for mcmc about corp2
for (i in c(3:4)) {
  jump <- seq(from=1, to=nrow(Samples), by=25) # thin by factor 25
  filter_allSamp <- Samples[jump,i]
  
  nr <- length(filter_allSamp)
  
  l <- quantile(filter_allSamp, probs = 0.025)
  h <- quantile(filter_allSamp, probs = 0.975)
  media <- mean(filter_allSamp)
  med <- quantile(filter_allSamp, probs = 0.5)
  v <- var(filter_allSamp)
  
  l_limit <- c(l_limit, l)
  h_limit <- c(h_limit, h)
  median <- c(median, med)
  mean <- c(mean, media)
  var <- c(var, v)
  
  plot(c(1:nr), filter_allSamp,
       type="l",
       xlab="iteration",
       ylab='lambda', main=title[i])
  abline(h=media, lwd=1, lty=2, col='red')
  legend(x = "topright", legend = c("Mean"),
         lty = c(2), col = c('red'), lwd = 1, box.lty=0) 
  minor.tick(nx=4, ny=4, tick.ratio=0.5) 
  
  postDen <- density(filter_allSamp, n=10**3)
  
  plot(postDen$x, postDen$y, type='l', 
       col='navy', lwd = 2,
       xlab='lambda',
       ylab="p(lambda | data)", main=title[i] )
  polygon (c(l, postDen$x[postDen$x>l & postDen$x<h], h),
           c(0, postDen$y[postDen$x>l & postDen$x<h], 0), col=rgb(1, 0, 0,0.5))
  abline(v=media, lwd=1, lty=2)
  legend(x = "topright", legend = c("Mean"),
         lty = c(2), col = c('black'), lwd = 1, box.lty=0)
  minor.tick(nx=4, ny=4, tick.ratio=0.5) 
}
out <- dev.off()

summary2 <- tibble(name_prior=c('Corp1_Unif', 'Corp1_Jeff', 'Corp2_Unif', 'Corp2_Jeff'), 
                   mean=mean, median=median, var=var, l_limit=l_limit, h_limit=h_limit)


cat(' #########  ESERCIZIO 2   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('See Es2_MCMC_corp1.pdf to visualize the result ', '\n')
cat('See Es2_MCMC_corp2.pdf to visualize the result ', '\n', '\n')
cat('The following table summarises the obtained results:', '\n')
print(summary2)
cat('----------------------------------------------------------', '\n')	

                                                           ### ESERCIZIO 3 ###
# Experiment 1
n <- 116 
y <- 11
H0 <- 0.1

# Frequentist approach
freq_est1 <- y/n   # it has binomial distribution
x <- c(0:25)
p1 <- dbinom(x, n, H0)
cum <- cumsum(p1)

# compute lower and upper limit at 95%
l <- tail(x[cum<=0.025], 1)
u <- head(x[cum>=0.975], 1)
# setting color for the plot in order to have red tails
color <- rep('grey', length(x))
color[x<=l] <- 'red'
color[x>u] <- 'red'

# plot frequentist results 
pdf('Es3_exp1_F.pdf')
barplot(p1, name = c(0:25), col = color, main='Frequentist HT experiment I', sub='Since y=11 is inside the acceptance regior we can not reject H0' , xlab='y', ylab='prob')
minor.tick(nx=0, ny=4, tick.ratio=0.5)
text(x=14.5, y=0.03, "Acceptance", col='black', font=4, cex=2)
text(x=14.5, y=0.02, "region", col='black', font=4, cex=2)
out <- dev.off() 





# Bayesian approach
# using a beta prior Beta(1, 10), since the likelihood is a binomial, I know that the posterior is a beta
# with alpha'=y+alpha and beta'=n-y+beta
prob <- seq(0, 1, length=10**4)
posterior1 <- dbeta(prob, 1+y, n-y+10)

percentile <- function(posterior, parameter, perc) {
  cumulative <- cumsum(posterior/10**4)
  index <- which.min(abs(cumulative - perc))
  return(parameter[index])
}

mean_posterior1 <- sum(prob*posterior1)/10**4
var_posterior1 <- sum((prob-mean_posterior1)**2*posterior1)/10**4
lower_limit_1 <- percentile(posterior1, prob, 0.025)
upper_limit_1 <- percentile(posterior1, prob, 0.975)

pdf('Es3_exp1_B.pdf')
par(mfrow = c(2,1))
# plot posterior
plot(prob, posterior1, type='l', main='Posterior experiment I', xlab='p', ylab='prob(p | data)', xlim=c(0, 0.3))
polygon (c(lower_limit_1, prob[prob>lower_limit_1 & prob<upper_limit_1], upper_limit_1),
         c(0, posterior1[prob>lower_limit_1 & prob<upper_limit_1], 0), col=rgb(1, 0, 0,0.5))
abline(v = mean_posterior1, col='green', lty=2)
legend(x = "topright", legend = c("Posterior", '95% cred. int.', "Mean"),
       lty = c(1, 1, 2), col = c('black', 'red','green'), lwd = 2) 
minor.tick(nx=4, ny=5, tick.ratio=0.5)

# plot Bayesian Ht result
plot(prob, posterior1, type='l', main='Two sides HT experiment I', xlab='p', ylab='prob(p | data)', xlim=c(0, 0.3),
     sub='Since p=0.1 is inside the 95% regior we can not reject H0' )
polygon (c(lower_limit_1, prob[prob>lower_limit_1 & prob<upper_limit_1], upper_limit_1),
         c(0, posterior1[prob>lower_limit_1 & prob<upper_limit_1], 0), col=rgb(1, 0, 0,0.5))
abline(v = 0.1, col='blue', lty=1)
legend(x = "topright", legend = c("Posterior", '95% cred. int.', "H0"),
       lty = c(1, 1, 1), col = c('black', 'red','blue'), lwd = 2) 
minor.tick(nx=4, ny=5, tick.ratio=0.5)
out <- dev.off() 

# results all in a tibble
summary3_1 <- tibble(mean=mean_posterior1, var=var_posterior1, lowerlimit=lower_limit_1, upperlimit=upper_limit_1)


# Experiment 2
n <- 165 
y <- 9
H0 <- 0.1

# Frequentist approach
freq_est2 <- y/n   # it has binomial distribution
x <- c(0:35)
p2 <- dbinom(x, n, H0)
cum2 <- cumsum(p2)

# compute lower and upper limit at 95%
l <- tail(x[cum2<=0.025], 1)
u <- head(x[cum2>=0.975], 1)
# setting color for the plot in order to have red tails
color <- rep('grey', length(x))
color[x<=l] <- 'red'
color[x>u] <- 'red'


# plot frequentist results 
pdf('Es3_exp2_F.pdf')
barplot(p2, name = c(0:35), col = color, main='Frequentist HT experiment II', sub='Since y=9 is inside the acceptance regior we can not reject H0' , xlab='y', ylab='prob')
minor.tick(nx=0, ny=4, tick.ratio=0.5)
text(x=20, y=0.02, "Acceptance", col='black', font=4, cex=2)
text(x=20, y=0.01, "region", col='black', font=4, cex=2)
out <- dev.off() 


# Bayesian approach
# using a beta prior, since the likelihood is a binomial, I know that the posterior is a beta
# with alpha'=y+alpha and beta'=n-y+beta
prob <- seq(0, 1, length=10**4)

# posterior with as prior: Beta(1, 10)
posterior2_beta <- dbeta(prob, 1+y, n-y+10)
# posterior with as prior the posterior of the previous test (that is a beta)
posterior2_post <- dbeta(prob, 1+11+y, n-y+116-11+10)

percentile <- function(posterior, parameter, perc) {
  cumulative <- cumsum(posterior/10**4)
  index <- which.min(abs(cumulative - perc))
  return(parameter[index])
}

# result obtained using beta(1,10) prior
mean_posterior2_beta <- sum(prob*posterior2_beta)/10**4
var_posterior2_beta <- sum((prob-mean_posterior2_beta)**2*posterior2_beta)/10**4
lower_limit_2_beta <- percentile(posterior2_beta, prob, 0.025)
upper_limit_2_beta <- percentile(posterior2_beta, prob, 0.975)

# result obtained using the first exp posterior as prior
mean_posterior2_post <- sum(prob*posterior2_post)/10**4
var_posterior2_post <- sum((prob-mean_posterior2_post)**2*posterior2_post)/10**4
lower_limit_2_post <- percentile(posterior2_post, prob, 0.025)
upper_limit_2_post <- percentile(posterior2_post, prob, 0.975)

# all results in a tibble
summary3_2 <- tibble(prior=c('beta(1,10)', 'posterior'), mean=c(mean_posterior2_beta, mean_posterior2_post),
                     var=c(var_posterior2_beta, var_posterior2_post), lowerlimit=c(lower_limit_2_beta, lower_limit_2_post),
                     upperlimit=c(upper_limit_2_beta, upper_limit_2_post))



pdf('Es3_exp2_B.pdf')
par(mfrow = c(2,2))

plot(prob, posterior2_beta, type='l', main='Posterior experiment II', xlab='p', ylab='prob(p | data)', xlim=c(0, 0.2))
polygon (c(lower_limit_2_beta, prob[prob>lower_limit_2_beta & prob<upper_limit_2_beta], upper_limit_2_beta),
         c(0, posterior2_beta[prob>lower_limit_2_beta & prob<upper_limit_2_beta], 0), col=rgb(1, 0, 0,0.5))
abline(v = mean_posterior2_beta, col='green', lty=2, lwd = 2)
legend(x = "topright", legend = c("Posterior", '95% cred. int.', "Mean"),
       lty = c(1, 1, 2), col = c('black', 'red','green'), lwd = 2, title='Beta(1,10) as prior: ', cex=0.6) 
minor.tick(nx=4, ny=5, tick.ratio=0.5)

plot(prob, posterior2_beta, type='l', main='Two sides HT experiment II', xlab='p', ylab='prob(p | data)', xlim=c(0, 0.2),
     sub='p=0.1 is outside the 95% regior we can reject H0' )
polygon (c(lower_limit_2_beta, prob[prob>lower_limit_2_beta & prob<upper_limit_2_beta], upper_limit_2_beta),
         c(0, posterior2_beta[prob>lower_limit_2_beta & prob<upper_limit_2_beta], 0), col=rgb(1, 0, 0,0.5))
abline(v = 0.1, col='blue', lty=1)
legend(x = "topright", legend = c("Posterior", '95% cred. int.', "H0"),
       lty = c(1, 1, 1), col = c('black', 'red','blue'), lwd = 2, title='Beta(1,10) as prior: ', cex=0.6)
minor.tick(nx=4, ny=5, tick.ratio=0.5)

plot(prob, posterior2_post, type='l', main='Posterior experiment II', xlab='p', ylab='prob(p | data)', xlim=c(0, 0.2))
polygon (c(lower_limit_2_post, prob[prob>lower_limit_2_post & prob<upper_limit_2_post], upper_limit_2_post),
         c(0, posterior2_post[prob>lower_limit_2_post & prob<upper_limit_2_post], 0), col=rgb(1, 0, 0,0.5))
abline(v = mean_posterior2_post, col='green', lty=2, lwd = 2)
legend(x = "topright", legend = c("Posterior", '95% cred. int.', "Mean"),
       lty = c(1, 1, 2), col = c('black', 'red','green'), lwd = 2, title='Exp1 post as prior: ', cex=0.6) 
minor.tick(nx=4, ny=5, tick.ratio=0.5)

plot(prob, posterior2_post, type='l', main='Two sides HT experiment II', xlab='p', ylab='prob(p | data)', xlim=c(0, 0.2),
     sub='p=0.1 is inside the 95% regior we can not reject H0' )
polygon (c(lower_limit_2_post, prob[prob>lower_limit_2_post & prob<upper_limit_2_post], upper_limit_2_post),
         c(0, posterior2_post[prob>lower_limit_2_post & prob<upper_limit_2_post], 0), col=rgb(1, 0, 0,0.5))
abline(v = 0.1, col='blue', lty=1)
legend(x = "topright", legend = c("Posterior", '95% cred. int.', "H0"),
       lty = c(1, 1, 1), col = c('black', 'red','blue'), lwd = 2, title='Exp1 post as prior: ', cex=0.6)
minor.tick(nx=4, ny=5, tick.ratio=0.5)
out <- dev.off() 

cat(' #########  ESERCIZIO 3   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('              Experiment I ', '\n', '\n')
cat('The frequentist approach can be visualized in: Es3_exp1_F.pdf', '\n')
cat('The Bayesian approach can be visualized in: Es3_exp1_B.pdf', '\n')
cat('Bayesian results: ', '\n')
print(summary3_1)

cat('\n', '\n')
cat('              Experiment II ', '\n', '\n')
cat('The frequentist approach can be visualized in: Es3_exp2_F.pdf', '\n')
cat('The Bayesian approach can be visualized in: Es3_exp2_B.pdf', '\n')
cat('Bayesian results: ', '\n')
print(summary3_2)

cat('----------------------------------------------------------', '\n')	

                                                           ### ESERCIZIO 4 ###
cat(' #########  ESERCIZIO 4   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('.bug file used to create the model:  ', '\n', '\n')
cat('model {', '\n', '# data likelihood' ,'\n', 'X ~ dbinom(p, Tot);','\n' ,'\n')
cat(' # a uniform prior for lambda', '\n', 'p ~ dbeta(1, 10);' ,'\n', '}' ,'\n','\n')
cat('----------------------------------------------------------', '\n')
cat('The plots about the reslts of the mcmc are inside: Es4_chain.pdf', '\n')
cat('The inference on p is shown in: Es4_Gibbs_hist.pdf', '\n', '\n')
cat('Informations about the mcmc and the resulting posterior : ', '\n', '\n')

data <- NULL
data$X <- 11 # number of counts
data$Tot <- 116

# contenuto file 'Model_es_4.5.bug':
#model {
#  # data likelihood
#  X ~ dbinom(p, Tot);

#  # a uniform prior for lambda
#  p ~ dbeta(1, 10);
#}

# import the model
model <- "Model_es_4.5.bug"
jm <- jags.model(model , data)

# burn-in
update(jm , 1000)

# create mcmc 
chain <- coda.samples(jm , c("p"), n.iter=10000)

pdf('Es4_chain.pdf')
plot(chain , col="navy")
out <- dev.off() 

chain.df <- as.data.frame(as.mcmc(chain))

print(summary(chain))

pdf('Es4_Gibbs_hist.pdf')
mean <- paste('Mean =  ', round(summary(chain)$statistics[1], digits = 3))
Sd <- paste('Sd =  ', round(summary(chain)$statistics[2], digits = 3))
ll <- paste('2.5% limit =  ', round(summary(chain)$quantiles[1], digits = 3))
ul <- paste('97.5% limit =  ', round(summary(chain)$quantiles[5], digits = 3))
hist(chain.df$p , nc=50, prob=TRUE , col='darkolivegreen2',
     xlab='p', ylab='prob(p | data)', main='Inference on p')
minor.tick(nx=4, ny=5, tick.ratio=0.5)
legend(x = "topright", legend = c(mean, Sd, ll, ul),
       lty = c(0, 0, 0,  0), col = c('', '','',''), lwd = 2, title='Summary: ', cex=0.6) 
out <- dev.off() 

cat('----------------------------------------------------------', '\n')	