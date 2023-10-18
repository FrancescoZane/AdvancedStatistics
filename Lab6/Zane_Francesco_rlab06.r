library(tibble)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(mvtnorm)
library(coda)
library(rjags)
library(dplyr)
library(lubridate)
library(tsibble) 
library(scales)

                                                                     ### ESERCIZIO 1 ###

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
      #cat(s1, "of", s2, "+", s3, s4, "\n")
    }
    
  }
  
  if(demo) {
    return(list(funcSamp=funcSamp, thetaPropAll=thetaPropAll))
  } else {
    return(funcSamp)
  }
}

# def function I want to sample
f <- function(x) {
  fun <- (0.5 * exp(-((x+3)**2)/2) + 0.5 * exp(-((x-3)**2)/2))
  return(fun)
}

# function as input for metro()
f_input <- function(x) {
  return(c(log10(f(x)), 0))
}

pdf('Es1_6_Nburnin.pdf')
# study the variation of the burn-in
par(mfrow = c(3,2))
for (i in c(0, 500, 2000)) {
  
  set.seed(12345)
  ThetaZero <- runif(1, -10, 10)           # random initial value
  sampleCov <- 1                           # fix sigma for sampling 
  allSamp <- metrop(func=f_input, thetaInit=ThetaZero, Nburnin=i, Nsamp=10**5, sampleCov=sampleCov , verbose=1e3)
  cat(str(allSamp))                        # run mcmc
  
  jump <- seq(from=1, to=nrow(allSamp), by=25)         # thin factor
  filter_allSamp <- allSamp[jump, 3]
  nr <- length(filter_allSamp)
  
  plot(c(1:nr), filter_allSamp,                        # plot the resulting sampled values of mcmc
       xlab="iterations",
       ylab='x', main=paste('Nburnin =',i), type='l')
  minor.tick(nx=4, ny=4, tick.ratio=0.5) 
  
  postDen <- density(filter_allSamp, n=300)
  x <- seq(-8, 8, length=10**4)
  normalization <- sum(f(x)*(16/10**4))               # calculate in order to compare the hist of the mcmc with the function
  
  plot(postDen$x, postDen$y, type='s',                # plot the 'hist' of the resulting sampled values of mcmc
       col='navy', lwd = 2,
       xlab='x',
       ylab="f(x)",  ylim=c(0, 0.28))
  minor.tick(nx=4, ny=4, tick.ratio=0.5)
  lines(x, f(x)/normalization, col='red')
  legend(x = "topright", legend = c("Function  f(x)", 'Density of sampled data'),
         lty = c(1, 1), col = c('red','blue'), lwd = 2, cex=0.6) 
}
o <- dev.off()

pdf('Es1_6_thin.pdf')
# study the variation of the thin factor
par(mfrow = c(4,3))
for (i in c(1,10,20,100)) {
  
  set.seed(1234)
  ThetaZero <- runif(1, -10, 10)              # random initial value
  sampleCov <- 1                              # fix sigma for sampling 
  allSamp <- metrop(func=f_input, thetaInit=ThetaZero, Nburnin=200, Nsamp=10**4*i, sampleCov=sampleCov , verbose=1e3)      # run mcmc
  # Nsamp proportional to 'i', in this way all chain has the same length
  
  jump <- seq(from=1, to=nrow(allSamp), by=i) # thin by factor i
  filter_allSamp <- allSamp[jump, 3]
  nr <- length(filter_allSamp)
  
  chain <- as.mcmc(filter_allSamp)
  l = seq(0,300,5)
  ACF <- autocorr(chain, lags=l)              # compute ACF at different steps of the mcmc
  
  plot(c(1:200), filter_allSamp[1:200],       # plot the fisrt 200 resulting sampled values of mcmc
       xlab="iterations",
       ylab='x', main=paste('thin by factor', i), ylim=c(-6, 6))
  lines(c(1:200), filter_allSamp[1:200], col='red')
  minor.tick(nx=4, ny=4, tick.ratio=0.5) 
  
  plot(l[1:20] , ACF[1:20], ylim=c(0,1),      # plotting the result of the ACF   
       pch=12, col='navy',
       xlab='chain length', ylab='ACF', cex=1.3)
  text(400,0.9, paste('Each', i))
  text(400,0.85,
       sprintf("effective size: %.2f",
               effectiveSize(chain)))
  
  postDen <- density(filter_allSamp, n=300)
  x <- seq(-8, 8, length=10**4)
  normalization <- sum(f(x)*(16/10**4))      # calculate in order to compare the hist of the mcmc with the function
  
  plot(postDen$x, postDen$y, type='s',       # plot the 'hist' of the resulting sampled values of mcmc
       col='navy', lwd = 2,
       xlab='x',
       ylab="f(x)",  ylim=c(0, 0.22), xlim=c(-9, 9))
  minor.tick(nx=4, ny=4, tick.ratio=0.5)
  lines(x, f(x)/normalization, col='red')
  legend(x='topright', legend = c("f(x)", 'mcmc'),
         lty = c(1, 1), col = c('red','blue'), lwd = 2, cex=0.6) 
}
o <- dev.off()

cat('\n', '\n')
cat(' #########  ESERCIZIO 1   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat(' Visualize the result varing the burnin: Es1_6_Nburnin.pdf', '\n', '\n')
cat(' Here we can notice that there are no differences. ', '\n')
cat(' In this case the burnin is useless since just after few iterations' , '\n', 'the algorithm is already in the interested region.', '\n')
cat('----------------------------------------------------------', '\n')	
cat(' Visualize the result varing the thin factor: Es1_6_thin.pdf', '\n', '\n')
cat(' Here we can notice that varing the thin factorn the ACF value,','\n','among the iterations, goes down more rapidly. ', '\n')
cat(' The plots on the left allow us to notice that increasing the thin', '\n')
cat(' factor the dependences among the sampled values of the mcmc decrease. ', '\n')
cat('----------------------------------------------------------', '\n', '\n')	


                                                                     ### ESERCIZIO 2 ###

cat(' #########  ESERCIZIO 2  #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('.bug file used to create the model:  ', '\n', '\n')
cat('model {', '\n', '# data likelihood' ,'\n', 'X ~ dbinom(p, Tot);','\n' ,'\n')
cat(' # a uniform prior for lambda', '\n', 'p ~ dbeta(1, 1);' ,'\n', '}' ,'\n','\n')
cat('----------------------------------------------------------', '\n')

set.seed(123456)

# contenuto file 'Model_es_2.6.bug':
#model {
#  # data likelihood
#  X ~ dbinom(p, Tot);
  
#  # a uniform prior for lambda
#  p ~ dbeta(1,1);
#}


# def function to analyse the vaccine data
analysis <- function(vaccine, placebo, name_vaccine) {
  
  # vaccine data
  data_v <- c()
  data_v$X <- vaccine[1]
  data_v$Tot <- vaccine[2]
  
  # import jags model
  model <- "Model_es_2.6.bug"
  jm <- jags.model(model , data_v)
  
  # run the burnin and the chain
  update(jm , 1000)
  chain_v <- coda.samples(jm , c("p"), n.iter=10000)
  chain_v.df <- as.data.frame(as.mcmc(chain_v))
  
  mean_v <- paste('Mean =  ', round(summary(chain_v)$statistics[1], digits = 3))
  Sd_v <- paste('Sd =  ', round(summary(chain_v)$statistics[2], digits = 3))
  ll_v <- paste('2.5% limit =  ', round(summary(chain_v)$quantiles[1], digits = 3))
  ul_v <- paste('97.5% limit =  ', round(summary(chain_v)$quantiles[5], digits = 3))
  
  # placebo data
  data_p <- c()
  data_p$X <- placebo[1]
  data_p$Tot <- placebo[2]
  
  # import jags model
  model <- "Model_es_2.6.bug"
  jm <- jags.model(model , data_p)
  
  # run the burnin and the chain
  update(jm , 1000)
  chain_p <- coda.samples(jm , c("p"), n.iter=10000)
  chain_p.df <- as.data.frame(as.mcmc(chain_p))
  
  mean_p <- paste('Mean =  ', round(summary(chain_p)$statistics[1], digits = 3))
  Sd_p <- paste('Sd =  ', round(summary(chain_p)$statistics[2], digits = 3))
  ll_p <- paste('2.5% limit =  ', round(summary(chain_p)$quantiles[1], digits = 3))
  ul_p <- paste('97.5% limit =  ', round(summary(chain_p)$quantiles[5], digits = 3))
  
  # diff rate 
  diff <- (1 - (chain_v.df$p/chain_p.df$p))*100
  postDen <- density(diff, n=500)
  
  media <- mean(diff)
  med <- quantile(diff, probs = 0.5)
  v <- var(diff)
  l <- quantile(diff, probs = 0.025)
  u <- quantile(diff, probs = 0.975)
  
  # plot the results
  pdf(paste(name_vaccine, '.pdf', sep=''))
  par(mfrow = c(3,1))
  
  # vaccine
  hist(chain_v.df$p, freq=FALSE, breaks=50, col='darkolivegreen2',
       xlab='percentage infected', ylab='density', main=paste(name_vaccine, 'vaccine'))
  minor.tick(nx=4, ny=5, tick.ratio=0.5)
  legend(x = "topright", legend = c(mean_v, Sd_v, ll_v, ul_v),
         lty = c(0, 0, 0,  0), col = c('', '','',''), lwd = 2, title='Summary: ', cex=0.6)
  
  #placebo
  hist(chain_p.df$p, freq=FALSE, breaks=50, col='darkolivegreen2',
       xlab='percentage infected', ylab='density', main=paste(name_vaccine, 'placebo'))
  minor.tick(nx=4, ny=5, tick.ratio=0.5)
  legend(x = "topright", legend = c(mean_p, Sd_p, ll_p, ul_p),
         lty = c(0, 0, 0,  0), col = c('', '','',''), lwd = 2, title='Summary: ', cex=0.6)
  
  # effectiveness
  plot(postDen$x, postDen$y, type='s', 
       col='navy', lwd = 2,
       xlab='vaccine effectiveness (%)',
       ylab="probability",
       main=paste('pdf', name_vaccine, 'vaccine effectiveness'),
       sub=paste('This study shows an effectiveness = ', round(media, digits=1), '±',
                 round(v, digits=1), '            (gauss app.)'))
  minor.tick(nx=5, ny=4, tick.ratio=0.5)
  legend(x = "topright", legend = c(paste('Mean =  ', round(media, digits=3)), 
                                    paste('Median = ', round(med, digits=3)),
                                    paste('Sd =  ', round(v, digits = 3)),
                                    paste('2.5% limit =  ', round(l, digits = 3)),
                                    paste('97.5% limit =  ', round(u, digits = 3))), 
         lty = c(2, 3, 0, 0, 0), col = c('black', 'black','','', ''), lwd = 2, title='Summary: ', cex=0.6)
  polygon (c(l, postDen$x[postDen$x>l & postDen$x<u], u),
           c(0, postDen$y[postDen$x>l & postDen$x<u], 0), col=rgb(1, 0, 0,0.5))
  abline(v = media, col='black', lty=2)
  abline(v = med, col='black', lty=3)
  o <- dev.off()
}

# Run analysis for Jcovden
analysis(c(116, 19630), c(348, 19691), 'Jcovden') 

# Run analysis for Spikevax
analysis(c(11, 14134), c(185, 14073), 'Spikevax') 

# Run analysis for Vaxzevria (I test)
analysis(c(64, 5258), c(154, 5210), 'Vaxzevria (I test)')

# Run analysis for Vaxzevria (II test)
analysis(c(73, 17662), c(130, 8550), 'Vaxzevria (II test)') 

cat('----------------------------------------------------------', '\n')	
cat('Vaccine Jcovden : Jcovden.pdf', '\n')
cat('Vaccine Spikevax : Spikevax.pdf', '\n')
cat('Vaccine Vaxzevria (I test) : Vaxzevria (I test).pdf', '\n')
cat('Vaccine Vaxzevria (II test) : Vaxzevria (II test).pdf.pdf', '\n')
cat('----------------------------------------------------------', '\n', '\n')	



                                                                     ### ESERCIZIO 3 ###

set.seed(123456)

World = read.csv('owid-covid-data.csv', header=TRUE)

Analysis <- function(place, cont) {

  if (place=='') {

    # data concerning a continent

    vaccine <- World[c('date', 'continent', 'new_vaccinations')] %>% filter(continent==cont) %>% filter(is.na(new_vaccinations)==FALSE) %>%
      group_by(date) %>% summarize(date=first(date), new_vaccinations=sum(new_vaccinations))
    vaccine['people_vaccinated'] <- cumsum(vaccine$new_vaccinations)
    vaccine['Week'] <- yearweek(as.Date(vaccine$date))

    death <- World[c('date', 'continent', 'new_deaths')] %>% filter(continent==cont) %>%filter(is.na(new_deaths)==FALSE) %>%
      group_by(date) %>% summarize(date=first(date), new_deaths=sum(new_deaths))
    death['total_deaths'] <- cumsum(death$new_deaths)
    death['Week'] <- yearweek(as.Date(death$date))

    death['Week'] <- yearweek(as.Date(death$date))
    area <- cont
  } else {

    # data concerning a single state

    vaccine <- World[c('date', 'location', 'new_vaccinations', 'people_vaccinated')] %>% filter(location==place)  %>% filter(is.na(people_vaccinated)==FALSE & is.na(new_vaccinations)==FALSE)
    vaccine['Week'] <- yearweek(as.Date(vaccine$date))
    death <- World[c('date', 'location', 'total_deaths', 'new_deaths')] %>% filter(location==place) %>% filter(is.na(total_deaths)==FALSE & is.na(new_deaths)==FALSE)
    death['Week'] <- yearweek(as.Date(death$date))
    area <- place
  }

  # compute the weekly average
  vacc_week <- vaccine %>% group_by(Week) %>% summarize(week_vacc=mean(new_vaccinations), date=first(date))
  death_week <- death %>% group_by(Week) %>% summarize(week_death=mean(new_deaths), date=first(date))

  # plot of cumulative death
  tot_death <- ggplot(death, aes(x=as.Date(date, "%Y-%m-%d"), y=total_deaths)) +
    geom_line() +
    labs(x='Date', y='', title=paste(area, 'cumulative number of death people')) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) + scale_x_date(date_breaks = "2 month", date_labels = "%b/%y") +
    scale_y_continuous(labels = label_number(scale_cut=cut_si('') ), breaks=seq(0, max(death$total_deaths), length=10))

  # plot of daily death
  new_death_day <- ggplot(death, aes(x=as.Date(date, "%Y-%m-%d"), y=new_deaths)) +
    geom_line() +
    labs(x='date', y='', title=paste(area, 'daily number of death people')) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) + scale_x_date(date_breaks = "2 month", date_labels = "%b/%y") +
    scale_y_continuous(labels = label_number(scale_cut=cut_si('') ), breaks=seq(0, max(death$new_deaths), length=10))

  # plot of daily death (week average)
  new_death_week <- ggplot(death_week, aes(x=as.Date(date, "%Y-%m-%d"), y=week_death)) +
    geom_line() +
    labs(x='date', y='', title=paste(area, 'weekly average number of death people')) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) + scale_x_date(date_breaks = "2 month", date_labels = "%b/%y") +
    scale_y_continuous(labels = label_number(scale_cut=cut_si('') ), breaks=seq(0, max(death_week$week_death), length=10))

  pdf(paste(area, '_death.pdf', sep=''))
  grid.arrange(new_death_day, new_death_week, tot_death, nrow = 3)
  o <- dev.off()

  # plot of cumulative vaccined people
  tot_vacc <- ggplot(vaccine, aes(x=as.Date(date, "%Y-%m-%d"), y=people_vaccinated)) +
    geom_line() +
    labs(x='Date', y='', title=paste(area, 'cumulative number of vaccinated people')) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) + scale_x_date(date_breaks = "2 month", date_labels = "%b/%y") +
    scale_y_continuous(labels = label_number(scale_cut=cut_si('') ), breaks=seq(0, max(vaccine$people_vaccinated), length=10))

  # plot of daily vaccined people
  new_vacc_day <- ggplot(vaccine, aes(x=as.Date(date, "%Y-%m-%d"), y=new_vaccinations)) +
    geom_line() +
    labs(x='date', y='', title=paste(area, 'daily number of vaccinated people')) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) + scale_x_date(date_breaks = "2 month", date_labels = "%b/%y") +
    scale_y_continuous(labels = label_number(scale_cut=cut_si('') ), breaks=seq(0, max(vaccine$new_vaccinations), length=10))

  # plot of daily vaccined people (week average)
  new_vacc_week <- ggplot(vacc_week, aes(x=as.Date(date, "%Y-%m-%d"), y=week_vacc)) +
    geom_line() +
    labs(x='date', y='', title=paste(area, 'weekly average number of vaccinated people')) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) + scale_x_date(date_breaks = "2 month", date_labels = "%b/%y") +
    scale_y_continuous(labels = label_number(scale_cut=cut_si('') ), breaks=seq(0, max(vacc_week$week_vacc), length=10))

  pdf(paste(area, '_vaccine.pdf', sep=''))
  grid.arrange(new_vacc_day, new_vacc_week, tot_vacc, nrow = 3)
  o <- dev.off()
}

Analysis('', 'North America')
Analysis('', 'Asia')
Analysis('', 'Europe')
Analysis('France', 'Europe')
Analysis('Italy', 'Europe')

cat(' #########  ESERCIZIO 3  #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('Vaccine data: ', '\n', '\n')
cat('Asia : Asia_vaccine.pdf', '\n')
cat('North America : North America_vaccine.pdf', '\n')
cat('Europe : Europe_vaccine.pdf', '\n')
cat('France : France_vaccine.pdf', '\n')
cat('Italy : Italy_vaccine.pdf', '\n')
cat('----------------------------------------------------------', '\n')	
cat('Death data: ', '\n', '\n')
cat('Asia : Asia_death.pdf', '\n')
cat('North America : North America_death.pdf', '\n')
cat('Europe : Europe_death.pdf', '\n')
cat('France : France_death.pdf', '\n')
cat('Italy : Italy_death.pdf', '\n')
cat('----------------------------------------------------------', '\n', '\n')	



