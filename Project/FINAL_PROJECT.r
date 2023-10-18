# In the same folder of this code other folders are needed

# Folders in which the code takes information:
# - Data : contains the .csv file 
# - bug : contains the bug file

# Folders in which the code saves graphs:
# - Month
# - Constant_model
# - Linear_model
# - ModelComparison
# - PlaceAnalysis
# - PlaceComparison
# - PaperComparison
# - Arima 

# path up to the folder that contains the code
path <- '/home/francesco/AdvaceStatistics/Project/'

library(Hmisc)
library(dplyr)
library(gplots)
library(rjags)
library(coda)
library(scales)
library(forecast)
library(tseries)
library(lubridate)
library(qpdf)

                                                              ### plot the month

# this function return a pdf with five plots, each plot represents the behavior of the selected month over the years for each place
plot_month <- function(mese) {
  stations <- c('Auronzo_2m', 'Roverchiara_2m', 'PortoTolle_2m', 'Auronzo', 'Castelfranco')
  
  name_month <- switch(mese, 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEPT', 'OCT', 'NOV', 'DEC')

  pdf(paste(path, '/Month/' ,name_month, '.pdf', sep=''))
  par(mfrow = c(2,3))
  
  for (place in stations) {
    
    # import file cancelling na rows
    file <- read.csv(paste(path, 'Data/',place,'.csv',sep=''), header=TRUE, sep=',')
    file <- file |> filter( !is.na(TMIN), !is.na(TMED), !is.na(TMAX) )
    
    # convert in dates
    file$date <- as.Date(file$date)
    file$month <- month(file$date)
    file$year <- year(file$date)
    
    # month average
    file_month <- file |> group_by(year, month) |> summarise(Tmin=mean(TMIN), Tmed=mean(TMED), Tmax=mean(TMAX))
    
    # select the month we want to study
    Temperature <- file_month |> filter(month==mese)
    
    # plot
    plot(Temperature$year, Temperature$Tmed, xlab='Years', ylab='Temperature [°C]', main=paste(place, '  -  ', name_month, sep=''),
         type='l', col='green',
         ylim=c(min(Temperature$Tmin - 1), max(Temperature$Tmax + 4)))
    lines(Temperature$year, Temperature$Tmin, col='blue')
    lines(Temperature$year, Temperature$Tmax, col='red')
    grid()
    
    legend("topright", 
           legend = c('Average', 'Minimum', 'Maximum'),
           lty = c(1, 1, 1), col = c('green', 'blue', 'red'), lwd = 2, title='Legend', cex=0.5)
    
    minor.tick(nx=5, ny=5, tick.ratio=0.5)
  }
  
  o <- dev.off()
}

# choose the month you want to study
m <- 12
plot_month(m)


### Prepare the data for the fit

stations <- c('Auronzo', 'Auronzo_2m', 'Castelfranco', 'PortoTolle_2m', 'Roverchiara_2m')

for (place in stations) {

  # import file cancelling na rows
  file <- read.csv(paste(path, 'Data/',place,'.csv',sep=''), header=TRUE, sep=',')
  file <- file |> filter( !is.na(TMIN), !is.na(TMED), !is.na(TMAX) )
  
  # convert in dates
  file$date <- as.Date(file$date)
  file$month <- month(file$date)
  file$Anno <- year(file$date)
  
  # year average
  file_min <- data.frame(Anno=file$Anno, Tmin=file$TMIN) |> group_by(Anno) |> summarise(Medio=mean(Tmin))
  file_med <- data.frame(Anno=file$Anno, Tmed=file$TMED) |> group_by(Anno) |> summarise(Medio=mean(Tmed))
  file_max <- data.frame(Anno=file$Anno, Tmax=file$TMAX) |> group_by(Anno) |> summarise(Medio=mean(Tmax))
  
  # create columns 
  file_min$info <- 'min'
  file_med$info <- 'ave'
  file_max$info <- 'max'
  
  file_min <- file_min[-1,]
  file_max <- file_max[-1,]
  file_med <- file_med[-1,]
  
  # final data frame
  final_file <- rbind(file_min, file_med, file_max)
  comment(final_file) <- place
  assign(place , final_file)
  
}


# This function returns 5 pdfs
#  - jags linear model
#  - linear model and cov
#  - jags constant model
#  - constant model and cov
#  - comparison between the two models

# A pdf file for each of them and a unique pdf file as summary
Analysis <- function(place, Temp) {
  
  # offset of the data
  min_years <- min(place$Anno)
  place$Anno <- place$Anno - min_years
  
  # Gibbs sampling from linear model
  model_l <- "/home/francesco/AdvaceStatistics/Project/bug/LinearModel.bug"
  jm_l <- jags.model(model_l, filter(place[, c('Anno', 'Medio', 'info')], info==Temp))
  update(jm_l , 1000)
  chain_l <- coda.samples(jm_l , c("Bf", "alpha", "beta", "sigma"), n.iter=20000, n.thin = 20)
  
  # Gibbs sampling from constant model
  model_c <- "/home/francesco/AdvaceStatistics/Project/bug/ConstantModel.bug"
  jm_c <- jags.model(model_c, filter(place[, c('Medio', 'info')], info==Temp))
  update(jm_c , 1000)
  chain_c <- coda.samples(jm_c , c("Bf", "alpha", "sigma"), n.iter=20000, n.thin = 20)
  
  # plot the two chains
  pdf(paste(path, 'Linear_model/JAGS_linearmodel_', comment(place), '_', Temp, '.pdf', sep=''))
  plot(chain_l)
  o <- dev.off()
  
  pdf(paste(path, 'Constant_model/JAGS_constantmodel_', comment(place), '_', Temp, '.pdf', sep=''))
  plot(chain_c)
  o <- dev.off()
  Bayes_factor_final <- c()
  result_fit_final <- data.frame()
  # take fit values
  alpha_l <- summary(chain_l)$statistics[2,1]
  beta_l <- summary(chain_l)$statistics[3,1]
  sd_l <- summary(chain_l)$statistics[4,1]
  Bf_l <- summary(chain_l)$statistics[1,1]
  
  alpha_c <- summary(chain_c)$statistics[2,1]
  sd_c <- summary(chain_c)$statistics[3,1]
  Bf_c <- summary(chain_c)$statistics[1,1]
  
  years <- filter(place, info==Temp)$'Anno' + min_years

  
  
  
                                              # plot linear model & correlations
  
  # correlations for linear model
  chain_l_df <- as.data.frame( as.mcmc(chain_l) )
  # extract chains
  chain_alpha_l <- chain_l_df[, 2]
  chain_beta_l <- chain_l_df[, 3]
  chain_sd_l <-  chain_l_df[, 4]
  # calculate correlations
  cor_alpha_beta_l <- cor(chain_alpha_l, chain_beta_l)
  cor_alpha_sigma_l <- cor(chain_alpha_l, chain_sd_l)
  cor_beta_sigma_l <- cor(chain_beta_l, chain_sd_l)
  
  
  grDevices::cairo_pdf(paste(path, 'Linear_model/LinearModel_', comment(place), '_', Temp, '.pdf', sep=''))
  area <- matrix( c( 1, 1, 1, 2, 3, 4), nrow=2, byrow=TRUE)
  layout(area)
  
      # linear model
  plotCI(x=years, y=filter(place, info==Temp)$Medio, uiw=sd_l, xlab='years', ylab='Temperature °C',  
         main=paste(comment(place), Temp, '  -', '  Linear model'))
  lines(years, alpha_l + beta_l * (years-min_years), col='red')
  minor.tick(nx=5, ny=2, tick.ratio=0.5)
  legend("topright", inset=c(0, -0.18),
         legend = c(paste('α [ °C ]:  ', round(alpha_l, 3), '±', round(summary(chain_l)$statistics[2,2], 3)),
                    paste('β [ °C / y ]:  ', round(beta_l, 3), '±', round(summary(chain_l)$statistics[3,2], 3)),
                    paste('σ [ °C  ]:  ', round(sd_l, 3), '±', round(summary(chain_l)$statistics[4,2], 3))),
         lty = c(0, 0, 0), col = c('', '', ''), lwd = 2,
         title=paste('Model:  α + β*(x - ', min_years, ' )', sep=''), cex=0.6, xpd=TRUE)
  grid()
  
       # chain variables and correlations
  plot(chain_alpha_l, chain_beta_l, pch=20, col=alpha('black', 0.1), xlab='α', ylab='β', main='cor(α,β)')
  legend("topright", legend = c('cor(α,β) :  ', round(cor_alpha_beta_l, 3)), cex=0.6)
  plot(chain_alpha_l, chain_sd_l, col=alpha('black', 0.1), xlab='α', ylab='σ', main='cor(α,σ)')
  legend("topright", legend = c('cor(α,σ) :  ', round(cor_alpha_sigma_l, 3)), cex=0.6)
  plot(chain_beta_l, chain_sd_l, col=alpha('black', 0.1), xlab='β', ylab='σ', main='cor(β,σ)')
  legend("topright", legend = c('cor(β,σ) :  ', round(cor_beta_sigma_l, 3)), cex=0.6)
  
  o <- dev.off()

  
  
  
                                              # plot constant model & correlations
  
  # correlations for linear model
  chain_c_df <- as.data.frame( as.mcmc(chain_c) )
  # extract chains
  chain_alpha_c <- chain_c_df[, 2]
  chain_sd_c <- chain_l_df[, 3]
  # calculate correlations
  cor_alpha_sigma_c <- cor(chain_alpha_c, chain_sd_c)
  
  
  grDevices::cairo_pdf(paste(path, 'Constant_model/ConstantModel_', comment(place), '_', Temp, '.pdf', sep=''))
  area <- matrix( c( 1, 1, 1, 1, 1, 0, 2, 2, 2, 0), nrow=2, byrow=TRUE)
  layout(area)
  
      # constant model
  plotCI(x=years, y=filter(place, info==Temp)$Medio, uiw=sd_c, xlab='years', ylab='Temperature °C', 
         main=paste(comment(place), Temp, '  -', '  Constant model'))
  lines(years, rep(alpha_c, length(years)), col='red')
  minor.tick(nx=5, ny=2, tick.ratio=0.5)
  legend("topright", inset=c(0, -0.145),
         legend = c(paste('α [ °C ]:  ', round(alpha_c, 3), '±', round(summary(chain_c)$statistics[2,2], 3)),
                    paste('σ [ °C ]:  ', round(sd_c, 3), '±', round(summary(chain_c)$statistics[3,2], 3))),
         lty = c(0, 0), col = c('', ''), lwd = 2, title='Model:  α ', cex=0.6, xpd=TRUE)
  grid()
  
      # chain variables and correlations
  plot(chain_alpha_c, chain_sd_c, pch=20, col=alpha('black', 0.1), xlab='α', ylab='σ', main='cor(α,σ)')
  legend("topright", legend = c('cor(α,σ) :  ', round(cor_alpha_sigma_c, 3)), cex=0.6)
  
  o <- dev.off()
  
  
  
  # plot comparison between the two models
  
  grDevices::cairo_pdf(paste(path, 'ModelComparison/ModelComparison_', comment(place), '_', Temp, '.pdf', sep=''))
  par(mfrow = c(2,1))
  
  # linear model
  plotCI(x=years, y=filter(place, info==Temp)$Medio, uiw=sd_l, xlab='years', ylab='Temperature °C',  
         main=paste(comment(place), Temp, '  -', '  Linear'))
  lines(years, alpha_l + beta_l * (years-min_years), col='red')
  minor.tick(nx=5, ny=2, tick.ratio=0.5)
  legend("topright", inset=c(0, -0.31),
         legend = c(paste('α [ °C ]  :  ', round(alpha_l, 3), '±', round(summary(chain_l)$statistics[2,2], 3)),
                    paste('β [ °C / y ]:  ', round(beta_l, 3), '±', round(summary(chain_l)$statistics[3,2], 3)),
                    paste('σ [ °C ]  :  ', round(sd_l, 3), '±', round(summary(chain_l)$statistics[4,2], 3))),
         lty = c(0, 0, 0), col = c('', '', ''), lwd = 2,
         title=paste('Model:  α + β*(x - ', min_years, ' )', sep=''), cex=0.5, xpd=TRUE)
  grid()
  
  # constant model
  plotCI(x=years, y=filter(place, info==Temp)$Medio, uiw=sd_c, xlab='years', ylab='Temperature °C', 
         main=paste(comment(place), Temp, '  -', '  Constant'))
  lines(years, rep(alpha_c, length(years)), col='red')
  minor.tick(nx=5, ny=2, tick.ratio=0.5)
  legend("topright", inset=c(0, -0.25),
         legend = c(paste('α [ °C ]:  ', round(alpha_c, 3), '±', round(summary(chain_c)$statistics[2,2], 3)),
                    paste('σ [ °C ]:  ', round(sd_c, 3), '±', round(summary(chain_c)$statistics[3,2], 3))),
         lty = c(0, 0), col = c('', ''), lwd = 2, title='Model:  α ', cex=0.5, xpd=TRUE)
  grid()
  o <- dev.off()

  
  pdf_combine(input = c(paste(path, 'Linear_model/JAGS_linearmodel_', comment(place), '_', Temp, '.pdf', sep=''),
                        paste(path, 'Linear_model/LinearModel_', comment(place), '_', Temp, '.pdf', sep=''),
                        paste(path, 'Constant_model/JAGS_constantmodel_', comment(place), '_', Temp,'.pdf', sep=''),
                        paste(path, 'Constant_model/ConstantModel_', comment(place), '_', Temp, '.pdf', sep=''),
                        paste(path, 'ModelComparison/ModelComparison_', comment(place), '_', Temp, '.pdf', sep='')), 
              output = paste(path, 'PlaceAnalysis/', comment(place), '_', Temp, '.pdf', sep=''))
  
  
  
  return( list(Bf_l/Bf_c, paste(comment(place), '_', Temp, sep=''), beta_l, summary(chain_l)$statistics[3,2], Temp))
}

Bayes_factor <- c()
result_fit <- data.frame()

                                                                          # fits 
for (zone in list(Auronzo, Auronzo_2m, Castelfranco, PortoTolle_2m, Roverchiara_2m)) {
  for (t in c('min', 'max', 'ave')) {
    v <- Analysis(zone, t)
    Bayes_factor <- c(Bayes_factor, v[[1]])
    result_fit <- rbind(result_fit, c(v[[2]],  v[[3]],  v[[4]], v[[5]]))
}
}

# creation fit result (linear coeff) table
colnames(result_fit) <- c('DATA', 'beta', 'sd_beta', 'info')
result_fit$beta <- as.numeric(result_fit$beta)
result_fit$sd_beta <- as.numeric(result_fit$sd_beta)

# creation Bayes factors table
Bayes_factor <- matrix(data = Bayes_factor, nrow = 3, ncol = 5)
colnames(Bayes_factor) <- c('Auronzo', 'Auronzo_2m', 'Castelfranco', 'PortoTolle_2m', 'Roverchiara_2m')
rownames(Bayes_factor) <- c('min','max','ave')
Bayes_factor <- as.table(Bayes_factor)

write.table(Bayes_factor, file=paste(path, 'PlaceAnalysis/Bayes_factor.ods', sep=''))

                                             # compare data from different place

# This function returns a pdf that contains 3 plot related to Bayes comparison between the two places
#  - min temperature
#  - ave temperature
#  - max temperature
Linear_Trend <- function(Dataset, place1, place2) {

  # plot
  pdf(paste(path, 'PlaceComparison/', comment(place1), '-', comment(place2), '.pdf', sep=''))
  par(mfrow = c(3,1))
  
  for (Temp in c('ave', 'min', 'max')) {
    
    posto_1 <- paste(comment(place1), Temp, sep='_')
    media_1 <- filter(Dataset, DATA==posto_1)$beta
    std_1 <- filter(Dataset, DATA==posto_1)$sd_beta
    
    posto_2 <- paste(comment(place2), Temp, sep='_')
    media_2 <- filter(Dataset, DATA==posto_2)$beta
    std_2 <- filter(Dataset, DATA==posto_2)$sd_beta
    
    # parameters for the difference in normal approximation
    mean <- media_1 - media_2
    sd <- sqrt(std_1**2 + std_2**2)
    lower_limit <- mean - 1.96*sd
    upper_limit <- mean + 1.96*sd
    
    x <- seq(mean-5*sd, mean+5*sd, length=1000)
    diff <- dnorm(x, mean, sd)
    
    plot(x, diff, type='l', xlab='Diff linear coeff.', ylab='Prob',
         main=paste(comment(place1), '-', comment(place2), '    raising coeff', Temp))
    legend("topright", 
           legend = c(paste('mean =  ', round(mean, 3)),
                      paste('sd =  ', round(sd, 3)),
                      paste('95% limits : ', '[',round(lower_limit,3) , ',', round(upper_limit,3) , ']')),
           lty = c(0, 0), col = c('', ''), lwd = 2, title='Results', cex=0.8)
    minor.tick(nx=4, ny=5, tick.ratio=0.5)
    grid()
    polygon (c(mean-2*sd, x[x>mean-2*sd & x<mean+2*sd], mean+2*sd),
             c(0, diff[x>mean-2*sd & x<mean+2*sd], 0), col=rgb(1, 0, 0,0.5))
    abline(v = 0, col="black", lwd=2, lty=2)
  }
  
  o<- dev.off()
}


Linear_Trend(result_fit, Auronzo, Castelfranco)

Linear_Trend(result_fit, Auronzo_2m, PortoTolle_2m)
Linear_Trend(result_fit, Auronzo_2m, Roverchiara_2m)
Linear_Trend(result_fit, PortoTolle_2m, Roverchiara_2m)


                                                 # comparison with paper results

comparison <- function(Dataset, Temp, PaperResult, place, city) {
  
  if (city==FALSE) {
    m1 <- filter(Dataset, DATA==paste('Auronzo_2m_', Temp, sep=''))$beta
    m2 <- filter(Dataset, DATA==paste('PortoTolle_2m_', Temp, sep=''))$beta
    m3 <- filter(Dataset, DATA==paste('Roverchiara_2m_', Temp, sep=''))$beta
    media <- c(m1, m2, m3)
    sd1 <- filter(Dataset, DATA==paste('Auronzo_2m_', Temp, sep=''))$sd_beta
    sd2 <- filter(Dataset, DATA==paste('PortoTolle_2m_', Temp, sep=''))$sd_beta
    sd3 <-filter(Dataset, DATA==paste('Roverchiara_2m_', Temp, sep=''))$sd_beta
    std <- c(sd1, sd2, sd3)
    media_dati <- sum(media/(std**2)) / sum(1/(std**2))
    std_dati <- 1/sqrt(sum(1/(std**2)))
  } else if (city==TRUE) {
    riga <- paste(comment(place), Temp, sep='_')
    media_dati <- filter(Dataset, DATA==riga)$beta
    std_dati <- filter(Dataset, DATA==riga)$sd_beta
  }
  
  # data converted in 10 years (variable transformation)
  media_10y <- media_dati * 10
  std_10y <- std_dati * 10
  
  # paper result
  media_paper <- PaperResult[1]
  std_paper <- PaperResult[2]
  
  # new distribution in norm approximation
  check_media <- media_paper - media_10y
  check_std <- sqrt(std_paper**2 + std_10y**2)
  lower_limit <- check_media - 1.96*check_std
  upper_limit <- check_media + 1.96*check_std
  
  # x variables for the next three plots
  x1 <- seq(media_10y-3.5*std_10y, media_10y+3.5*std_10y, length=1000)
  x2 <- seq(media_paper-3.5*std_paper, media_paper+3.5*std_paper, length=1000)
  x3 <- seq(check_media-3.5*check_std, check_media+3.5*check_std, length=1000)
  
  if (city==FALSE) {
    title <- paste('Comparison', '_', Temp, '.pdf', sep='')
  } else if (city==TRUE) {
    title <- paste('Comparison', '_', comment(place), '_', Temp, '.pdf', sep='')
  }
  
  grDevices::cairo_pdf(paste(path, 'PaperComparison/', title, sep=''))
  area <- matrix( c( 1, 1, 0, 2, 2, 0, 3, 3, 3, 0), nrow=2, byrow=TRUE)
  layout(area)
  
  # plot Arpav result
  plot(x1, dnorm(x1, media_10y, std_10y), type='l', xlab='ΔT/10ys [°C]', ylab='prob', main='Arpav Data')
  legend("topleft",
         legend = c(paste('μ :  ', round(media_10y, 3)),
                    paste('σ :  ', round(std_10y, 3))),
         lty = c(0, 0), col = c('', ''), lwd = 2, title='Gaussian', cex=0.7, xpd=TRUE)
  minor.tick(nx=2, ny=4, tick.ratio=0.5)
  
  # plot paper result
  plot(x2, dnorm(x2, media_paper, std_paper), type='l', xlab='ΔT/10ys [°C]', ylab='prob', main='Paper Data')
  legend("topleft",
         legend = c(paste('μ :  ', round(media_paper, 3)),
                    paste('σ :  ', round(std_paper, 3))),
         lty = c(0, 0), col = c('', ''), lwd = 2, title='Gaussian', cex=0.7, xpd=TRUE)
  minor.tick(nx=2, ny=4, tick.ratio=0.5)
  
  # plot comparison result
  plot(x3, dnorm(x3, check_media, check_std), type='l', xlab='diff Papar - Arpav', ylab='prob', main='Results comparison')
  legend("topleft",
         legend = c(paste('μ :  ', round(check_media, 3)),
                    paste('σ :  ', round(check_std, 3)),
                    paste('95% limits : ', '[',round(lower_limit,3) , ',', round(upper_limit,3) , ']')),
         lty = c(0, 0, 0), col = c('', '', ''), lwd = 2, title='Gaussian', cex=0.7, xpd=TRUE)
  polygon (c(lower_limit, x3[x3>lower_limit & x3<upper_limit], upper_limit),
           c(0, dnorm(x3, check_media, check_std)[x3>lower_limit & x3<upper_limit], 0), col=rgb(1, 0, 0,0.5))
  abline(v=0, col="black", lwd=2, lty=2)
  minor.tick(nx=4, ny=4, tick.ratio=0.5)
  
  o<- dev.off()
  
}

A <- 'ave'
M <- 'max'
m <- 'min'

comparison(result_fit, A, c(0.38, 0.05), Auronzo, FALSE)
comparison(result_fit, M, c(0.42, 0.06), Auronzo, FALSE)
comparison(result_fit, m, c(0.34, 0.04), Auronzo, FALSE)

for (town in list(Auronzo_2m, PortoTolle_2m, Roverchiara_2m)) {
  comparison(result_fit, A, c(0.38, 0.05), town, TRUE)
  comparison(result_fit, M, c(0.42, 0.06), town, TRUE)
  comparison(result_fit, m, c(0.34, 0.04), town, TRUE)
}


                                                  # future prediction with Arima 
Arima_fit <- function(place, n) {
  
  for (Temp in c('ave', 'min', 'max')) {
    
    # prepare data to analyse
    v <- place[place$info==Temp, ]$Medio
    tmp <- data.frame(t=v)
    tmp$n <- (1:nrow(tmp))%/%n
    
    # compute the average over n years
    tmp2 <- tmp |> group_by(n) |> summarise(Medio=mean(t))
    v <- tmp2$Medio
    
    
    pdf(paste(path, 'Arima/ArimaFit_', comment(place), '_', Temp, '.pdf', sep=''))
    par(mfrow = c(2,2))
    
    # create time series and plot it
    tempo <- ts(v, start=min(place$Anno), frequency=1/n)
    plot(tempo, type='p', main=paste(comment(place), Temp, sep='_'),  xlab='Years', ylab='Temperature °C')
    lines(tempo)
    minor.tick(nx=4, ny=4, tick.ratio=0.5)
    
    # plot acf function of the time series
    acf(tempo, lag.max=length(v)%/%3)
    
    # create Arima model
    Arima_model <- auto.arima(tempo)
    
    # plot Arima fit result over data
    plot(tempo, type='l', main='Fit with ARIMA model', xlab='Years', ylab='Temperature °C')
    lines(fitted(Arima_model), col='blue')
    minor.tick(nx=4, ny=4, tick.ratio=0.5)

    # use Arima model to forecast the temperature
    future <- forecast(Arima_model, h=max(1, 5%/%n))
    
    plot(future, xlab='years', ylab='Temperature °C')
    grid()
    minor.tick(nx=4, ny=4, tick.ratio=0.5)
    
    o <- dev.off()
  }
}

for (zone in list(Auronzo, Castelfranco, Auronzo_2m, PortoTolle_2m, Roverchiara_2m)) {
  Arima_fit(zone, 4)
}
