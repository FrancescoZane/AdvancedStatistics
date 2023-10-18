library(Hmisc)
library(tibble)
library(coda)
library(rjags)
library("Rlab")   


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


cat('.bug file used to create the model:  ', '\n', '\n')
cat('model {', '\n', '# data likelihood' ,'\n', 'X ~ dbinom(p, Tot);','\n' ,'\n')
cat(' # a uniform prior for lambda', '\n', 'p ~ dbeta(1, 10);' ,'\n', '}' ,'\n','\n')
cat('----------------------------------------------------------', '\n')
cat('The plots about the reslts of the mcmc are inside: Es4_chain.pdf', '\n')
cat('The inference on p is shown in: Es4_Gibbs_hist.pdf', '\n', '\n')
cat('Informations about the mcmc and the resulting posterior : ', '\n')
print(summary(chain))

cat('----------------------------------------------------------', '\n')	



