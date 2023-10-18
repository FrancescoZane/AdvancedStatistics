library(tibble)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(plot.matrix)
set.seed(17052001)

cat(' #########  ESERCIZIO 1   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('1.1: See Es_1_unif_prior.pdf', '\n', '\n')

percentile <- function(k, input, output) {			
		# k = value of the cumulative distribution
		# output = probability distribution
		# input = relative cumulative distribution
		index <- which.min(abs(input - k))	
		return(output[index])
		}
		
post_function <- function(prior, likelihood) {
		posterior <- likelihood * prior
		k <- 50/10**4*sum(posterior)
		posterior <- posterior/k
		return(posterior)
		}

# data
N = c(5, 8, 4, 6, 11, 6, 6, 5, 6, 4)

# compute interval possible value of mu
mu <- seq(0.001, 50, length.out = 10**4)

# compute the uniform prior
prior_unif <- rep(1, 10**4)

# compute the likelihood
likelihood <- 1
for (num in N) { likelihood <- likelihood * dpois(num, mu) }

# compute the posterior
post_unif <- post_function(prior_unif, likelihood)
		   	
# compute summarized results	   
mean_unif <- 50/10**4 * sum(mu * post_unif)
var_unif <- 50/10**4 * sum((mu - mean_unif)**2 * post_unif)
median_unif <- percentile(0.5, cumsum(post_unif/10**4*50), mu)
low_limit_unif <- percentile(0.025, cumsum(post_unif/10**4*50), mu)
upper_limit_unif <- percentile(0.975, cumsum(post_unif/10**4*50), mu)

# set terminal output
cat('With uniform prior we obtain: ', '\n',
    ' Mean:        ', mean_unif, '\n',
    ' Median:      ', median_unif, '\n',
    ' Var:         ', var_unif, '\n',
    ' Low limit:   ', low_limit_unif, '\n',
    ' Upper limit: ', upper_limit_unif, '\n')

pdf("Es_1_unif_prior.pdf")
plot(mu, post_unif, type='l', main="Es 1 : POSTERIOR with UNIFORM PRIOR", 
     xlab="mu", ylab="Probability     P(mu | D)", col=1, xlim=c(3, 10))
abline(v=mean_unif, col=2, lty=1)
abline(v=median_unif, col=4, lty=1)
abline(v=c(low_limit_unif, upper_limit_unif), col=3, lty=2)
legend(x = "topright", legend = c("Posterior", "Mean", "Median", "95% limit"), lty = c(1, 1, 1, 2),           
       col = c(1, 2, 4, 3), lwd = 2)                
minor.tick(nx=5, ny=2, tick.ratio=0.5)
out<-dev.off()

	
cat('----------------------------------------------------------', '\n')	
cat('1.2: See Es_1_Jeff_prior.pdf', '\n', '\n')

# compute Jeff prior
prior_jeff <- 1/sqrt(mu)
	
# compute the posterior
post_jeff <- post_function(prior_jeff, likelihood)
	
# compute summarized results	   	   
mean_jeff <- 50/10**4 * sum(mu * post_jeff)
var_jeff <- 50/10**4 * sum((mu - mean_unif)**2 * post_jeff)
median_jeff <- percentile(0.5, cumsum(post_jeff/10**4*50), mu)
low_limit_jeff <- percentile(0.025, cumsum(post_jeff/10**4*50), mu)
upper_limit_jeff <- percentile(0.975, cumsum(post_jeff/10**4*50), mu)

# set terminal output
cat('With Jeffrey prior we obtain: ', '\n',
    ' Mean:        ', mean_jeff, '\n',
    ' Median:      ', median_jeff, '\n',
    ' Var:         ', var_jeff, '\n',
    ' Low limit:   ', low_limit_jeff, '\n',
    ' Upper limit: ', upper_limit_jeff, '\n')	

pdf("Es_1_Jeff_prior.pdf")
plot(mu, post_jeff, type='l', main="Es 1 : POSTERIOR with JEFFREYS' PRIOR", 
     xlab="mu", ylab="Probability     P(mu | D)", col=1, xlim=c(3, 10))
abline(v=mean_jeff, col=2, lty=1)
abline(v=median_jeff, col=4, lty=1)
abline(v=c(low_limit_jeff, upper_limit_jeff), col=3, lty=2)
legend(x = "topright", legend = c("Posterior", "Mean", "Median", "95% limit"), lty = c(1, 1, 1, 2),           
       col = c(1, 2, 4, 3), lwd = 2)                
minor.tick(nx=5, ny=2, tick.ratio=0.5)
out<-dev.off()

cat('----------------------------------------------------------', '\n')	
cat('1.3: See the following table', '\n', '\n')

# UNIFORM NORMAL APPROXIAMTION

# compute my norm approximation using the calculated mean and sd
unif_norm <- dnorm(mu, mean_unif, sqrt(var_unif))

mean_unif_approx <- 50/10**4 * sum(mu * unif_norm)
var_unif_approx <- 50/10**4 * sum((mu - mean_unif_approx)**2 * unif_norm)
median_unif_approx <- percentile(0.5, cumsum(unif_norm/10**4*50), mu)
low_limit_unif_approx <- percentile(0.025, cumsum(unif_norm/10**4*50), mu)
upper_limit_unif_approx <- percentile(0.975, cumsum(unif_norm/10**4*50), mu)

# JEFFREY NORMAL APPROXIAMTION

# compute my norm approximation using the calculated mean and sd
jeff_norm <- dnorm(mu, mean_jeff, sqrt(var_jeff))

mean_jeff_approx <- 50/10**4 * sum(mu * jeff_norm)
var_jeff_approx <- 50/10**4 * sum((mu - mean_jeff_approx)**2 * jeff_norm)
median_jeff_approx <- percentile(0.5, cumsum(jeff_norm/10**4*50), mu)
low_limit_jeff_approx <- percentile(0.025, cumsum(jeff_norm/10**4*50), mu)
upper_limit_jeff_approx <- percentile(0.975, cumsum(jeff_norm/10**4*50), mu)

# put all results in a summary table
summary1 <- tibble(distribution=c('Uniform','Jeffrey','Norm_unif','Norm_Jeff'),
		   mean=c(mean_unif, mean_jeff, mean_unif_approx, mean_jeff_approx),
		   median=c(median_unif, median_jeff, median_unif_approx, median_jeff_approx),
		   var=c(var_unif, var_jeff, var_unif_approx, var_jeff_approx),
		   low_limit=c(low_limit_unif, low_limit_jeff, low_limit_unif_approx, low_limit_jeff_approx),
		   upper_limit=c(upper_limit_unif, upper_limit_jeff, upper_limit_unif_approx, upper_limit_jeff_approx))

print(summary1)
cat('----------------------------------------------------------', '\n', '\n')	


cat(' #########  ESERCIZIO 2   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('2.1: ', '\n', '\n')

prob_fail <- 0.15

cat('The probability distribution of y, the number of times the new method', '\n',
    'fails to detect the disease, is a binomial distribution with:', '\n') 
cat(' p = the prob of fail of the new test', '\n',
    'y = the number of fails', '\n',
    'n = the number of total trials', '\n')

cat('----------------------------------------------------------', '\n')	
cat('2.2: ', '\n', '\n')	

n <- 75
y <- 6 
p_f <- y/n

cat('The frequentist estimator is         p_f = y/n' , '\n')
p_f <- y/n
cat('It is un unbiased estimator that follows a binomial distribution with:', '\n')
cat('Mean = p_f = ', p_f, '\n')
cat('Variance = p_f * (1-p_f) = ', p_f * (1-p_f), '\n')

cat('----------------------------------------------------------', '\n')	
cat('2.3: See Es_2_beta_prior.pdf', '\n')	

# value for my beta
mean <- 0.15
std <- 0.14 

# parameter of my beta to obtain the above values
a <- mean * ((mean*(1-mean))/(std**2)-1)
b <- (1-mean)/(mean)*a

# compute interval possible value of p
prob <- seq(from=1/10**4, to=1, length=10**4)

# compute beta prior
prior <- dbeta(prob, a, b)

# compute the likelihood
likelihood <- dbinom(y, n, prob)

# compute the posterior and normalize it
posterior <- prior*likelihood
Z <- sum(posterior/10**4)
posterior <- posterior/Z

post_mean <- sum(prob*posterior/10**4)
post_var <- sum((prob-post_mean)**2*posterior/10**4)


pdf("Es_2_beta_prior.pdf")
plot(prob, posterior, type='l', main="Es 2 : POSTERIOR FOR THE NEW METHOD", 
     xlab="p", ylab="Probability     P(p | n,y)", col=1, lwd = 2)
abline(v=c(post_mean-sqrt(post_var), post_mean+sqrt(post_var)), col=2, lty=2, lwd = 2)
abline(v=post_mean, col=5, lty=1, lwd = 2)
legend(x = "topright", legend = c("Posterior", "Mean +- sigma", "Mean"), lty = c(1, 2, 1),           
       col = c(1, 2, 5), lwd = 2)                
out<-dev.off()



cat('----------------------------------------------------------', '\n')	
cat('2.4: ', '\n','My NULL hypothesis is       H0 : p > 0.15  (the new test is worse)', '\n', '\n')	

percentile <- function(k, input, output) {			
		# k = value of the cumulative distribution
		# output = probability distribution
		# input = relative cumulative distribution
		index <- which.min(abs(input - k))	
		return(output[index])
		}

# using the computed posterior calculate the prob that p is bigger than 0.15
# computing the integral with discretize method
prob_b_bigger_0.15 <- sum(posterior[prob>0.15]/10**4)

cat('P(H0 | data) = P(y>0.15 | data) =', prob_b_bigger_0.15, '\n')
cat('Since this probability is smaller than 0.05 we reject the null hypothesis.', '\n')

cat('We can visualize this statement also graphically: see Es_2_Bayesian.pdf', '\n')
limit <- percentile(0.95, cumsum(posterior/10**4), prob)

pdf("Es_2_Bayesian.pdf")
plot (prob, posterior, type='l', main="Es 2 : BAYESIAN HT", 
     xlab="p", ylab="Probability     P(p | n,y)")
polygon ( c(limit,prob[prob>limit]) , c(0,posterior[prob>limit]), col=rgb(1, 0, 0,0.5))
abline(v = 0.15, col='blue', lty=2)
legend(x = "topright", legend = c("Posterior", "5% limit area", "p = 0.15"), lty = c(1, 1, 2),           
       col = c(1, 'red', 'blue'), lwd = 2)                
out<-dev.off()

cat('----------------------------------------------------------', '\n')	
cat('2.5: ', '\n','My NULL hypothesis is       H0 : p > 0.15  (the new test is worse)', '\n', '\n')

# compute the likelihood in considering H0 true
true_null <- dbinom(x = c(0:n), size = n, prob = 0.15)

# put the obtained value in a table just for number between 0 and 10 (interested region)
freq <- tibble(
	number=c(0:10) ,  
	prob=true_null[1:11]
	)

# compute the p-value by discretized method
p_value<-tail(cumsum(freq$prob[freq$number<=6]),1)
	
cat('p-value =', p_value, '\n')
cat('Since the p-value is bigger than 0.05 we can not reject the null hypothesis.', '\n')

cat('We can visualize this statement also graphically: see Es_2_Frequentist.pdf', '\n')

freq['cumulative'] <- cumsum(freq$prob)

pdf("Es_2_Frequentist.pdf")
ggplot(freq, aes(x=number, y=cumulative))+ 
geom_bar(stat = "identity")+
geom_hline(yintercept=0.05, color = "red")+
scale_x_continuous(breaks = seq(0,10,2))+
scale_y_continuous(breaks = seq(0,1,0.05))+
labs(x='Number fail detections', y='F(y)', title='Es 2 : FREQUENTIST HT', colour = "Legend")              
out<-dev.off()

cat('----------------------------------------------------------', '\n', '\n')	


cat(' #########  ESERCIZIO 3   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	

# Function to simulate uor data
simulation <- function(n, a, b) {
	      
	      tetha <- pi*(runif(n)-0.5)
	      x_k <- b * tan(tetha) + a
	      return(x_k)
	      }

# Choose the prior 
unif_prior <- function(a, b) {return(1)}  

# Function that creates the prior 
prior_function <- function(alpha, beta, f) {
		matrice <- matrix(data=NA, nrow=length(beta), ncol=length(alpha))
		
		for (a in 1:length(alpha)) {
			for (b in 1:length(beta)) {
				matrice[b,a]<-f(a,b)
				}	
			}
		
		return(matrice)
		} 
		
# Function that compute the likelihood for one measure
# k is my measure
single_likelihood_function <- function(k, alpha, beta, delta_alpha, delta_beta) {
				like <- matrix(data=NA, nrow=length(beta), ncol=length(alpha))
				
				for (a in 1:length(alpha)) {
					for (b in 1:length(beta)) {
						like[b,a] <- 1/pi * (beta[b])/(beta[b]**2 + (k-alpha[a])**2)
						}	
					}
					
				Z<-sum(like*delta_alpha*delta_beta)
				like<-like/Z
				
				return(like)	
				}
			
# Function that compute the likelihood for all measures combinig the likelihood for the signle measure				
total_likelihood_function <- function(measures, alpha, beta, delta_alpha, delta_beta, f) {
				like <- matrix(data=1, nrow=length(beta), ncol=length(alpha))
				
				for (k in measures) {
					like <- like * f(k, alpha, beta, delta_alpha, delta_beta)
					Z<-sum(like*delta_alpha*delta_beta)
					like<-like/Z
					}

				return(like)
				}

# True values
alpha_true <- 0
beta_true <- 1.5

# create simulated data
measures <- simulation(100, alpha_true, beta_true)

# compute interval possible value of alpha and beta
alpha <- seq(alpha_true-2, alpha_true+2, length.out = 100)
delta_alpha <- (max(alpha)-min(alpha))/(length(alpha)-1)
beta <- seq(0, beta_true*2, length.out = 100)
delta_beta <- (max(beta)-min(beta))/(length(beta)-1)

# compute my 2d prior
prior <- prior_function(alpha, beta, unif_prior)

# compute the total likelihood
likelihood <- total_likelihood_function(measures, alpha, beta, delta_alpha, delta_beta, single_likelihood_function)

# compute my 2d posterior and rescale it in order to obain values between 0 and 1
posterior <- prior * likelihood
posterior <- posterior/max(posterior)

pdf("Es_3_alpha-beta_posterior.pdf")
contour(alpha, beta, posterior, nlevels=5, main='Es 3 : POSTERIOR FOR ALPHA AND BETA', xlab='alpha [Km]', ylab='beta [Km]')
abline(v=alpha_true, h=beta_true, col='grey')              
out<-dev.off()


# Marginalization to obtain alpha posterior and normalize it
alpha_posterior <- apply(posterior, 1, sum)
k <- sum(alpha_posterior*delta_alpha) 
alpha_posterior <- alpha_posterior / k

percentile <- function(k, input, output) {			
		# k = value of the cumulative distribution
		# output = probability distribution
		# input = relative cumulative distribution
		index <- which.min(abs(input - k))	
		return(output[index])

		}
alpha_MPV <- alpha[alpha_posterior==max(alpha_posterior)]		
alpha_low_limit <- percentile(0.025, cumsum(alpha_posterior*delta_alpha), alpha)
alpha_upper_limit <- percentile(0.975, cumsum(alpha_posterior*delta_alpha), alpha)

pdf("Es_3_alpha.pdf")
plot(alpha, alpha_posterior, type='l', main='Es 3 : POSTERIOR FOR ALPHA', xlab='alpha [Km]', ylab='p(alpha | data)')
polygon (c(alpha_low_limit, alpha[alpha>alpha_low_limit & alpha<alpha_upper_limit], alpha_upper_limit),
         c(0, alpha_posterior[alpha>alpha_low_limit & alpha<alpha_upper_limit], 0), col=rgb(1, 0, 0,0.5))
abline(v = alpha_true, col='green', lty=2)
abline(v = alpha_MPV, col='blue', lty=2)
legend(x = "topright", legend = c("Posterior for alpha", '95% cred. int.', "True alpha", "Alpha prediction"),
       lty = c(1, 1, 2, 2), col = c(1, 'red','green', 'blue'), lwd = 2) 
minor.tick(nx=5, ny=4, tick.ratio=0.5)               
out<-dev.off()


# Marginalization to obtain beta posterior
beta_posterior <- apply(posterior, 2, sum)
k <- sum(beta_posterior*delta_beta) 
beta_posterior <- beta_posterior / k

beta_MPV <- beta[beta_posterior==max(beta_posterior)]
beta_low_limit <- percentile(0.025, cumsum(beta_posterior*delta_beta), beta)
beta_upper_limit <- percentile(0.975, cumsum(beta_posterior*delta_beta), beta)

pdf("Es_3_beta.pdf")
plot(beta, beta_posterior, type='l', main='Es 3 : POSTERIOR FOR BETA', xlab='beta [Km]', ylab='p(beta | data)')
polygon (c(beta_low_limit, beta[beta>beta_low_limit & beta<beta_upper_limit], beta_upper_limit),
         c(0, beta_posterior[beta>beta_low_limit & beta<beta_upper_limit], 0), col=rgb(1, 0, 0,0.5))
abline(v = beta_true, col='green', lty=2)
abline(v = beta_MPV, col='blue', lty=2)
legend(x = "topright", legend = c("Posterior for beta", '95% cred. int.', "True beta", "Beta prediction"),
       lty = c(1, 1, 2, 2), col = c(1, 'red','green', 'blue'), lwd = 2) 
minor.tick(nx=5, ny=5, tick.ratio=0.8)               
out<-dev.off()

# put the results in a summarized table
summary3 <- tibble(Variable=c('Alpha', 'Beta'),
		   Estimation_MPV=c(alpha_MPV, beta_MPV),
		   Low_limit=c(alpha_low_limit, beta_low_limit),
		   Upper_limit=c(alpha_upper_limit, beta_upper_limit))


cat('Real value of alpha: ', alpha_true, '\n')
cat('               beta: ', beta_true, '\n', '\n')

cat('What we obtain with inference is: ', '\n', '\n')
print(summary3)

cat('\n')
cat('To visualize these results see the following plots:', '\n',
    ' - for the 2d posterior     ->   Es_3_alpha-beta_posterior.pdf', '\n',
    ' - for the alpha posterior  ->   Es_3_alpha.pdf', '\n',
    ' - for the beta posterior   ->   Es_3_beta.pdf', '\n')

cat('----------------------------------------------------------', '\n', '\n')	


cat(' #########  ESERCIZIO 4   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	

# this function, given the true value of the parameter, returns to us:
# - mean=c(mean_a, mean_b) -> the predicted values for A and B
# - sd=c(sd_a, sd_b) -> the predicted values for sigma(A) and sigma(B) (of their posteriors)
# - cov=cov_ab -> the comvariance between A and B
# - rho=rho_ab -> the rho value
# - post=totalpost -> the 2d plot of the posterior
# - A_post=Apost -> the plot of the marginalized A posterior
# - B_post=Bpost -> the plot of the marginalized B posterior

analysis <- function(x_0, w, A, B, delta_t, segnale, v) {
	
	x <- seq(from=-7, to=7, by=0.5)
	S <- segnale(x , A, B, x_0, w, delta_t)          # biuld the real signal
	ddat <- rpois(length(S), S)	 		 # each point is a Possonian dist with mean S
 
	simulation_data <- tibble(x=x, signal=S, data=ddat)

	alim <- c(0, 4)
	blim <- c(0, 2)
	Nsamp <- 100
	uniGrid <- seq(from=1/(2*Nsamp),
	to=1-1/(2*Nsamp), by=1/Nsamp)
	delta_a <- diff(alim )/ Nsamp
	delta_b <- diff(blim )/ Nsamp
	a <- alim[1] + diff(alim )* uniGrid    # create sequence of value for a
	b <- blim[1] + diff(blim )* uniGrid    # create sequence of value for b
	
	
	# Log posterior   (at fixed A, B and varing the data)
	log.post <- function(d, x, A, B, x_0, w, t) {
		if(A<0 || B<0) {return(-Inf )} # the effect of the prior
		sum(dpois(d, lambda=segnale(x, A, B, x_0, w, t), log=TRUE)) # somma di prob di avere d come
									    # risultato sapendo che segue una
									    # dist poissoniana con media lambda
		}
	
	z <- matrix(data=NA , nrow=length(a), ncol=length(b))
	
	# creating my log posterior
	for(j in 1:length(a)) {
		for(k in 1:length(b)) {
			z[j,k] <- log.post(ddat , x , a[j], b[k], x_0, w, delta_t)
			}
		}
	
	z <- z - max(z) # set maximum to zero
	
	# Plot unnormalized 2D posterior as contours
	points <- expand.grid(a, b)
        posterior <- tibble(a=points[,1], b=points[,2], post=as.vector(exp(z)))
	
	if (v==1) {h_text <- paste( 'Value of w: ' , w, sep=' ')}
	else {h_text <- paste( 'The ratio A/B: ' , A/B, sep=' ')}
	
	# creating the plot for 2d posterior
        totalpost <- ggplot(posterior, aes(a, b, z=post))+
                  geom_contour_filled()+
                  theme(legend.key.size = unit(0.2, 'cm'), legend.title=element_text(size=5), legend.text=element_text(size=5))+
                  #ylim(0,2)+
                  #xlim(0,4)+
                  geom_hline(yintercept=B, color='red', linetype=2)+
                  geom_vline(xintercept=A, color='red', linetype=2)+
                  labs(x='Amplitude, A', y='Background, B', title=h_text)
	
	
	# Computing the marginalied A posterior
	p_a_D <- apply(exp(z), 1, sum)         # the sum function is applied on the rows of exp(z)
	p_a_D <- p_a_D/( delta_a*sum(p_a_D))   # normalize the posterior
	
	# Computing the marginalied B posterior
	p_b_D <- apply(exp(z), 2, sum)         # the sum function is applied on the columns of exp(z)
	p_b_D <- p_b_D/( delta_b*sum(p_b_D))   # normalize the posterior
	

	# creating plots for A and B marginalized posterior
	A_post <- tibble(A=a, post=p_a_D)
	Apost <- ggplot(data=A_post, mapping=aes(x=A, y=post)) +
	geom_line()+
	labs(x='Amplitude, A', y='P(A | data)', title=h_text)+
	geom_vline(xintercept=A, color='red')
	
	B_post <- tibble(B=b, post=p_b_D)
	Bpost <- ggplot(data=B_post, mapping=aes(x=B, y=post)) +
	geom_line()+
	labs(x='Background, B', y='P(B | data)', title=h_text)+
	geom_vline(xintercept=B, color='red')
	
	
	# Compute mean, standard deviation, covariance, correlation of/between A and B
	mean_a <- delta_a * sum(a * p_a_D)
	mean_b <- delta_b * sum(b * p_b_D)
	
	sd_a <- sqrt( delta_a * sum((a-mean_a)**2 * p_a_D) )
	sd_b <- sqrt( delta_b * sum((b-mean_b)**2 * p_b_D) )
	
	cov_ab <- 0
	for(j in 1:length(a)) {
		for(k in 1:length(b)) {
			cov_ab <- cov_ab + (a[j]-mean_a)*(b[k]-mean_b)*exp(z[j,k])
			}
		}
		
	cov_ab <- cov_ab / sum(exp(z))
	rho_ab <- cov_ab / (sd_a * sd_b)
	
	my_list <- list('mean'=c(mean_a, mean_b),
			'sd'=c(sd_a, sd_b),
			'cov'=cov_ab,
			'rho'=rho_ab,
			'post'=totalpost,
			'A_post'=Apost,
			'B_post'=Bpost)	
	return(my_list)
	}

segnale <- function(x, a, b, x0, w, t) {
		S <- t * (a*exp(-(x-x0)**2/(2*w**2)) + b)
		return(S)
		}


# Parametri modello 
x_0 <- 0
w <- 1
A <- 2
B <- 1
delta_t <- 10


cat('                      Varing the parameter w', '\n', '\n')

cat(' - 2d posterior varing w:                 See Es_4_2d_varing_w.pdf', '\n')
cat(' - A marginalized posterior varing w:     See Es_4_Apost_varing_w.pdf', '\n')
cat(' - B marginalized posterior varing w:     See Es_4_Bpost_varing_w.pdf', '\n', '\n')

# here I'll put all the objects
vettore <- c()

for (w in c(0.1, 0.25, 1, 2, 3 )) {
	
	set.seed(17052001)
	name <- paste('tabella_', w, sep='')
	tabella <- c(analysis(x_0, w, A, B, delta_t, segnale,1))
	vettore <- c(vettore, tabella)

	}


pdf("Es_4_2d_varing_w.pdf")
grid.arrange(vettore[5]$post, vettore[12]$post, vettore[19]$post, vettore[26]$post, vettore[33]$post, nrow = 3)        
out<-dev.off()


pdf("Es_4_Apost_varing_w.pdf")
grid.arrange(vettore[6]$A_post, vettore[13]$A_post, vettore[20]$A_post, vettore[27]$A_post, vettore[34]$A_post, nrow = 2)             
out<-dev.off()
 

pdf("Es_4_Bpost_varing_w.pdf")
grid.arrange(vettore[7]$B_post, vettore[14]$B_post, vettore[21]$B_post, vettore[28]$B_post, vettore[35]$B_post, nrow = 2)              
out<-dev.off()


meanA <- c(vettore[1]$mean[1], vettore[8]$mean[1], vettore[15]$mean[1], vettore[22]$mean[1], vettore[29]$mean[1])
meanB <- c(vettore[1]$mean[2], vettore[8]$mean[2], vettore[15]$mean[2], vettore[22]$mean[2], vettore[29]$mean[2])
sdA <- c(vettore[2]$sd[1], vettore[9]$sd[1], vettore[16]$sd[1], vettore[23]$sd[1], vettore[30]$sd[1])
sdB <- c(vettore[2]$sd[2], vettore[9]$sd[2], vettore[16]$sd[2], vettore[23]$sd[2], vettore[30]$sd[2])
cov <- c(vettore[3]$cov, vettore[10]$cov, vettore[17]$cov, vettore[24]$cov, vettore[31]$cov)
rho <- c(vettore[4]$rho, vettore[11]$rho, vettore[18]$rho, vettore[25]$rho, vettore[32]$rho)

varing_w <- tibble(w=c(0.1, 0.25, 1, 2, 3 ), mean_A=meanA, mean_B=meanB, SdA=sdA, SdB=sdB, COV=cov, rho=rho)

cat('Summarized results:', '\n')
print(varing_w)

cat('\n')
cat('To visualize these data:     See Es_4_2d_varing_w.pdf', '\n')

varing_w['ER_A'] <- sdA/meanA
varing_w['ER_B'] <- sdB/meanB

ma <- ggplot(varing_w, aes(x=w, y=mean_A))+
      geom_point()+
      geom_line()+
      labs(x='w', y='Mean parameter A')+
      geom_hline(yintercept=A, color='red')

mb <- ggplot(varing_w, aes(x=w, y=mean_B))+
      geom_point()+
      geom_line()+
      labs(x='w', y='Mean parameter B')+
      geom_hline(yintercept=B, color='red')

sta <- ggplot(varing_w, aes(x=w))+
       geom_point(aes(y=ER_A))+
       geom_line(aes(y=ER_A))+
       labs(x='w', y='Sigma(A) / A')

stb <- ggplot(varing_w, aes(x=w))+
       geom_point(aes(y=ER_B))+
       geom_line(aes(y=ER_B))+
       labs(x='w', y='Sigma(B) / B')

pdf("Es_4_results_varing_w.pdf")      
grid.arrange(ma, mb, sta, stb, nrow = 2)             
out<-dev.off() 
       
cat('----------------------------------------------------------', '\n')
cat('                      Varing the ratio A/B', '\n', '\n')

cat(' - 2d posterior varing the ratio A/B:                 See Es_4_2d_varing_B.pdf', '\n')
cat(' - A marginalized posterior varing the ratio A/B:     See Es_4_Apost_varing_B.pdf', '\n')
cat(' - B marginalized posterior varing the ratio A/B:     See Es_4_Bpost_varing_B.pdf', '\n', '\n')

vettore <- c()

for (k in c(1, 2, 3, 5, 10)) {
	
	set.seed(17052001)
	name <- paste('tabella_', w, sep='')
	tabella <- c(analysis(x_0, w, A, A/k, delta_t, segnale,2))
	vettore <- c(vettore, tabella)

	}

pdf("Es_4_2d_varing_B.pdf")
grid.arrange(vettore[5]$post, vettore[12]$post, vettore[19]$post, vettore[26]$post, vettore[33]$post, nrow = 3)             
out<-dev.off()

pdf("Es_4_Apost_varing_B.pdf")
grid.arrange(vettore[6]$A_post, vettore[13]$A_post, vettore[20]$A_post, vettore[27]$A_post, vettore[34]$A_post, nrow = 2) 
out<-dev.off()

pdf("Es_4_Bpost_varing_B.pdf")
grid.arrange(vettore[7]$B_post, vettore[14]$B_post, vettore[21]$B_post, vettore[28]$B_post, vettore[35]$B_post, nrow = 2) 
out<-dev.off()

meanA <- c(vettore[1]$mean[1], vettore[8]$mean[1], vettore[15]$mean[1], vettore[22]$mean[1], vettore[29]$mean[1])
meanB <- c(vettore[1]$mean[2], vettore[8]$mean[2], vettore[15]$mean[2], vettore[22]$mean[2], vettore[29]$mean[2])
sdA <- c(vettore[2]$sd[1], vettore[9]$sd[1], vettore[16]$sd[1], vettore[23]$sd[1], vettore[30]$sd[1])
sdB <- c(vettore[2]$sd[2], vettore[9]$sd[2], vettore[16]$sd[2], vettore[23]$sd[2], vettore[30]$sd[2])
cov <- c(vettore[3]$cov, vettore[10]$cov, vettore[17]$cov, vettore[24]$cov, vettore[31]$cov)
rho <- c(vettore[4]$rho, vettore[11]$rho, vettore[18]$rho, vettore[25]$rho, vettore[32]$rho)

varing_B <- tibble(A_over_B=c(1, 2, 3, 5, 10), mean_A=meanA, mean_B=meanB, SdA=sdA, SdB=sdB, COV=cov, rho=rho)

cat('Summarized results:', '\n')
print(varing_B)

cat('\n')
cat('To visualize these data:     See Es_4_results_varing_B.pdf', '\n')

varing_B['ER_A'] <- sdA/meanA
varing_B['ER_B'] <- sdB/meanB

ma <- ggplot(varing_B, aes(x=A_over_B, y=mean_A))+
      geom_point()+
      geom_line()+
      labs(x='A / B', y='Mean parameter A')+
      geom_hline(yintercept=A, color='red')

mb <- ggplot(varing_B, aes(x=A_over_B, y=mean_B))+
      geom_point()+
      geom_line()+
      labs(x='A / B', y='Mean parameter B')+
      geom_hline(yintercept=B, color='red')

sta <- ggplot(varing_B, aes(x=A_over_B))+
       geom_point(aes(y=ER_A))+
       geom_line(aes(y=ER_A))+
       labs(x='A / B', y='Sigma(A) / A')

stb <- ggplot(varing_B, aes(x=A_over_B))+
       geom_point(aes(y=ER_B))+
       geom_line(aes(y=ER_B))+
       labs(x='A / B', y='Sigma(B) / B')
       
pdf("Es_4_results_varing_B.pdf")      
grid.arrange(ma, mb, sta, stb, nrow = 2)   
out<-dev.off()

cat('----------------------------------------------------------', '\n', '\n')	

