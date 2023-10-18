library(tibble)
library(dplyr)
library(ggplot2)
library(gridExtra)
set.seed(170501)

cat(' #########  ESERCIZIO 1   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('1.1: See Prior.pdf', '\n', '\n')
y = 7
n = 20

stepfunction_prior <- function(x) {							# create g() function
			ifelse(x<=0.2, x,
					ifelse(x<=0.3, 0.2,
							ifelse(x<=0.5, 0.5-x, 0							
							)
						)
				)
			}

# Creating a table containing all discrete values of our pdf 
PRIOR <- tibble(probability=seq(1/10**4, 1-1/10**4, length.out = 10**4))       
PRIOR['unif'] <- dbeta(PRIOR$probability,1,1)
PRIOR['Jeff'] <- dbeta(PRIOR$probability,0.5,0.5)
PRIOR['step'] <- stepfunction_prior(PRIOR$probability)

pdf("Prior.pdf")
ggplot(PRIOR, aes(x=probability))+
geom_line(aes(y=step, color='g() function'))+
geom_line(aes(y=Jeff, color='Jeffrey prior'))+
geom_line(aes(y=unif, color='Uniform'))+
ylim(0, 3)+
labs(x='p', y='Probability', title='PRIOR', colour = "Priors")
out<-dev.off()

PRIOR['likelihood'] <- dbinom(x=y, size=n, prob=PRIOR$probability)
cat('\n')
cat('The warning is due to the fact the I impose ylim(0, 3) and for this reason some', '\n')
cat('points belonging to the g() function are not visible in the plot.', '\n')
cat('----------------------------------------------------------', '\n')	
cat('1.2: See Posterior.pdf', '\n', '\n')

post <- function(prior) {						# function that given the prior gives us the normalized posterior
	posterior <- PRIOR$likelihood * prior
	k <- 1/10**4*sum(posterior)					# calculate the intergral of the posterior
	posterior <- 1/k * posterior					# normalize the posterior
	return(posterior)
	}
	
PRIOR['post.unif'] <- post(PRIOR$unif)
PRIOR['post.Jeff'] <- post(PRIOR$Jeff)
PRIOR['post.step'] <- post(PRIOR$step)

pdf("Posteriors.pdf")
ggplot(PRIOR, aes(x=probability))+
geom_line(aes(y=post.step, color='g() function'))+
geom_line(aes(y=post.Jeff, color='Jeffrey prior'))+
geom_line(aes(y=post.unif, color='Uniform'))+
labs(x='p', y='Probability', title='POSTERIOR', colour = "Posteriors")
out<-dev.off()

mean_post.unif <- 1/10**4*sum(PRIOR$probability*PRIOR$post.unif)
mean_post.Jeff <- 1/10**4*sum(PRIOR$probability*PRIOR$post.Jeff)
mean_post.step <- 1/10**4*sum(PRIOR$probability*PRIOR$post.step)

var_post.unif <- 1/10**4*sum((PRIOR$probability-mean_post.unif)**2*PRIOR$post.unif)
var_post.Jeff <- 1/10**4*sum((PRIOR$probability-mean_post.Jeff)**2*PRIOR$post.Jeff)
var_post.step <- 1/10**4*sum((PRIOR$probability-mean_post.step)**2*PRIOR$post.step)

# Creating a summaring table
summary <- tibble(Associated_Prior=c('Uniforn', 'Jeffrey', 'Step g()'),
		  Mean=c(mean_post.unif, mean_post.Jeff, mean_post.step),
		  Var=c(var_post.unif, var_post.Jeff, var_post.step))
cat('Summary table: ', '\n')
print(summary)

cat('----------------------------------------------------------', '\n')	
cat('1.3: Summary table', '\n', '\n')

# build the cumulative distributions for the three posteriors to find the credibility intervals
PRIOR['cum_post.unif'] <- cumsum(PRIOR$post.unif/10**4)			
PRIOR['cum_post.Jeff'] <- cumsum(PRIOR$post.Jeff/10**4)
PRIOR['cum_post.step'] <- cumsum(PRIOR$post.step/10**4)

# function that allows us to find the probability releated to a given value of the cumulative distribution
percentile <- function(k, input, output) {			
		# k = value of the cumulative distribution
		# output = probability distribution
		# input = relative cumulative distribution
		index <- which.min(abs(input - k))	
		return(output[index])
		}
		
per_2.5_unif <- percentile(0.025, PRIOR$cum_post.unif, PRIOR$probability)
per_2.5_Jeff <- percentile(0.025, PRIOR$cum_post.Jeff, PRIOR$probability)
per_2.5_step <- percentile(0.025, PRIOR$cum_post.step, PRIOR$probability)

per_97.5_unif <- percentile(0.975, PRIOR$cum_post.unif, PRIOR$probability)
per_97.5_Jeff <- percentile(0.975, PRIOR$cum_post.Jeff, PRIOR$probability)
per_97.5_step <- percentile(0.975, PRIOR$cum_post.step, PRIOR$probability)

summary['perc_2.5'] <- c(per_2.5_unif, per_2.5_Jeff, per_2.5_step)
summary['perc_97.5'] <- c(per_97.5_unif, per_97.5_Jeff, per_97.5_step)

print(summary)

cat('----------------------------------------------------------', '\n')	
cat('1.4: See 95_confidence.pdf', '\n')

plt_unif <- ggplot(PRIOR, aes(x=probability))+
	    geom_line(aes(y=post.unif))+
	    geom_vline(xintercept = per_2.5_unif, linetype="dotted", color = "blue")+
	    geom_vline(xintercept = per_97.5_unif, linetype="dotted", color = "blue")+
	    geom_vline(xintercept = mean_post.unif, linetype="dashed", color = "red")+
	    labs(x='p', y='Probability', title='Uniform')

plt_Jef <- ggplot(PRIOR, aes(x=probability))+
	   geom_line(aes(y=post.Jeff))+
	   geom_vline(xintercept = per_2.5_Jeff, linetype="dotted", color = "blue")+
	   geom_vline(xintercept = per_97.5_Jeff, linetype="dotted", color = "blue")+
	   geom_vline(xintercept = mean_post.Jeff, linetype="dashed", color = "red")+
	   labs(x='p', y='Probability', title='Jeffrey')

plt_beta <- ggplot(PRIOR, aes(x=probability))+
	    geom_line(aes(y=post.step))+
	    geom_vline(xintercept = per_2.5_step, linetype="dotted", color = "blue")+
	    geom_vline(xintercept = per_97.5_step, linetype="dotted", color = "blue")+
	    geom_vline(xintercept = mean_post.step, linetype="dashed", color = "red")+
	    labs(x='p', y='Probability', title='Step g()')
	    
pdf("95_confidence.pdf")
grid.arrange(plt_unif, plt_Jef, plt_beta, nrow = 1)
out<-dev.off()
cat('----------------------------------------------------------', '\n', '\n')	
	

cat(' #########  ESERCIZIO 2   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('2.1: See Water.pdf', '\n', '\n')

y = 17
n = 116

# Creating a table containing all discrete values of our pdf 
PRIOR2 <- tibble(probability=seq(1/10**4, 1-1/10**4, length.out = 10**4))
PRIOR2['Unif'] <- dbeta(PRIOR2$probability,1,1)
PRIOR2['Beta'] <- dbeta(PRIOR2$probability,1,4)

PRIOR2['likelihood'] <- dbinom(x=y, size=n, prob=PRIOR2$probability)

post <- function(prior) {						# function that given the prior gives us the normalized posterior
	posterior <- PRIOR2$likelihood * prior
	k <- 1/10**4*sum(posterior)					# calculate the intergral of the posterior
	posterior <- 1/k * posterior					# normalize the posterior
	return(posterior)
	}
	
PRIOR2['post.unif'] <- post(PRIOR2$Unif)
PRIOR2['post.beta'] <- post(PRIOR2$Beta)


prior2 <- ggplot(PRIOR2, aes(x=probability))+
	 geom_line(aes(y=Unif, color='Uniform'))+
	 geom_line(aes(y=Beta, color='Beta(1,4)'))+
	 labs(x='p', y='Probability', title='PRIOR', colour = "Priors")
	 
like2 <- ggplot(PRIOR2, aes(x=probability))+
	 geom_line(aes(y=likelihood, color='Binomial'))+
	 labs(x='p', y='Probability', title='LIKELIHOOD', colour = "lLikelihood")
	 
posterior2 <- ggplot(PRIOR2, aes(x=probability))+
	     geom_line(aes(y=post.unif, color='Uniform'))+
	     geom_line(aes(y=post.beta, color='Beta(1,4)'))+
	     labs(x='p', y='Probability', title='POSTERIOR', colour = "Posteriors")

pdf("Water.pdf")
grid.arrange(prior2, like2, posterior2, nrow = 3)
out<-dev.off()


mean_post.unif <- 1/10**4*sum(PRIOR2$probability*PRIOR2$post.unif)
var_post.unif <- 1/10**4*sum((PRIOR2$probability-mean_post.unif)**2*PRIOR2$post.unif)

mean_post.beta <- 1/10**4*sum(PRIOR2$probability*PRIOR2$post.beta)
var_post.beta <- 1/10**4*sum((PRIOR2$probability-mean_post.beta)**2*PRIOR2$post.beta)

# Creating a summaring table
summary2 <- tibble(Associated_Prior=c('Uniforn', 'Beta(1,4)'),
		   Mean=c(mean_post.unif, mean_post.beta),
		   Var=c(var_post.unif, var_post.beta))

cat('Summary table: ', '\n')
print(summary2)   
		   
cat('----------------------------------------------------------', '\n')	
cat('2.2: For the normal approximations I take the means and the square of the variances', '\n', '\n')	
cat('     See Norm_approximation.pdf', '\n', '\n')	   

PRIOR2['norm_unif'] <- dnorm(PRIOR2$probability, mean_post.unif, sqrt(var_post.unif))
PRIOR2['norm_beta'] <- dnorm(PRIOR2$probability, mean_post.beta, sqrt(var_post.beta))

unif2 <- ggplot(PRIOR2, aes(x=probability))+
	 geom_line(aes(y=post.unif, color='Uniform'))+
	 geom_line(aes(y=norm_unif, color='Norm approx.'))+
	 labs(x='p', y='Probability', title='UNIFORM', colour = "")
	 
beta2 <- ggplot(PRIOR2, aes(x=probability))+
	 geom_line(aes(y=post.beta, color='Beta(1,4)'))+
	 geom_line(aes(y=norm_beta, color='Norm approx.'))+
	 labs(x='p', y='Probability', title='Beta (1,4)', colour = "")
	 
pdf("Norm_approximation.pdf")
grid.arrange(unif2, beta2, nrow = 2)
out<-dev.off()

cat('----------------------------------------------------------', '\n')	
cat('2.2: See the following table', '\n', '\n')

# function that allows us to find the probability releated to a given value of the cumulative distribution
percentile <- function(k, input, output) {			
		# k = value of the cumulative distribution
		# output = probability distribution
		# input = relative cumulative distribution
		index <- which.min(abs(input - k))	
		return(output[index])
		}
		
per_2.5_unif <- percentile(0.025, cumsum(PRIOR2$post.unif/10**4), PRIOR2$probability)
per_2.5_beta <- percentile(0.025, cumsum(PRIOR2$post.beta/10**4), PRIOR2$probability)
per_2.5_norm_unif <- percentile(0.025, cumsum(PRIOR2$norm_unif/10**4), PRIOR2$probability)
per_2.5_norm_beta <- percentile(0.025, cumsum(PRIOR2$norm_beta/10**4), PRIOR2$probability)

per_97.5_unif <- percentile(0.975, cumsum(PRIOR2$post.unif/10**4), PRIOR2$probability)
per_97.5_beta <- percentile(0.975, cumsum(PRIOR2$post.beta/10**4), PRIOR2$probability)
per_97.5_norm_unif <- percentile(0.975, cumsum(PRIOR2$norm_unif/10**4), PRIOR2$probability)
per_97.5_norm_beta <- percentile(0.975, cumsum(PRIOR2$norm_beta/10**4), PRIOR2$probability)

summary2 <- tibble(Associated_Prior=c('Uniforn', 'Beta(1,4)', 'Norm Unif', 'Norm Beta'),
		  Mean=c(mean_post.unif, mean_post.beta, PRIOR2$probability[which.max(PRIOR2$norm_unif)], PRIOR2$probability[which.max(PRIOR2$norm_beta)]),
		  perc_2.5=c(per_2.5_unif, per_2.5_beta, per_2.5_norm_unif, per_2.5_norm_beta),
		  perc_97.5=c(per_97.5_unif, per_97.5_beta, per_97.5_norm_unif, per_97.5_norm_beta))

print(summary2)

cat('----------------------------------------------------------', '\n')	
cat('2.2: See Final_plot2.pdf', '\n', '\n')

unif <- ggplot(PRIOR2, aes(x=probability))+
	geom_line(aes(y=post.unif))+
	labs(x='p', y='Probability', title='Uniform')+
	geom_vline(xintercept = per_2.5_unif, linetype="dotted", color = "blue")+
	geom_vline(xintercept = per_97.5_unif, linetype="dotted", color = "blue")
	
beta <- ggplot(PRIOR2, aes(x=probability))+
	geom_line(aes(y=post.beta))+
	labs(x='p', y='Probability', title='Beta(1,4)')+
	geom_vline(xintercept = per_2.5_beta, linetype="dotted", color = "blue")+
	geom_vline(xintercept = per_97.5_beta, linetype="dotted", color = "blue")
	
nor_unif <- ggplot(PRIOR2, aes(x=probability))+
	    geom_line(aes(y=norm_unif))+
	    labs(x='p', y='Probability', title='Normal uniform approximation')+
	    geom_vline(xintercept = per_2.5_norm_unif, linetype="dotted", color = "blue")+
	    geom_vline(xintercept = per_97.5_norm_unif, linetype="dotted", color = "blue")
	  
nor_beta <- ggplot(PRIOR2, aes(x=probability))+
	    geom_line(aes(y=norm_beta))+
	    labs(x='p', y='Probability', title='Normal beta(1,4) approximation')+
	    geom_vline(xintercept = per_2.5_norm_beta, linetype="dotted", color = "blue")+
	    geom_vline(xintercept = per_97.5_norm_beta, linetype="dotted", color = "blue")		

pdf("Final_plot2.pdf")
grid.arrange(unif, beta, nor_unif, nor_beta, nrow = 2)	
out<-dev.off()
cat('----------------------------------------------------------', '\n', '\n')

cat(' #########  ESERCIZIO 3   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('3.1: See Toss_a_coin.pdf', '\n')

testa = 15
croce = 15

# Creating a table containing all discrete values of our pdf 	
PRIOR3 <- tibble(probability=seq(1/10**4, 1-1/10**4, length.out = 10**4))
PRIOR3['unif'] <- dbeta(PRIOR3$probability,1,1)
PRIOR3['beta'] <- dbeta(PRIOR3$probability,2,2)

prior3 <- ggplot(PRIOR3, aes(x=probability))+
	 geom_line(aes(y=beta, color='Beta(2, 2)'))+
	 geom_line(aes(y=unif, color='Uniform'))+
	 labs(x='p', y='Probability', title='PRIOR', colour = "Priors")

PRIOR3['likelihood'] <- dbinom(x=testa, size=testa+croce, prob=PRIOR3$probability)

like3 <- ggplot(PRIOR3, aes(x=probability))+
	 geom_line(aes(y=likelihood, color='Binomial'))+
	 labs(x='p', y='Probability', title='LIKELIHOOD', colour = "Likelihood")
		
post <- function(prior) {						# function that given the prior gives us the normalized posterior
	posterior <- PRIOR3$likelihood * prior
	k <- 1/10**4*sum(posterior)					# calculate the intergral of the posterior
	posterior <- 1/k * posterior					# normalize the posterior
	return(posterior)
	}		
	
PRIOR3['post.unif'] <- post(PRIOR3$unif)
PRIOR3['post.beta'] <- post(PRIOR3$beta)	

posterior3 <- ggplot(PRIOR3, aes(x=probability))+
	     geom_line(aes(y=post.beta, color='Beta(2, 2)'))+
	     geom_line(aes(y=post.unif, color='Uniform'))+
	     labs(x='p', y='Probability', title='POSTERIOR', colour = "Posteriors")

pdf("Toss_a_coin.pdf")
grid.arrange(prior3, like3, posterior3, nrow = 3)
out<-dev.off()

cat('----------------------------------------------------------', '\n')	
cat('3.2: Summary table', '\n')

# function that allows us to find the probability releated to a given value of the cumulative distribution
percentile <- function(k, input, output) {			
		# k = value of the cumulative distribution
		# output = probability distribution
		# input = relative cumulative distribution
		index <- which.min(abs(input - k))	
		return(output[index])
		}
		
PRIOR3['cum_post.unif'] <- cumsum(PRIOR3$post.unif/10**4)
PRIOR3['cum_post.beta'] <- cumsum(PRIOR3$post.beta/10**4)

per_2.5_unif <- percentile(0.025, PRIOR3$cum_post.unif, PRIOR3$probability)
per_2.5_beta <- percentile(0.025, PRIOR3$cum_post.beta, PRIOR3$probability)

per_97.5_unif <- percentile(0.975, PRIOR3$cum_post.unif, PRIOR3$probability)
per_97.5_beta <- percentile(0.975, PRIOR3$cum_post.beta, PRIOR3$probability)
	
# Creating a summaring table	
summary3 <- tibble(Associated_Prior=c('Uniforn', 'Beta(2,2)'),
		  MPV=c(PRIOR3$probability[which.max(PRIOR3$post.unif)], PRIOR3$probability[which.max(PRIOR3$post.beta)]),
		  perc_2.5=c(per_2.5_unif, per_2.5_beta),
		  perc_97.5=c(per_97.5_unif, per_97.5_beta))
		
print(summary3)
		
cat('----------------------------------------------------------', '\n')	
cat('3.3: Summary sequential table', '\n')		
		
prior <- dbeta(PRIOR3$probability,1,1)
toss <- c('T', 'T', 'T', 'T', 'T', 'H', 'T', 'T', 'H', 'H', 'T', 'T', 'H', 'H', 'H', 'T', 'H', 'T', 'H', 'T', 'H', 'H', 'T', 'H', 'T', 'H', 'T', 'H', 'H', 'H')
lancio <- 0

# The following elements are list, after each toss the informations about the toss will be append to these lists.
tentativo <- 0 
most_prob_value <- 0.5
low_limit <- 0
upper_limit <- 1
estrazione <- '-'
 
for (r in toss) {
	lancio <- lancio+1 
	if (r=='T') {T<-1}
	else {T<-0}
	
	# Here I compute the posterior and then I normalize it
	lh <- dbinom(x=T, size=1, prob=PRIOR3$probability)
	posterior <- lh * prior
	k <- 1/10**4*sum(posterior)
	posterior <- 1/k * posterior
	
	# Find the MPV 
	most_probable_value <- PRIOR3$probability[which.max(posterior)]	
	
	cum_seq_post <- cumsum(posterior/10**4)
	
	# Find the 2.5% and 97.5% limits
	per_2.5_seq <- percentile(0.025, cum_seq_post, PRIOR3$probability)
	per_97.5_seq <- percentile(0.975, cum_seq_post, PRIOR3$probability)
	
	# Update the lists with the information releated to this toss 
	tentativo <- c(tentativo, lancio)
	estrazione <- c(estrazione, r)
	most_prob_value <- c(most_prob_value, most_probable_value)
	low_limit <- c(low_limit, per_2.5_seq)
	upper_limit <- c(upper_limit, per_97.5_seq)
	
	# In this way the current posterior will become the prior for the next toss
	prior <- posterior
	}

# Summaring table
SEQUENTIAL <- tibble(Toss_number=tentativo, result=estrazione, MPV=most_prob_value, perc_2.5=low_limit, perc_97.5=upper_limit)
print(SEQUENTIAL, n=31)

cat('----------------------------------------------------------', '\n')	
cat('3.4: If we compare', '\n', '\n')
print(tail(SEQUENTIAL, 1))
cat('\n', 'with', '\n', '\n')
print(summary3)
cat('\n', 'we can notice that we obtain the same result.', '\n')
cat('----------------------------------------------------------', '\n', '\n')
	

cat(' #########  ESERCIZIO 4   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('4.1 - 4.2: Use the function simulation() to choose a random box and to sample from it.', '\n')

simulation <- function(n) {
	q <- 6*runif(1)
	# q variable allows me to select one box and then
	# we define the value of the variable prob according to the box selected 
	if (q<=1) {				# H0
		prob<-1
		box<-'H0  (B B B B B)'}			
	else if (1<=q & q<=2) {			# H1
		prob<-0.8
		box<-'H1  (B B B B W)'}			
	else if (2<=q & q<=3) {			# H2
		prob<-0.6
		box<-'H2  (B B B W W)'}			
	else if (3<=q & q<=4) {			# H3
		prob<-0.4
		box<-'H3  (B B W W W)'}			
	else if (4<=q & q<=5) {			# H4
		prob<-0.2
		box<-'H4  (B W W W W)'}			
	else if (5<=q) {			# H5
		prob=0
		box<-'H5  (W W W W W)'}			
	
	# now we can simulate the extraction of n stone from the ramdomly selected box
	samples <- x <- rbinom(n,1,prob)	# 1 rappresenta la pallina nera
	cat('The chosen box is :', box, '\n')
	cat('The extraction from the box is:', samples, '(1=black, 0=white)', '\n')
	return(samples) 	
	}
	
cat('----------------------------------------------------------', '\n')	
cat('4.3: ', '\n')	
analisi <- function(sampled) {
	# the following varaibles will be update after each extraction
	W <- 0
	B <- 0
	# the following variables are list, after each extraction the new informations will be append to these lists.
	P_H0 <- 1/6
	P_H1 <- 1/6
	P_H2 <- 1/6
	P_H3 <- 1/6
	P_H4 <- 1/6
	P_H5 <- 1/6
	for (s in sampled) {
		# updating the number of white and black stones extracted from the box
		if (s==1) {B<-B+1}
		else {W<-W+1}
		
		tot <- B+W
		
		# computing the prior for each box given the number of black and white stone extracted
		P_E_H0 <- dbinom(x=B, size=tot, prob=1)
		P_E_H1 <- dbinom(x=B, size=tot, prob=0.8)
		P_E_H2 <- dbinom(x=B, size=tot, prob=0.6)
		P_E_H3 <- dbinom(x=B, size=tot, prob=0.4)
		P_E_H4 <- dbinom(x=B, size=tot, prob=0.2)
		P_E_H5 <- dbinom(x=B, size=tot, prob=0)
		
		evidence <- (P_E_H0+P_E_H1+P_E_H2+P_E_H3+P_E_H4+P_E_H5)*(1/6)
		
		# computing the posterior and updating the previous variable 	
		P_H0 <- c(P_H0, P_E_H0*(1/6)/evidence)
		P_H1 <- c(P_H1, P_E_H1*(1/6)/evidence)
		P_H2 <- c(P_H2, P_E_H2*(1/6)/evidence)
		P_H3 <- c(P_H3, P_E_H3*(1/6)/evidence)
		P_H4 <- c(P_H4, P_E_H4*(1/6)/evidence)
		P_H5 <- c(P_H5, P_E_H5*(1/6)/evidence)
		}
	
	# pur all the result in a tibble
	result <- tibble(Draw=seq(0,length(sampled)))
	result['Result'] <- c('-', sampled)
	result['Box_H0'] <- P_H0
	result['Box_H1'] <- P_H1
	result['Box_H2'] <- P_H2
	result['Box_H3'] <- P_H3
	result['Box_H4'] <- P_H4
	result['Box_H5'] <- P_H5
	return(result)
	}


num_extraction = 30

lanci <- simulation(num_extraction)
risultati <- analisi(lanci)

cat('The result of the simulation of', num_extraction, 'extractions is:' ,'\n', '\n')
print(risultati, n=length(lanci)+1)

cat('----------------------------------------------------------', '\n')	
cat('4.4: See Extraction_simulation.pdf', '\n')

pH0 <- ggplot(risultati, aes(x=Draw))+
       geom_point(aes(y=Box_H0))+
       geom_line(aes(y=Box_H0))+
       ylim(0,1)+
       labs(x='number of extraction', y='Probability', title='Box H0')
     
pH1 <- ggplot(risultati, aes(x=Draw))+
       geom_point(aes(y=Box_H1))+
       geom_line(aes(y=Box_H1))+
       ylim(0,1)+
       labs(x='number of extraction', y='Probability', title='Box H1')

pH2 <- ggplot(risultati, aes(x=Draw))+
       geom_point(aes(y=Box_H2))+
       geom_line(aes(y=Box_H2))+
       ylim(0,1)+
       labs(x='number of extraction', y='Probability', title='Box H2')
       
pH3 <- ggplot(risultati, aes(x=Draw))+
       geom_point(aes(y=Box_H3))+
       geom_line(aes(y=Box_H3))+
       ylim(0,1)+
       labs(x='number of extraction', y='Probability', title='Box H3')
       
pH4 <- ggplot(risultati, aes(x=Draw))+
       geom_point(aes(y=Box_H4))+
       geom_line(aes(y=Box_H4))+
       ylim(0,1)+
       labs(x='number of extraction', y='Probability', title='Box H4')
       
pH5 <- ggplot(risultati, aes(x=Draw))+
       geom_point(aes(y=Box_H5))+
       geom_line(aes(y=Box_H5))+
       ylim(0,1)+
       labs(x='number of extraction', y='Probability', title='Box H5')

pdf("Extraction_simulation.pdf")
grid.arrange(pH0, pH1, pH2, pH3, pH4, pH5, nrow = 2)
out<-dev.off()
cat('----------------------------------------------------------', '\n')	

