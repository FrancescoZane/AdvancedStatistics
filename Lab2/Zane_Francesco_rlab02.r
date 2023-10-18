library(tibble)
library(dplyr)
library(ggplot2)
set.seed(061217)
	
	
	
											###   ESERCIZIO 1   ###
cat(' #########  ESERCIZIO 1   #########', '\n', '\n')	
cat('----------------------------------------------------------', '\n')	
cat('1.1 - 1.2 See ddisc.pdf for the pdf and pdisc.pdf for the cdf', '\n')								
# definisco la pdf	
ddisc <- function(x) {
	result <- ifelse((x>=0 & x<6), (x%/%1)/15, 0)
	return(result)
	}

PDF <- tibble(step=seq(0,7,1), value=ddisc(seq(0,7,1)))

pdf("ddisc.pdf")
ggplot(PDF, aes(x=step, y=value))+
geom_step()+
geom_point()+
labs(x='x', y='Probability', title='Probability density function')
out<-dev.off()

# definisco la cdf
pdisc <- function(x) {
	y = 0
	x[x<0] <- 0
	for (k in seq(1,max(x),1)) {y <- c(y, sum(ddisc(seq(0,k,1))))}   # qui sto creando il vettore [1, 2, 3, ..., max(x)] e k varia all'interno di questo vettore
									 # per ogni valore di k aggiungo un nuovo elemento al vettore y pari alla cumulativa della pdf fino al valore k
									 # in questo modo gli elementi del vettore y sono le cumulative dei vari elementi del vettore x
	return(y[x+1])							 # utilizzo il valore x come indice (per il vettore y) per ottenere la giusta cmulativa a lui associata
	}
	

CDF <- tibble(step=seq(0,7,1), value=pdisc(seq(0,7,1)))

pdf("pdisc.pdf")
ggplot(CDF, aes(x=step, y=value))+
geom_step()+
geom_point()+
labs(x='x', y='Probability', title='Cumulative density function')
out<-dev.off()
cat('----------------------------------------------------------', '\n')	

# calcolo la media utilizzando la pdf e quindi come somma dei vari valori discreti moltiplicati per la loro probabilita' di uscire
mean <- 0
for (i in seq(1,5,1)) {
	mean <- mean + i*ddisc(i)
	}
cat('1.3', 'The mean: ', mean, '\n')
	
# analogamente calcolo anche la varianza
var <- 0
for (i in seq(1,5,1)) {
	var <- var + (i-mean)**2*ddisc(i)
	}
cat('    The variance: ', var, '\n')
cat('----------------------------------------------------------', '\n')	

fun <- 0
for (i in seq(1,5,1)) {
	fun <- fun + (i*(6-i))*ddisc(i)
	}
cat('1.4 The expected value for the function is: ', fun, '\n')
cat('----------------------------------------------------------', '\n')	

samples_1 <- function(n) {
	sample(seq(0,5), size=n, replace = TRUE, prob = ddisc(seq(0,5)))
	}	

sample_data <- tibble(value=samples_1(10**5))
cat('1.5 Create function called: samples_1', '\n')
cat('----------------------------------------------------------', '\n')

cat('1.6 See simalation.pdf', '\n')	
pdf("simulation.pdf")
ggplot()+
geom_histogram(data=sample_data, aes(x=value, after_stat(density)), binwidth=1, color="darkblue", fill="lightblue")+
geom_point(data=PDF, aes(x=step, y=value, color='theoretical pdf'))+
labs(x='x', y='Probability', title='Simulated samples')
out<-dev.off()
cat('----------------------------------------------------------', '\n', '\n')


										###   ESERCIZIO 2   ###
cat('\n', '#########  ESERCIZIO 2   #########', '\n', '\n')
cat('----------------------------------------------------------', '\n')	
# definisco la pdf											
dtrian <- function(a, b, c, x) {
	result <- ifelse((x>=a & x<=b),
		ifelse((x>=a & x<c),
			2*(x-a)/((b-a)*(c-a)),
			2*(b-x)/((b-a)*(b-c)))
		, 0)
	return(result)
	}
	
a <- 0
b <- 4
c <- 1
x <- seq(a-0.5, b+0.5, 0.001)

PDF <- tibble(x=x, y=dtrian(a, b, c, x))
pdf("pdf_triang.pdf")
ggplot(PDF, aes(x=x, y=y))+
geom_line()+
labs(x='x', y='Probability', title='Probability density function')
out<-dev.off()
cat('2.1 See pdf_triang.pdf', '\n')
cat('----------------------------------------------------------', '\n')	

samples_2 <- function(a, b, c, n) {			# per campionare i dati dalla distribuzione triandolare utilizzo l'inverse method
	th <- (c-a)**2/((b-a)*(c-a))			# essendo la pdf una funzione composta anche la cdf che volgio inventire e' una funzione composta, th e' quindi il valore della cdf calcolata nel punto
							# c e mi serve per decidere quale delle due funzioni utilizzare quando devo invertire la cdf
	v <- runif(n)
	result <- ifelse(v<=th, a + sqrt(v * (b-a) * (c-a)), b - sqrt((1-v) * (b-a) * (b-c)) )      # qui sto invertendo la cdf e per l'appunto utilizzo th come discriminante
	return(result)
	}
cat('2.2 The algorithm is called: samples_2', '\n')	
cat('----------------------------------------------------------', '\n')	


sample_data <- tibble(value=samples_2(a,b,c,10**4))

pdf("Simulation_triang.pdf")
ggplot()+
geom_histogram(data=sample_data, aes(x=value, after_stat(density)), binwidth=0.1, color="darkblue", fill="lightblue")+
geom_line(data=PDF, aes(x=x, y=y, color='theoretical pdf'))+
labs(x='x', y='Probability', title='Simulated samples')
out<-dev.off()
cat('2.3 See Simulation_triang.pdf', '\n')
cat('----------------------------------------------------------', '\n', '\n')	


										###   ESERCIZIO 3   ###
cat('\n', '#########  ESERCIZIO 3   #########', '\n', '\n')
cat('----------------------------------------------------------', '\n')	
lambda = 1/30  # in minutes

# as pdf for exp distribution I use:    dexp(x, rate = 1, log = FALSE)

cdf <- function(t) {
	y <- 1 - exp(-lambda * t) 
	return(y)
	}

samples_3 <- function(n) {              # per campionare i dati dalla distribuzione triandolare utilizzo l'inverse method
	v <- runif(n)
	k <- numeric()
	k <- - (log(v)) / (lambda)    # since for uniform random numbers generated 1-u=u
	return(k)
	}

sample_data <- tibble(value=samples_3(60))
interval <- seq(0, max(sample_data$value), 0.1)
PDF <- tibble(time=interval, prob=dexp(interval, rate = lambda, log = FALSE))

pdf("Simulation_doctor.pdf")
ggplot()+
geom_histogram(data=sample_data, aes(x=value, after_stat(density)), binwidth=5, color="darkblue", fill="lightblue")+
geom_line(data=PDF, aes(x=time, y=prob, color='theoretical pdf'))+
labs(x='Time [minute]', y='Probability', title='Simulated samples')
out<-dev.off()
cat('3.1 See Simulation_doctor.pdf', '\n')
cat('----------------------------------------------------------', '\n')	

p_time_less_12 <- pexp(12, rate=lambda, log = FALSE)
cat('3.2 The probability that a person will wait for less than 12 minutes: ', p_time_less_12, '\n')
cat('----------------------------------------------------------', '\n')	

f <- function(t) {
	y <- t * lambda * exp(-lambda * t) 
	return(y)
	}
	
manipulation <- integrate(f, 0, 1000, subdivisions = 10**6)$value
exp_mean <- mean(sample_data$value)

cat('3.3 Average waiting time (in minutes):', '\n')
cat('    - theory = 1/lambda = 30', '\n')
cat('    - calculation from the pdf = ', manipulation, '\n')
cat('    - mean of sampled data = ', exp_mean, '\n')
cat('----------------------------------------------------------', '\n')	

p_time_more_60 <- 1-pexp(60, rate=lambda, log = FALSE)
cat('3.4 The probability for waiting more than one hour before being received: ', p_time_more_60, '\n')
cat('----------------------------------------------------------', '\n', '\n')	


										###   ESERCIZIO 4   ###
cat('\n', '#########  ESERCIZIO 4   #########', '\n', '\n')
cat('----------------------------------------------------------', '\n')	
# legend 
# sk = student knows 				sdk = student dosen't know
# ac = answer is correct			aw = answer is wrong
# p__...  means prob that ...			p__..._~~  means prob that ... happens after ~~

p__sk <- 0.7
p__ac_sk <- 1
p__ac_sdk <- 0.2



prior <- p__sk
likehood <- p__ac_sk
evidence <- p__ac_sk*p__sk+ p__ac_sdk*(1-p__sk)

p__sk_ac <- prior * likehood / evidence
cat('Once a correct answer is given, the probability that the student really knew the correct answer is: ', p__sk_ac, '\n')
cat('----------------------------------------------------------', '\n', '\n')	

										###   ESERCIZIO 5   ###
cat('\n', '#########  ESERCIZIO 5   #########', '\n', '\n')
cat('----------------------------------------------------------', '\n')	
a <- 0
b <- 60
sam <- a+runif(10**6)*(b-a)        # here I simulate the arrival of 10**6 passengers, the numbers represent the minute in which tjey arrive after the 10:45 (then for example 10 menas 10:55)
passengers <- tibble(time=sam )

hist(passengers$time, breaks=seq(0,60,1), xlab='Arrival time (minute)', ylab='Counts', main='Passengers arrival time')

delay <- function(x) {		   # this function links the moment in which a passenger arrives with the number of minutes he has to wait for the next train
	result <- ifelse(x<=15,
			-x+15,
			ifelse(x<=45,
				-x+45,
				-x+75))
	return(result)
	}
		
passengers['ritardo'] <- delay(passengers$time)

isto<-hist(passengers$ritardo, breaks=seq(0,30,1), xlab='Delay time (minute)', ylab='Counts', main='Time passengers have to wait')
# ora questo grafico verra' usato come pdf sperimentale per estrarre le informazioni richieste 


percent <- function(x, v) {                         # presa una caratteristica/valore per ogni colonna dell'istogramma, questa funzione permetter di sommare questi valori a partire dalla prima colonna fino 
						    # a quella desiderata ed indicata dal valore x, nel nostro caso quello che andremo a sommare sara' la densita' ovvero la probabilita'
	s <- 0
	for (i in seq(1,x,1)) {
		s <- s+v[i]
		}
	return(s)
	}
	
wait_most_10 <- percent(10,isto$density)
cat('5.1 Probability that she has to wait at most 10 minutes: ', wait_most_10, '\n')
cat('----------------------------------------------------------', '\n')	

wait_least_15 <- 1-percent(15,isto$density)
cat('5.2 Probability that she has to wait at least 15 minutes: ', wait_least_15, '\n')
cat('----------------------------------------------------------', '\n')	

average <- 0                                    
for (i in seq(1,30,1)) {   			    # qui calcolo la media usando l'istogramma prima creato come distribuzione di probabilita' discreta
	average <- average+i*isto$density[i]
	}						

cat('5.3 The average time (in minute) spent waiting is: ', average, '\n')
cat('----------------------------------------------------------', '\n', '\n')	

										###   ESERCIZIO 6   ###
cat('\n', '#########  ESERCIZIO 6   #########', '\n', '\n')
cat('----------------------------------------------------------', '\n')	
mean <- 0.10
sd <- 0.12
n <- 200
price <- 85

tot <- n * price
mean_ex <- mean*tot
sd_ex <- sd*tot

prob_atleast_800 <- 1 - pnorm(800, mean_ex, sd_ex)
cat('6.1 The probability that after a year his net profit from the investment is at least 800 euro is:', prob_atleast_800, '\n')
cat('----------------------------------------------------------', '\n')	




