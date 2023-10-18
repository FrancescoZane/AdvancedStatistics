library(tibble)
library(dplyr)
library(lubridate)
library(ggplot2)

										# 		EXERCISE 1
										
										
									# read the data and import them in a data.frame or tibble structure
american_airline <- read.table('american_airline_empl.txt', sep = '\t', header = TRUE)
delta_airline <- read.table('delta_airline_empl.txt', sep = '\t', header = TRUE)
federal_express <- read.table('federal_express_empl.txt', sep = '\t', header = TRUE)
united_airline <- read.table('united_airline_empl.txt', sep = '\t', header = TRUE)

american_airline$Company = 'american_airline'
delta_airline$Company = 'delta_airline'
federal_express$Company = 'federal_express'
united_airline$Company = 'united_airline'

									# merge the four data tibble in a common tibble
final_data = as_tibble(bind_rows(american_airline, delta_airline, federal_express, united_airline))              # I bind them in vertical way

									# produce a plot of the behaviour of the employees as a function of time for all four companies
final_data$Date = paste(as.character(final_data$Year), as.character(final_data$Month), sep='-') |> ym()          # convert into character and then create dates

final_data$Full.time<-as.numeric(gsub(",", "", final_data$Full.time))
final_data$Part.time<-as.numeric(gsub(",", "", final_data$Part.time))
final_data$Grand.Total<-as.numeric(gsub(",", "", final_data$Grand.Total))

pdf("Part-time.pdf")
ggplot(final_data, aes(x=Date, y=Part.time, color=Company))+
geom_line()+
labs(x='Year', y='Number of employees', title='Part-time employees')+
scale_color_manual(values=c('red', 'yellow', 'green', 'blue'))
dev.off()

pdf("Full-time.pdf")
ggplot(final_data, aes(x=Date, y=Full.time, color=Company))+
geom_line()+
labs(x='Year', y='Number of employees', title='Full-time employees')+
scale_color_manual(values=c('red', 'yellow', 'green', 'blue'))
dev.off()

									# when did each company reach the minimum and maximum number of employess ?
max_emp <- aggregate(Grand.Total ~ Company, data=final_data, FUN=max)   # for each companies I find the max number of employees
min_emp <- aggregate(Grand.Total ~ Company, data=final_data, FUN=min)   # for each companies I find the min number of employees

max_emp_date <- merge(max_emp, final_data, by=c('Grand.Total', 'Company'))  
min_emp_date <- merge(min_emp, final_data, by=c('Grand.Total', 'Company'))
print('When companies reach the max number of employees: ')
max_emp_date
print('When companies reach the min number of employees: ')
min_emp_date

									# plot the fraction of part-time worker over the total employess as a function of time
final_data$ratio <- final_data$Part.time/final_data$Grand.Total

pdf("Part-time over all workers.pdf")
ggplot(final_data, aes(x=Date, y=ratio, color=Company))+
geom_line()+
labs(x='Year', y='Ratio', title='Part-time over Grand-total employees')+
scale_color_manual(values=c('red', 'yellow', 'green', 'blue'))
dev.off()

									# did the COVID-19 pandemic have any influence in the employed workers of the airline companies
covid_time <- final_data |> filter( Year>=2019, Year<=2023 )

pdf("Covid.pdf")
ggplot(covid_time, aes(x=Date, y=Grand.Total, color=Company))+
geom_line()+
labs(x='Year', y='Number of employees', title='Total number of employees during covid')+
scale_color_manual(values=c('red', 'yellow', 'green', 'blue'))
dev.off()

print('I can visualise that for most of companies the number of employees dropped when the covid pandemic started')







										# 		EXERCISE 2

library('nycflights13')

							# Plot the total number of flights departed from each of the three NYC airports as a function of time
flights
number <- flights |> group_by(month, day, origin) |> tally()
number$year <- 2013
number$Date <- paste(as.character(number$year), as.character(number$month), as.character(number$day), sep='-') |> ymd()

pdf("Departures.pdf")
ggplot(number, aes(x=Date, y=n, color=origin))+
geom_line()+
labs(x='Date', y='Number of flights', title='Number of flights over the year')+
scale_color_manual(values=c('red', 'green', 'blue'))
dev.off()


							# Plot the average number of flights computed over the first five working days of each week as a function of the
							# week number of the year. Produce the same plot for the flights departing over the weekend (Saturdays and Sundays).
voli <- number |> group_by(month, day) |> summarise(numero=sum(n))
voli$year <- 2013
voli$Date <- paste(as.character(voli$year), as.character(voli$month), as.character(voli$day), sep='-') |> ymd()
voli$w <- isoweek(voli$Date)
voli$d <- wday(voli$Date, week_start=1)
voli$week <- 'Working day'
voli$week[voli$d==7 | voli$d==6] <- 'Week-end'

week_fly <- voli |> group_by(w, week) |> summarise(voli=mean(numero))

pdf("Flightsoverweek.pdf")
ggplot(week_fly, aes(x=w, y=voli, color=week))+
geom_line()+
labs(x='Weeks', y='Number of flights', title='Number of flights over the weeks')+
scale_color_manual(values=c('red', 'blue'))
dev.off()



							# extract the following pieces of information (separately for each NYC airport):
							# - min, max and average delay for each day of the year (show the data in corresponding plots)
not_cancelled <- flights |> filter( !is.na(dep_delay), !is.na(arr_delay) )   # delate the cancelled flights that are caractired by nan values
delay <- not_cancelled |> group_by(month, day, origin) |> summarize(media=mean(dep_delay), max=max(dep_delay), min=min(dep_delay))
delay$year <- 2013
delay$Date <- paste(as.character(delay$year), as.character(delay$month), as.character(delay$day), sep='-') |> ymd()


pdf("Mean_delay.pdf")
ggplot(delay, aes(x=Date, y=media, color=origin))+
geom_line()+
labs(x='Date', y='Minute', title='Mean delay over the year')+
scale_color_manual(values=c('red', 'blue', 'green'))
dev.off()

pdf("Max_delay.pdf")
ggplot(delay, aes(x=Date, y=max, color=origin))+
geom_line()+
labs(x='Date', y='Minute', title='Max delay over the year')+
scale_color_manual(values=c('red', 'blue', 'green'))
dev.off()

pdf("Min_delay.pdf")
ggplot(delay, aes(x=Date, y=min, color=origin))+
geom_line()+
labs(x='Date', y='Minute', title='Min delay over the year')+
scale_color_manual(values=c('red', 'blue', 'green'))
dev.off()


							# assuming the distance flew by the plane is, at first approximation, the distance between the two
							# connecting airports (as given in the data frame), compute the average speed of each plane. Produce
							# a plot of the average plane speed as a function of departure day of the year
speed <- filter(not_cancelled) |> group_by(month, day) |> summarize(speed_mean=mean(distance/air_time)*96.56)
speed$year <- 2013
speed$Date <- paste(as.character(speed$year), as.character(speed$month), as.character(speed$day), sep='-') |> ymd()

pdf("Speed.pdf")
ggplot(speed, aes(x=Date, y=speed_mean))+
geom_line()+
labs(x='Date', y='[Km/h]', title='Average plane speed over the year')
dev.off()

							# analyze the flights offered by each airline company and determine:
							# - the airline companies offering the largest two numbers of flights per day and per week;
flights_companies <- flights |> group_by(carrier, month, day) |> tally()
print('The airline companies offering the largest two numbers of flights per day are: ')
flights_companies |> arrange(desc(n)) |> group_by(month, day) |> slice(1:2)							
							
flights_companies$year <- 2013
flights_companies$week <- paste(as.character(flights_companies$year), as.character(flights_companies$month), as.character(flights_companies$day)) |> ymd() |> isoweek()
week_flights <- flights_companies |> group_by(carrier, week) |> summarise(total=sum(n))
print('The airline companies offering the largest two numbers of flights per week are: ')
week_flights |> arrange(desc(total)) |> group_by(week) |> slice(1:2)							
							
							
							# - the airline company offering the smallest number of flight per month;
month_flights <- flights_companies |> group_by(carrier, month) |> summarise(total=sum(n))
print('The airline company offering the smallest number of flight per month is: ')
month_flights |> arrange(desc(total)) |> group_by(month) |> slice(which.min(total))							
							
							
							# - the airline company offering the longest distance flight per month;
print('The airline company offering the longest distance flight per month is: ')
flights |> group_by(carrier, month) |> summarise(dist=max(distance)) |> group_by(month) |> top_n(1, dist)
