\name{micSim}
\alias{micSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Run microsimulation (sequentially)
}
\description{
 Performs a continuous-time microsimulation run (sequentially, i.e., using only one CPU core).
}
\usage{
micSim(initPop, immigrPop = NULL, transitionMatrix, absStates = NULL, 
  initStates = c(), initStatesProb = c(), maxAge = 99, simHorizon, 
  fertTr = c(), dateSchoolEnrol="09/01")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{initPop}{
    Data frame comprising the starting population of the simulation.
}
  \item{immigrPop}{
    Data frame comprising information about the immigrants entering the population across simulation time. 
}
  \item{transitionMatrix}{
  A matrix indicating the transition pattern and the names of the functions determining the respective transition rates.
}
  \item{absStates}{
 A vector indicating the absorbing states of the model.
}
  \item{initStates}{
 A vector comprising all initial states that newborns might enter.
}
  \item{initStatesProb}{
 A vector comprising the probabilities corresponding to \code{initStates}. In sum, these probabilities have to be one.
}
  \item{maxAge}{
 A scalar indicating the maximal age which an individual can reach during simulation. \code{maxAge} has to be greater than zero
}
  \item{simHorizon}{
 A vector comprising the starting and ending date of the simulation. Both dates have to be chron objects. The starting date has to precede the ending date.
}
  \item{fertTr}{
A vector indicating all transitions triggering a child birth event during simulation, that is, the creation of a new individual. 
}
  \item{dateSchoolEnrol}{
A string of the form "month/day" indicating the general enrollment date for elementary school, e.g., "09/01" for September 1st. The default setting is "09/01". 
}
}
\details{
All nonabsorbing states considered during simulation have to be defined as composite states. In more detail, they consist of labels indicating values of state variables. Within states, labels are separated by a forward slash "/". Possible state variables are, for example, gender, number of children ever born, and educational attainment. Corresponding values are, for example, "m" and "f" (gender), "0","1","2", and "3+" (number of children ever born), "no", "low", "med", and "high" (educational attainment). Possible examples of states are "m/0/low" for a childless male with elementary education  or "f/1/high" for a female with one child and a higher secondary school degree. All state variables considered plus accordant value labels have to be provided by the user. The only exception is gender which is predefined by labels "m" and "f" indicating male and female individuals. The label values "no" and "low" are reserved for enrolment events to elementary school (see below).

Nonabsorbing states have to be given as strings such as "dead" for being dead or "rest" for emigrated.

\code{micSim} is able to conduct enrollment events to elementary school such that they take place on \code{dateSchoolEnrol} of a particular year. For this purpose, a state variable defining educational attainment has to be created first. Then, labels of possible values have to be defined such that "no" describes no education and "low" describes elementary education. Finally, the transition function determining the transition rate for the respective enrollment event has to be defined to return "Inf" for the age x at which children should be enrolled (e.g., at age seven) and zero otherwise. That way, an event "school enrollment on \code{dateSchoolEnrol} of the year in which a child turns x years old" is enforced.

If educational attainment is not considered, \code{dateSchoolEnrol} can let be unspecified: \code{dateSchoolEnrol=c()}. 

The starting population \code{initPop} has to be given in the form of a data frame. Each row of the data frame corresponds to one individual. \code{initPop} has to comprise the following information: 
unique numerical person identifier (ID), birth date, and initial state (i.e., the state occupied by the individual when entering the synthetic population). Birth dates have to be chron objects. 

Information about immigrants has to be given in the form of a data frame (\code{immigrPop}). Each row of the data frame corresponds to one immigrant.
\code{immigrPop} contains the following data: unique numerical person identifier (ID), immigration date, birth date, and initial state (i.e., the state occupied by the immigrant when entering the simulated population). Immigration dates and birth dates have to be chron objects. 

For each transition that should be considered during simulation accordant transition rates have to be provided. \code{micSim} requires these rates in form of functions which are handed over via the transition matrix \code{transitionMatrix} (described in the subsequent paragraph). The \code{MicSim} package allows rates to depend on three time scales: age, calendar time, and the time that has elapsed since the last change of a particular state variable (e.g., the time elapsed since wedding). 
In accordance therewith, \code{micSim} requires transition rates functions to feature three input parameters, namely \code{age}, \code{calTime}, and \code{duration}.
Via \code{age} the age of an individual is handed over, via \code{caltime} the calendar time, and via \code{duration} the time that has elapsed since the last change of the affected state variable.
All three input parameters might vary, or only one or two of them.
Also none of the input parameters can be specified to vary, i.e., a transition rate can be defined to be constant. 
If rates are assumed to be independent of a specific time scale, the corresponding input argument can simply be ignored within the body of the rates function (i.e., is not used to determine a specific rate value).
For illustration, see the examples in the example section. Note that allowing transition rates to vary along the time elapsed since a last transition facilitates modelling gestation gaps after a delivery: For a period of nine or ten months transition rates for higher order parities are simply set to zero (cf., the complex example in the example section).  

The transition matrix \code{transitionMatrix} has as many rows as the simulation model comprises nonabsorbing states and as many columns as the simulation model comprises absorbing and nonabsorbing states. The rows of \code{transitionMatrix} mark starting states of transitions and the columns mark arrival states. At positions of \code{transitionMatrix} indicating impossible transitions, the matrix contains zeros. Otherwise the name of the function determining the respective transition rates has to be given. The function \link{buildTransitionMatrix} supports the construction of \code{transitionMatrix}.

If, during simulation, an individual reaches \code{maxAge}, he/she stays in his/her current state until simulation ending date is reached, that is, the respective individual is no longer at risk of experiencing any events and his/her ongoing episode will be censored at simlation ending date.  

It is recommended to set \code{simHorizon} using the function \link{setSimHorizon}. 

Each element of \code{fertTr} has to be of the form "A->B", that is, "A" indicates the starting attribute of the transition and "B" the arrival attribute. ("->" is the placeholder defined to mark a transition.) For example, "0" (childless) gives the starting point of the transition marking a first birth event and "1" (first child) its arrival point. All fertility attributes given in \code{fertTr} have to be part of the state variable specifiying fertility in the state space. That is, if there is none, \code{fertTr} is empty: \code{fertTr=c()}.
}
\value{
 The data frame \code{pop} contains the whole synthetic population considered during simulation including all events generated. In more detail, \code{pop} contains as many rows as there are transitions performed by the individuals. (Also, "entering the population" is considered as an event. In general, individuals can enter the simulation via three channels: by being part of the starting population, by immigration, and by being born during simulation). 
 
The function \link{convertToLongFormat} reshapes the microsimulation output into long format. 
}
\note{
  Concerning run times \code{micSim} is not very performant. That is because it is purely implemented in R, i.e., it does not incorporate routines implemented in a high level programming language like Java, C++ or Python. This is work in progress. For the meantime, if a computer cluster is accessible, we recommend parallel computing using \link{micSimParallel}.
}
\author{
Sabine Zinn
}
\examples{
######################################################################################
# 1. Simple example only dealing with mortality events
######################################################################################
# Clean workspace 
rm(list=ls())

# Defining simulation horizon
simHorizon <- setSimHorizon(startDate="01/01/2000", endDate="31/12/2100")

# Seed for random number generator
set.seed(234)

# Definition of maximal age
maxAge <- 120

# Defintion of nonabsorbing and absorbing states
sex <- c("m","f")
stateSpace <- sex
attr(stateSpace,"name") <- "sex"
absStates <- "dead"

# Definition of an initial population 
dts <- c("31/12/1930","03/04/1999","15/10/1956","11/11/1991","01/01/1965")
birthDates <- chron(dates=dts,format=c(dates="d/m/Y"))
initStates <- c("f","m","f","m","m")
initPop <- data.frame(ID=1:5,birthDate=birthDates,initState=initStates)

# Definition of mortality rates (Gompertz model).
mortRates <- function(age, calTime, duration){
  a <- 0.00003
  b <- ifelse(calTime<=2020, 0.1, 0.097)
  rate <- a*exp(b*age)
  return(rate)
}

# Transition pattern and assignment of functions specifying transition rates
absTransitions <- c("dead","mortRates")
transitionMatrix <- buildTransitionMatrix(allTransitions=NULL,
  absTransitions=absTransitions, stateSpace=stateSpace)

# Execute microsimulation (sequentially, i.e., using only one CPU)
pop <- micSim(initPop=initPop, transitionMatrix=transitionMatrix, absStates=absStates, 
  maxAge=maxAge, simHorizon=simHorizon)
  
  
######################################################################################
# 2. More complex example dealing with mortality, changes in the fertily and the marital 
# status, in the educational attainment, as well as dealing with migration
######################################################################################
# Clean workspace 
rm(list=ls())

# Defining simulation horizon
simHorizon <- setSimHorizon(startDate="01/01/2014", endDate="31/12/2024")

# Seed for random number generator
set.seed(234)

# Definition of maximal age 
maxAge <- 100  

# Defintion of nonabsorbing and absorbing states
sex <- c("m","f")                     
fert <- c("0","1+")           
marital <- c("NM","M","D","W")        
edu <- c("no","low","med","high")   
stateSpace <- expand.grid(sex=sex,fert=fert,marital=marital,edu=edu)
absStates <- c("dead","rest")   

# General date of enrollment to elementary school
dateSchoolEnrol <- "09/01"

# Definition of an initial population (for illustration purposes, create a random population)
N = 100                                                       
initBirthDatesRange <- chron(dates=c("31/12/1950","01/01/2014"), format=c(dates="d/m/Y"), 
  out.format=c(dates="d/m/year"))
birthDates <- dates(initBirthDatesRange[1] + runif(N, min=0, max=diff(initBirthDatesRange)))
getRandInitState <- function(birthDate){
  age <- trunc(as.numeric(simHorizon[1] - birthDate)/365.25)
  s1 <- sample(sex,1)
  s2 <- ifelse(age<=18, fert[1], sample(fert,1))
  s3 <- ifelse(age<=18, marital[1], ifelse(age<=22, sample(marital[1:3],1), 
    sample(marital,1)))
  s4 <- ifelse(age<=7, edu[1], ifelse(age<=18, edu[2], ifelse(age<=23, sample(edu[2:3],1), 
    sample(edu[-1],1))))
  initState <- paste(c(s1,s2,s3,s4),collapse="/")
  return(initState)
}
initPop <- data.frame(ID=1:N, birthDate=birthDates, 
  initState=sapply(birthDates, getRandInitState))

# Definition of immigrants entering the population (for illustration purposes, create immigrants 
# randomly)
M = 20                                                           
immigrDatesRange <- as.numeric(simHorizon)
immigrDates <- dates(chron(immigrDatesRange[1] + runif(M, min=0,max=diff(immigrDatesRange)), 
  format=c(dates="d/m/Y", times="h:m:s"), out.format=c(dates="d/m/year",times="h:m:s")))
immigrAges <- runif(M, min=15*365.25, max=70*365.25)
immigrBirthDates <- dates(chron(as.numeric(immigrDates) - immigrAges, 
  format=c(dates="d/m/Y", times="h:m:s"), out.format=c(dates="d/m/year", times="h:m:s")))
IDmig <- max(as.numeric(initPop[,"ID"]))+(1:M)
immigrPop <- data.frame(ID = IDmig, immigrDate = immigrDates, birthDate=immigrBirthDates, 
  immigrInitState=sapply(immigrBirthDates, getRandInitState))  

# Definition of initial states for newborns 
initStates <- rbind(c("m","0","NM","no"),c("f","0","NM","no")) 
# Definition of related occurrence probabilities
initStatesProb <- c(0.515,0.485)                              

# Definition of (possible) transition rates  
# (1) Fertility rates (Hadwiger mixture model)
fert1Rates <- function(age, calTime, duration){  # parity 1
  b <- ifelse(calTime<=2020, 3.9, 3.3)
  c <- ifelse(calTime<=2020, 28, 29)
  rate <-  (b/c)*(c/age)^(3/2)*exp(-b^2*(c/age+age/c-2))
  rate[age<=15 | age>=45] <- 0
  return(rate)
}
fert2Rates <- function(age, calTime, duration){  # partiy 2+
  b <- ifelse(calTime<=2020, 3.2, 2.8)
  c <- ifelse(calTime<=2020, 32, 33)
  rate <-  (b/c)*(c/age)^(3/2)*exp(-b^2*(c/age+age/c-2))
  rate[age<=15 | age>=45 | duration<0.75] <- 0
  return(rate)
}
# (2) Rates for first marriage (normal density)
marriage1Rates <- function(age, calTime, duration){  
  m <- ifelse(calTime<=2020, 25, 30)
  s <- ifelse(calTime<=2020, 3, 3)
  rate <- dnorm(age, mean=m, sd=s)
  rate[age<=16] <- 0
  return(rate)
}
# (3) Remariage rates (log-logistic model)
marriage2Rates <- function(age, calTime, duration){  
  b <- ifelse(calTime<=2020, 0.07, 0.10)
  p <- ifelse(calTime<=2020, 2.7,2.7)
  lambda <- ifelse(calTime<=1950, 0.04, 0.03)
  rate <- b*p*(lambda*age)^(p-1)/(1+(lambda*age)^p)
  rate[age<=18] <- 0
  return(rate)
}
# (4) Divorce rates (normal density)
divorceRates <- function(age, calTime, duration){
  m <- 40
  s <- ifelse(calTime<=2020, 7, 6)
  rate <- dnorm(age,mean=m,sd=s)
  rate[age<=18] <- 0
  return(rate)
}
# (5) Widowhood rates (gamma cdf)
widowhoodRates <- function(age, calTime, duration){
  rate <- ifelse(age<=30, 0, pgamma(age-30, shape=6, rate=0.06))
  return(rate)
}
# (6) Rates to change educational attainment
# Set rate to `Inf' to make transition for age 7 deterministic.
noToLowEduRates <- function(age, calTime, duration){
  rate <- ifelse(age==7,Inf,0) 
  return(rate)
}
lowToMedEduRates <- function(age, calTime, duration){
  rate <- dnorm(age,mean=16,sd=1)
  rate[age<=15 | age>=25] <- 0
  return(rate)
}
medToHighEduRates <- function(age, calTime, duration){
  rate <- dnorm(age,mean=20,sd=3)
  rate[age<=18 | age>=35] <- 0
  return(rate)
}
# (7) Mortality rates (Gompertz model)
mortRates <- function(age, calTime, duration){
  a <- .00003
  b <- ifelse(calTime<=2020, 0.1, 0.097)
  rate <- a*exp(b*age)
  return(rate)
}
# (8) Emigration rates 
emigrRates <- function(age, calTime, duration){
  rate <- ifelse(age<=18,0,0.0025)
  return(rate)
}

# Transition pattern and assignment of functions specifying transition rates
fertTrMatrix <- cbind(c("0->1+","1+->1+"),                         
  c("fert1Rates", "fert2Rates"))
maritalTrMatrix <- cbind(c("NM->M","M->D","M->W","D->M","W->M"),              
  c("marriage1Rates","divorceRates","widowhoodRates",
 "marriage2Rates","marriage2Rates"))
eduTrMatrix <- cbind(c("no->low","low->med","med->high"),
  c("noToLowEduRates","lowToMedEduRates","medToHighEduRates")) 
allTransitions <- rbind(fertTrMatrix, maritalTrMatrix, eduTrMatrix)
absTransitions <- rbind(c("dead","mortRates"),c("rest","emigrRates"))
transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions,
  absTransitions=absTransitions, stateSpace=stateSpace)

# Define transitions triggering a birth event
fertTr <- fertTrMatrix[,1]

# Execute microsimulation (sequentially, i.e., using only one CPU core)
pop <- micSim(initPop=initPop, immigrPop=immigrPop, 
  transitionMatrix=transitionMatrix, absStates=absStates, 
  initStates=initStates, initStatesProb=initStatesProb, 
  maxAge=maxAge, simHorizon=simHorizon, fertTr=fertTr, 
  dateSchoolEnrol=dateSchoolEnrol)
}