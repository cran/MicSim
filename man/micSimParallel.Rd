\name{micSimParallel}
\alias{micSimParallel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Run microsimulation (parallel computing)
}
\description{
The function \code{micSimParallel} is a parallelized version of the function \link{micSim}. That is, it runs a continuous-time microsimulation simulation distributed, i.e., using more than one CPU core.}
\usage{
micSimParallel(initPop, immigrPop = NULL, initPopList = c(), immigrPopList = c(), 
  transitionMatrix, absStates = NULL, varInitStates = c(), initStatesProb = c(), 
  fixInitStates = c(), maxAge = 99, simHorizon, fertTr = c(), 
  monthSchoolEnrol=c(), cores=1, seeds=1254)
}
\arguments{
 \item{initPop}{Either an initial population has to be given as a whole or splitted to be run at the distinct cores, see arguments \code{initPopList}. If it is given as a whole it is automatically splitted by MicSim such that the population parts run on the distinct cores are approx. equally sized.
}
 \item{immigrPop}{Optionally, a population of migrants entering the virtual population along simulation time can be given, see arguments \code{immigrPopList}. This migrants population can either be given as a whole or splitted to be run at the distinct cores. If it is given as a whole it is automatically splitted by MicSim such that the population parts run on the distinct cores are approx. equally sized.   
}
 \item{transitionMatrix}{
}
 \item{absStates}{
}
 \item{varInitStates}{
}
 \item{initStatesProb}{
}
 \item{initPopList}{
 Optional: A list containing the initial population split for the distinct cores. 
}
 \item{immigrPopList}{
 Optional: A list containing the immigration population split for the distinct cores. 
}
 \item{fixInitStates}{
}
 \item{maxAge}{
}
 \item{simHorizon}{
}
 \item{fertTr}{
}
 \item{monthSchoolEnrol}{See \link{micSim}.
}
 \item{cores}{
  Number of CPUs to be used.
}
  \item{seeds}{
  Seeds for pseudo number generators used for parallel computing.
}
}
\details{
The argument \code{cores} must not exceed the number of cores of the computer (cluster) used. 

In \code{seeds} as many seeds should be given as cores are used. If less are given, the latter are repeated to complete the set of seeds. 


}
\value{
 The data frame \code{pop} contains the whole synthetic population considered during simulation including all events generated. For more details, see \link{micSim}.
}
\author{
Sabine Zinn
}
\examples{
\donttest{

# Clean workspace 
rm(list=ls())

# Defining simulation horizon
startDate <- 20140101 # yyyymmdd
endDate   <- 20241231 # yyyymmdd
simHorizon <- c(startDate=startDate, endDate=endDate)

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

# General month of enrollment to elementary school
monthSchoolEnrol <- 9

# Definition of an initial population (for illustration purposes, create a random population)
N = 10000                                                       
birthDates <- runif(N, min=getInDays(19500101), max=getInDays(20131231)) 
getRandInitState <- function(birthDate){
  age <- trunc((getInDays(simHorizon[1]) - birthDate)/365.25) 
  s1 <- sample(sex,1)
  s2 <- ifelse(age<=18, fert[1], sample(fert,1))
  s3 <- ifelse(age<=18, marital[1], ifelse(age<=22, sample(marital[1:3],1), 
                                           sample(marital,1)))
  s4 <- ifelse(age<=7, edu[1], ifelse(age<=18, edu[2], ifelse(age<=23, sample(edu[2:3],1), 
                                                              sample(edu[-1],1))))
  initState <- paste(c(s1,s2,s3,s4),collapse="/")
  return(initState)
}
initPop <- data.frame(ID=1:N, birthDate=birthDates, initState=sapply(birthDates, getRandInitState))
initPop$birthDate <- getInDateFormat(initPop$birthDate)
range(initPop$birthDate)

# Definition of immigrants entering the population (for illustration purposes, create immigrants 
# randomly)
M = 2000                                                           
immigrDates <- runif(M, min=getInDays(20140101), max=getInDays(20241231)) 
immigrAges <- runif(M, min=15*365.25, max=70*365.25)
immigrBirthDates <- immigrDates - immigrAges
IDmig <- max(as.numeric(initPop[,"ID"]))+(1:M)
immigrPop <- data.frame(ID = IDmig, immigrDate = immigrDates, birthDate=immigrBirthDates, 
                        immigrInitState=sapply(immigrBirthDates, getRandInitState))  
immigrPop$birthDate <- getInDateFormat(immigrPop$birthDate)
immigrPop$immigrDate <- getInDateFormat(immigrPop$immigrDate)

# Definition of initial states for newborns 
varInitStates <- rbind(c("m","0","NM","no"),c("f","0","NM","no")) 
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

# Run microsimulation on cluster with three cores (settings depend on cluster used)
\dontrun{
cores <- 3
seeds <- c(1233,1245,265)
initPopList <- list(initPop[1:5000,], initPop[5001:8000,],initPop[8001:nrow(initPop),])
immigrPopList <- list(immigrPop[1:1000,], immigrPop[1001:1500,],immigrPop[1501:nrow(immigrPop),])

pop <- micSimParallel(initPopList=initPopList, immigrPopList=immigrPopList, 
  transitionMatrix=transitionMatrix, absStates=absStates, varInitStates=varInitStates, 
  initStatesProb=initStatesProb, maxAge=maxAge, simHorizon=simHorizon, 
  fertTr=fertTr, monthSchoolEnrol=monthSchoolEnrol, 
  cores=cores, seeds=seeds)
}
}
}