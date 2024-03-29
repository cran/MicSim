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
micSim(initPop, immigrPop=NULL, transitionMatrix, absStates=NULL, 
        fixInitStates = c(), varInitStates=c(), initStatesProb=c(), 
        maxAge=99, simHorizon, fertTr=c(), monthSchoolEnrol=c())
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
  A matrix indicating the transition pattern and the names of the functions determining the respective transition rates (with rates to be returned as vectors, i.e. for input age 0 to 10 eleven rate values have to be returned).
}
  \item{absStates}{
 A vector indicating the absorbing states of the model.
}
  \item{fixInitStates}{
 (Vector of) Indices of substates determining the attributes/substates that a newborn will be taken over from the mother.
 If empty or not defined, no attributes will be inherited.
}
  \item{varInitStates}{
 (A vector comprising the) Substates / attributes that are assigned to a newborn randomly according to the probabilities \code{initStatesProb}, i.e. that are not inherited from the mother. 
}
  \item{initStatesProb}{
 A vector comprising the probabilities corresponding to \code{varInitStates}. 
 If \code{fixInitStates} are given (i.e. attributes from the mother are inherited), these probabilities have to sum to one   conditioned on the inherited attributes, i.e. for  each (set of) inherited attribute(s) separately. Otherwise, the sum of \code{initStatesProb} has to be one.
}
  \item{maxAge}{
 A scalar indicating the exact maximal age (i.e., sharp 100.00 years) which an individual can reach during simulation. \code{maxAge} has to be greater than zero. 
}
  \item{simHorizon}{
 A vector comprising the starting and ending date of the simulation. Both dates have to be given as strings in the format `yyyymmdd'. The starting date has to precede the ending date.
}
  \item{fertTr}{
A vector indicating all transitions triggering a child birth event during simulation, that is, the creation of a new individual. 
}
  \item{monthSchoolEnrol}{
The month (as numeric value from 1 to 12) indicating the general enrollment month for elementary school, e.g., 9 for September. If transition to elementary school is not defined (see below under `details') and no such month is given school enrollment to elementary school is not modelled / simulated. 
}

}
\details{
All nonabsorbing states considered during simulation have to be defined as composite states. In more detail, they consist of labels indicating values of state variables. Within states, labels are separated by a forward slash "/". Possible state variables are, for example, gender, number of children ever born, and educational attainment. Corresponding values are, for example, "m" and "f" (gender), "0","1","2", and "3+" (number of children ever born), "no", "low", "med", and "high" (educational attainment). Possible examples of states are "m/0/low" for a childless male with elementary education  or "f/1/high" for a female with one child and a higher secondary school degree. All state variables considered plus accordant value labels have to be provided by the user. The only exception is gender which is predefined by labels "m" and "f" indicating male and female individuals. The label values "no" and "low" are reserved for enrolment events to elementary school (see below).

Nonabsorbing states have to be given as strings such as "dead" for being dead or "rest" for emigrated.

\code{micSim} is able to conduct enrollment events to elementary school such that they take place on the first day of the \code{monthSchoolEnrol}th month of a particular year. For this purpose, a state variable defining educational attainment has to be created first. Then, labels of possible values have to be defined such that "no" describes no education and "low" describes elementary education. Finally, the transition function determining the transition rate for the respective enrollment event has to be defined to return "Inf" for the age x at which children should be enrolled (e.g., at age seven) and zero otherwise. That way, an event "school enrollment on \code{dateSchoolEnrol} of the year in which a child turns x years old" is enforced. A related illustration is given below in the second example.

If school enrollment is not of interest to the modeller, \code{monthSchoolEnrol} can let be unspecified. Then during simulation that feature is ignored. 

The starting population \code{initPop} has to be given in the form of a data frame. Each row of the data frame corresponds to one individual. \code{initPop} has to comprise the following information: 
unique numerical person identifier (ID), birth date, and initial state (i.e., the state occupied by the individual when entering the synthetic population). Birth dates have to be given as strings in the format `yyyymmdd', e.g. `20220815' for Aug 15th 2022. Be aware that at simulation starting date all individuals in the initial population have already to be born and younger than \code{maxAge}. Otherwise, \code{micSim} throws an error message pointing to this issue. 

Information about immigrants has to be given in the form of a data frame (\code{immigrPop}). Each row of the data frame corresponds to one immigrant.
\code{immigrPop} contains the following data: unique numerical person identifier (ID), immigration date, birth date, and initial state (i.e., the state occupied by the immigrant when entering the simulated population). Immigration dates and birth dates have to provided as strings in the format `yyyymmdd', e.g. `20220815' for Aug 15th 2022. Immigration dates have to be specified to occur after simulation starting date and before simulation stopping date. Immigrants must be born when they migrate. Otherwise, \code{micSim} throws error messages pointing to this issues.   

For each transition that should be considered during simulation accordant transition rates have to be provided. Since MicSim's model is a continuous-time multi-state model these rates are transition intensities (as also used for defining time-inhomogeneuous Markov models) and not probabilities. Palloni (2000) illustrates very well the difference between both concepts. Zinn (2011) describes methods for estimating rates for MicSim's model. A crude way of transforming transition probabilities to rates is assuming that the rates \eqn{lambda_{ij}} (for leaving state \eqn{i} to enter state \eqn{j}) are constant in the time interval (of length \eqn{t}) captured by a corresponding probability \eqn{p_{ij}}:

\eqn{p_{ij} = 1 - exp(- lambda_{ij} * t)} which yields \eqn{lambda_{ij} = - 1/t * ln(1-p_{ij})}.

Be aware that this is only an approximation since this formula belongs to a time-homogeneuous Markov model and not to the more flexible time-inhomogeneuous Markov model (as used by MicSim). Thus, here for the time interval covered by \eqn{p_{ij}} a time-homogeneuous Markov model is assumed. Many users may have annual transition probabilities at hand, i.e., \eqn{t=1}. 

\code{micSim} requires these rates in form of functions which are handed over via the transition matrix \code{transitionMatrix} (described in the subsequent paragraph). The \code{MicSim} package allows rates to depend on three time scales: age, calendar time, and the time that has elapsed since the last change of a particular state variable (e.g., the time elapsed since wedding). 
In accordance therewith, \code{micSim} requires transition rates functions to feature three input parameters, namely \code{age}, \code{calTime}, and \code{duration}.
Via \code{age} the age of an individual is handed over, via \code{caltime} the calendar time, and via \code{duration} the time that has elapsed since the last change of the affected state variable. 
All three input parameters might vary, or only one or two of them.
Also none of the input parameters can be specified to vary, i.e., transition rates can be defined to be constant. 
Since \code{micSim} computes integrals of rates along simulation procedure, the rates functions must deliver vector of rates for vectors of inputs, i.e. for an input vector of ages (e.g. ages [0,1,2,3]) the rates functions have to given as many rate values as is the length of the age vector (in the example, four rate values). More details on this are given in the examples below or in the vignette to this package. 
If rates are assumed to be independent of a specific time scale, the corresponding input argument can simply be ignored within the body of the rates function (i.e., is not used to determine a specific rate value).
For illustration, see the examples in the example section. 
Beware that rates for age have to be delivered at maximal only until \code{maxAge}. If *more* rates are given, this does not cause an error but they are not used.  
Note that allowing transition rates to vary along the time elapsed since a last transition facilitates modelling gestation gaps after a delivery: For a period of nine or ten months transition rates for higher order parities are simply set to zero (e.g., see the complex example in the example section).  

The transition matrix \code{transitionMatrix} has as many rows as the simulation model comprises nonabsorbing states and as many columns as the simulation model comprises absorbing and nonabsorbing states. The rows of \code{transitionMatrix} mark starting states of transitions and the columns mark arrival states. At positions of \code{transitionMatrix} indicating impossible transitions, the matrix contains zeros. Otherwise the name of the function determining the respective transition rates has to be given. The function \link{buildTransitionMatrix} supports the construction of \code{transitionMatrix}.

If, during simulation, an individual reaches \code{maxAge}, he/she stays in his/her current state until simulation ending date is reached, that is, the respective individual is no longer at risk of experiencing any events and his/her ongoing episode will be censored at simlation ending date.  

Each element of \code{fertTr} has to be of the form "A->B", that is, "A" indicates the starting attribute of the transition and "B" the arrival attribute. ("->" is the placeholder defined to mark a transition.) For example, "0" (childless) gives the starting point of the transition marking a first birth event and "1" (first child) its arrival point. All fertility attributes given in \code{fertTr} have to be part of the state variable specifiying fertility in the state space. That is, if there is none, \code{fertTr} is empty: \code{fertTr=c()}.
}
\value{
 The data frame \code{pop} contains the whole synthetic population considered during simulation including all events generated. In more detail, \code{pop} contains as many rows as there are transitions performed by the individuals. Also, "entering the population" is considered as an event. In general, individuals can enter the simulation via three channels: by being part of the starting population, by immigration, and by being born during simulation. If fertility events are part of the model's specification (i.e., \code{fertTr} is not empty), \code{pop} contains an additional column indicating the ID of the mother for individuals born during simulation. For all other individuals, the ID of the mother is unknown (i.e., set to `NA').
 
The function \link{convertToLongFormat} reshapes the microsimulation output into long format, while the function \link{convertToWideFormat} gives the microsimulation in wide format. 
}
\note{
For large-scale models and simulation, I recommend parallel computing using \link{micSimParallel}. This speeds up execution times considerably.
However, before running an extensive simulation on multiple cores, the package user should definitely check whether the input for the simulation fits. 
This can best be achieved by first running a short and less extensive simulation with only one core (e.g., running only a one percent sample of the initial population).}
\author{
Sabine Zinn
}

\references{
Palloni, A. (2001). Increment-Decrement Life Tables. In: Preston, S., Heuveline, P., & Guillot, M. (eds). Demography: measuring and modeling population processes. Malden, MA: Blackwell Publishers.

Zinn, S. (2011). Preparation of required input data. In: Zinn, S. A Continuous-Time Microsimulation and First Steps Towards a Multi-Level Approach in Demography, Disseration, Chapter 3, \url{https://rosdok.uni-rostock.de/file/rosdok_derivate_0000004766/Dissertation_Zinn_2011.pdf}

}

\examples{
######################################################################################
# 1. Simple example only dealing with mortality events
######################################################################################
\donttest{
# Clean workspace 
rm(list=ls())

# Defining simulation horizon
startDate <- 20000101 # yyyymmdd
endDate   <- 21001231 # yyyymmdd
simHorizon <- c(startDate=startDate, endDate=endDate)

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
birthDates <- c("19301231","19990403","19561015","19911111","19650101")
initStates <- c("f","m","f","m","m")
initPop <- data.frame(ID=1:5,birthDate=birthDates,initState=initStates)

# Definition of mortality rates (Gompertz model)
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
# 2. More complex, but only illustrative example dealing with mortality, changes in 
# fertily, and with the inheritance of attributes of the mother
######################################################################################  

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
nat <- c("DE","AT","IT") # nationality
fert <- c("0","1")           
stateSpace <- expand.grid(sex=sex,nat=nat,fert=fert)
absStates <- "dead" 

# Definition of an initial population (for illustration purposes, create a random population)
N = 100   
birthDates <- runif(N, min=getInDays(19500101), max=getInDays(20131231)) 
getRandInitState <- function(birthDate){
  age <- trunc((getInDays(simHorizon[1]) - birthDate)/365.25) 
  s1 <- sample(sex,1)
  s2 <- sample(nat,1)
  s3 <- ifelse(age<=18, fert[1], sample(fert,1))
  initState <- paste(c(s1,s2,s3),collapse="/")
  return(initState)
}
initPop <- data.frame(ID=1:N, birthDate=birthDates, initState=sapply(birthDates, getRandInitState))
initPop$birthDate <- getInDateFormat(initPop$birthDate)

# Definition of initial states for newborns
# To have possibility to define distinct sex ratios for distinct nationalities, 
# inherit related substate from the mother 
fixInitStates <- 2 # give indices for attribute/substate that will be taken over 
                   # from the mother, here: nat 
varInitStates <- rbind(c("m","DE","0"), c("f","DE","0"),
                       c("m","AT","0"), c("f","AT","0"), 
                       c("m","IT","0"), c("f","IT","0")) 
initStatesProb <- c(0.515,0.485, 
                    0.515,0.485,
                    0.515,0.485)
# Mind: depending on the inherited attribute nat="DE", nat="AT", or nat="IT"
# initials probabilites must sum to one                       

# Definition of (possible) transition rates  
# Fertility rates (Hadwiger mixture model)
fertRates <- function(age, calTime, duration){ 
  b <- ifelse(calTime<=2020, 3.5, 3.0)
  c <- ifelse(calTime<=2020, 28, 29)
  rate <-  (b/c)*(c/age)^(3/2)*exp(-b^2*(c/age+age/c-2))
  rate[age<=15 | age>=45] <- 0
  return(rate)
}
# Mortality rates (Gompertz model)
mortRates <- function(age, calTime, duration){
  a <- .00003
  b <- ifelse(calTime<=2020, 0.1, 0.097)
  rate <- a*exp(b*age)
  return(rate)
}

fertTrMatrix <- cbind(c("f/DE/0->f/DE/1", "f/AT/0->f/AT/1", "f/IT/0->f/IT/1"),
                      c(rep("fertRates",3))) 
allTransitions <- fertTrMatrix

absTransitions <- cbind(c("f/DE/dead", "f/AT/dead", "f/IT/dead", 
                          "m/DE/dead", "m/AT/dead", "m/IT/dead"),
                       c(rep("mortRates",6)))                          

transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions,
                                          absTransitions=absTransitions, 
                                          stateSpace=stateSpace)

# Define transitions triggering a birth event
fertTr <- fertTrMatrix[,1]

# Execute microsimulation 
pop <- micSim(initPop=initPop,
               transitionMatrix=transitionMatrix, absStates=absStates,
               varInitStates=varInitStates, initStatesProb=initStatesProb,
               fixInitStates=fixInitStates,
               maxAge=maxAge, simHorizon=simHorizon,fertTr=fertTr)  
  
######################################################################################
# 3. Complex example dealing with mortality, changes in the fertily and the marital 
# status, in the educational attainment, as well as dealing with migration
######################################################################################

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
N = 100                                                       
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
M = 20                                                           
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

# Execute microsimulation 
pop <- micSim(initPop=initPop, immigrPop=immigrPop, 
              transitionMatrix=transitionMatrix, 
              absStates=absStates, 
              varInitStates=varInitStates, 
              initStatesProb=initStatesProb, 
              maxAge=maxAge, 
              simHorizon=simHorizon, 
              fertTr=fertTr, 
              monthSchoolEnrol=monthSchoolEnrol)  
  
}
}