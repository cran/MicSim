#####################################################################################
#####################################################################################
### FUNCTION EXECUTING MICROSIMULATION                                             ##
### SZ, FEB 2023                                                                  ##
#####################################################################################
#####################################################################################
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# I. Execute microsimulation as single thread
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
#'
#' micSim
#'
#' Run microsimulation (sequentially)
#'
#' @usage micSim(initPop, immigrPop=NULL, transitionMatrix, absStates=NULL,
#'              fixInitStates = c(), varInitStates=c(), initStatesProb=c(),
#'              maxAge=99, simHorizon, fertTr=c(), monthSchoolEnrol=c())
#'
#' @description Performs a continuous-time microsimulation run (sequentially, i.e., using only one CPU core).
#'
#' @param initPop Data frame comprising the starting population of the simulation.
#' @param immigrPop Data frame comprising information about the immigrants entering the population across simulation time.
#' @param transitionMatrix A matrix indicating the transition pattern and the names of the functions determining the respective transition rates (with rates to be returned as vectors,
#' i.e. for input age 0 to 10 eleven rate values have to be returned).
#' @param absStates A vector indicating the absorbing states of the model.
#' @param fixInitStates (Vector of) Indices of subStates determining the attributes/subStates that a newborn will be taken over from the mother. If empty or not defined, no attributes will be inherited.
#' @param varInitStates (A vector comprising the) SubStates / attributes that are assigned to a newborn randomly according to the probabilities \code{initStatesProb}, i.e. that are not inherited from the mother.
#' @param initStatesProb A vector comprising the probabilities corresponding to \code{varInitStates}.
#' If \code{fixInitStates} are given (i.e. attributes from the mother are inherited), these probabilities have to sum to one conditioned on the inherited attributes,
#' i.e. for each (set of) inherited attribute(s) separately. Otherwise, the sum of \code{initStatesProb} has to be one.
#' @param maxAge A scalar indicating the exact maximal age (i.e., sharp 100.00 years) which an individual can reach during simulation. \code{maxAge} has to be greater than zero.
#' @param simHorizon A vector comprising the starting and ending date of the simulation. Both dates have to be given in the format 'yyyymmdd'. The starting date has to precede the ending date.
#' @param fertTr A vector indicating all transitions triggering a child birth event during simulation, that is, the creation of a new individual.
#' @param monthSchoolEnrol The month (as numeric value from 1 to 12) indicating the general enrollment month for elementary school, e.g., 9 for September.
#' If transition to elementary school is not defined (see below under 'details') and no such month is given school enrollment to elementary school is not modeled / simulated.
#'
#' @details All non-absorbing  states considered during simulation have to be defined as composite states. In more detail, they consist of labels indicating values of state variables.
#' Within states, labels are separated by a forward slash "/". Possible state variables are, for example, gender, number of children ever born, and educational attainment.
#' Corresponding values are, for example, "m" and "f" (gender), "0","1","2", and "3+" (number of children ever born), "no", "low", "med", and "high" (educational attainment).
#' Possible examples of states are "m/0/low" for a childless male with elementary education  or "f/1/high" for a female with one child and a higher secondary school degree.
#' All state variables considered plus accordant value labels have to be provided by the user. The only exception is gender which is predefined by labels "m" and "f" indicating male and female individuals.
#' The label values "no" and "low" are reserved for enrollment events to elementary school (see below).
#'
#' Non-absorbing states have to be given as strings such as "dead" for being dead or "rest" for emigrated.
#'
#' \code{micSim} is able to conduct enrollment events to elementary school such that they take place on the first day of the \code{monthSchoolEnrol}th month of a particular year.
#' For this purpose, a state variable defining educational attainment has to be created first.
#' Then, labels of possible values have to be defined such that "no" describes no education and "low" describes elementary education.
#' Finally, the transition function determining the transition rate for the respective enrollment event has to be defined to return "Inf" for the age x at which children should be enrolled (e.g., at age seven) and zero otherwise.
#' That way, an event "school enrollment on \code{dateSchoolEnrol} of the year in which a child turns x years old" is enforced. A related illustration is given below in the second example.
#'
#' If school enrollment is not of interest to the modeler, \code{monthSchoolEnrol} can let be unspecified. Then during simulation that feature is ignored.
#'
#' The starting population \code{initPop} has to be given in the form of a data frame. Each row of the data frame corresponds to one individual. \code{initPop} has to comprise the following information:
#' unique numerical person identifier (ID), birth date, and initial state (i.e., the state occupied by the individual when entering the synthetic population).
#' Birth dates have to be given as strings in the format 'yyyymmdd', e.g. '20220815' for Aug 15th 2022. Be aware that at simulation starting date all individuals in the initial population have already to be born and younger than \code{maxAge}.
#' Otherwise, \code{micSim} throws an error message pointing to this issue.
#'
#' Information about immigrants has to be given in the form of a data frame (\code{immigrPop}). Each row of the data frame corresponds to one immigrant.
#' \code{immigrPop} contains the following data: unique numerical person identifier (ID), immigration date, birth date, and initial state (i.e., the state occupied by the immigrant when entering the simulated population).
#' Immigration dates and birth dates have to provided as strings in the format 'yyyymmdd', e.g. '20220815' for Aug 15th 2022.
#' Immigration dates have to be specified to occur after simulation starting date and before simulation stopping date. Immigrants must be born when they migrate. Otherwise, \code{micSim} throws error messages pointing to this issues.
#'
#' For each transition that should be considered during simulation accordant transition rates have to be provided.
#' Since MicSim's model is a continuous-time multi-state model these rates are transition intensities (as also used for defining time-inhomogeneous Markov models) and not probabilities.
#' Palloni (2000) illustrates very well the difference between both concepts. Zinn (2011) describes methods for estimating rates for MicSim's model.
#' A crude way of transforming transition probabilities to rates is assuming that the rates \eqn{lambda_{ij}} (for leaving state \eqn{i} to enter state \eqn{j}) are constant in the time interval (of length \eqn{t}) captured by a corresponding probability \eqn{p_{ij}}:
#'
#' \eqn{p_{ij} = 1 - exp(- lambda_{ij} * t)} which yields \eqn{lambda_{ij} = - 1/t * ln(1-p_{ij})}.
#'
#' Be aware that this is only an approximation since this formula belongs to a time-homogeneous Markov model and not to the more flexible time-inhomogeneous Markov model (as used by MicSim).
#' Thus, here for the time interval covered by \eqn{p_{ij}} a time-homogeneous Markov model is assumed. Many users may have annual transition probabilities at hand, i.e., \eqn{t=1}.
#'
#' \code{micSim} requires these rates in form of functions which are handed over via the transition matrix \code{transitionMatrix} (described in the subsequent paragraph).
#' The \code{MicSim} package allows rates to depend on three time scales: age, calendar time, and the time that has elapsed since the last change of a particular state variable (e.g., the time elapsed since wedding).
#' In accordance therewith, \code{micSim} requires transition rates functions to feature three input parameters, namely \code{age}, \code{calTime}, and \code{duration}.
#' Via \code{age} the age of an individual is handed over, via \code{caltime} the calendar time, and via \code{duration} the time that has elapsed since the last change of the affected state variable, and via \code{duration} the time that has elapsed since the last change of the affected state variable.
#' All three input parameters might vary, or only one or two of them.
#' Also none of the input parameters can be specified to vary, i.e., transition rates can be defined to be constant.
#' Since \code{micSim} computes integrals of rates along simulation procedure, the rates functions must deliver vector of rates for vectors of inputs,
#' i.e. for an input vector of ages (e.g. ages (0,1,2,3)) the rates functions have to given as many rate values as is the length of the age vector (in the example, four rate values).
#' More details on this are given in the examples below or in the vignette to this package.
#' If rates are assumed to be independent of a specific time scale, the corresponding input argument can simply be ignored within the body of the rates function (i.e., is not used to determine a specific rate value).
#' For illustration, see the examples in the example section.
#' Beware that rates for age have to be delivered at maximal only until \code{maxAge}. If *more* rates are given, this does not cause an error but they are not used.
#' Note that allowing transition rates to vary along the time elapsed since a last transition facilitates modelling gestation gaps after a delivery:
#' For a period of nine or ten months transition rates for higher order parities are simply set to zero (e.g., see the complex example in the example section).
#'
#' The transition matrix \code{transitionMatrix} has as many rows as the simulation model comprises non absorbing states and as many columns as the simulation model comprises absorbing and non-absorbing states.
#' The rows of \code{transitionMatrix} mark starting states of transitions and the columns mark arrival states. At positions of \code{transitionMatrix} indicating impossible transitions, the matrix contains zeros.
#' Otherwise the name of the function determining the respective transition rates has to be given. The function \link{buildTransitionMatrix} supports the construction of \code{transitionMatrix}.
#'
#' If, during simulation, an individual reaches \code{maxAge}, he/she stays in his/her current state until simulation ending date is reached, that is,
#' the respective individual is no longer at risk of experiencing any events and his/her ongoing episode will be censored at simulation ending date.
#'
#' Each element of \code{fertTr} has to be of the form "A->B", that is, "A" indicates the starting attribute of the transition and "B" the arrival attribute. ("->" is the placeholder defined to mark a transition).
#' For example, "0" (childless) gives the starting point of the transition marking a first birth event and "1" (first child) its arrival point.
#' All fertility attributes given in \code{fertTr} have to be part of the state variable specifying fertility in the state space. That is, if there is none, \code{fertTr} is empty: \code{fertTr=c()}.
#'
#' NOTE: For large-scale models and simulation, I recommend parallel computing using \link{micSimParallel}. This speeds up execution times considerably.
#' However, before running an extensive simulation on multiple cores, the package user should definitely check whether the input for the simulation fits.
#' This can best be achieved by first running a short and less extensive simulation with only one core (e.g., running only a one percent sample of the initial population).
#'
#' @returns  The data frame \code{pop} contains the whole synthetic population considered during simulation including all events generated.
#' In more detail, \code{pop} contains as many rows as there are transitions performed by the individuals. Also, "entering the population" is considered as an event.
#' In general, individuals can enter the simulation via three channels: by being part of the starting population, by immigration, and by being born during simulation.
#' If fertility events are part of the model's specification (i.e., \code{fertTr} is not empty), \code{pop} contains an additional column indicating the ID of the mother for individuals born during simulation.
#' For all other individuals, the ID of the mother is unknown (i.e., set to 'NA').
#'
#' The function \link{convertToLongFormat} reshapes the microsimulation output into long format, while the function \link{convertToWideFormat} gives the microsimulation in wide format.
#'
#' @examples
#'
#' \dontrun{
#' ######################################################################################
#' # 1. Simple example only dealing with mortality events
#' ######################################################################################
#' # Clean workspace
#' rm(list=ls())
#'
#' # Defining simulation horizon
#' startDate <- 20000101 # yyyymmdd
#' endDate   <- 21001231 # yyyymmdd
#' simHorizon <- c(startDate=startDate, endDate=endDate)
#'
#' # Seed for random number generator
#' set.seed(234)
#'
#' # Definition of maximal age
#' maxAge <- 120
#'
#' # Definition of non-absorbing and absorbing states
#' sex <- c("m","f")
#' stateSpace <- sex
#' attr(stateSpace,"name") <- "sex"
#' absStates <- "dead"
#'
#' # Definition of an initial population
#' birthDates <- c("19301231","19990403","19561015","19911111","19650101")
#' initStates <- c("f","m","f","m","m")
#' initPop <- data.frame(ID=1:5,birthDate=birthDates,initState=initStates)
#'
#' # Definition of mortality rates (Gompertz model)
#' mortRates <- function(age, calTime){
#'   a <- 0.00003
#'   b <- ifelse(calTime<=2020, 0.1, 0.097)
#'   rate <- a*exp(b*age)
#'   return(rate)
#' }
#'
#' # Transition pattern and assignment of functions specifying transition rates
#' absTransitions <- c("dead","mortRates")
#' transitionMatrix <- buildTransitionMatrix(allTransitions=NULL,
#'   absTransitions=absTransitions, stateSpace=stateSpace)
#'
#' # Execute microsimulation (sequentially, i.e., using only one CPU)
#' pop <- micSim(initPop=initPop, transitionMatrix=transitionMatrix, absStates=absStates,
#'   maxAge=maxAge, simHorizon=simHorizon)
#'
#' ######################################################################################
#' # 2. More complex, but only illustrative example dealing with mortality, changes in
#' # fertility, and with the inheritance of attributes of the mother
#' ######################################################################################
#'
#' # Clean workspace
#' rm(list=ls())
#'
#' # Defining simulation horizon
#' startDate <- 20140101 # yyyymmdd
#' endDate   <- 20241231 # yyyymmdd
#' simHorizon <- c(startDate=startDate, endDate=endDate)
#'
#' # Seed for random number generator
#' set.seed(234)
#'
#' # Definition of maximal age
#' maxAge <- 100
#'
#' # Definition of non-absorbing and absorbing states
#' sex <- c("m","f")
#' nat <- c("DE","AT","IT") # nationality
#' fert <- c("0","1")
#' stateSpace <- expand.grid(sex=sex,nat=nat,fert=fert)
#' absStates <- "dead"
#'
#' # Definition of an initial population (for illustration purposes, create a random population)
#' N = 100
#' birthDates <- runif(N, min=getInDays(19500101), max=getInDays(20131231))
#' getRandInitState <- function(birthDate){
#'   age <- trunc((getInDays(simHorizon[1]) - birthDate)/365.25)
#'   s1 <- sample(sex,1)
#'   s2 <- sample(nat,1)
#'   s3 <- ifelse(age<=18, fert[1], sample(fert,1))
#'   initState <- paste(c(s1,s2,s3),collapse="/")
#'   return(initState)
#' }
#' initPop <- data.frame(ID=1:N, birthDate=birthDates, initState=sapply(birthDates, getRandInitState))
#' initPop$birthDate <- getInDateFormat(initPop$birthDate)
#'
#' # Definition of initial states for newborns
#' # To have possibility to define distinct sex ratios for distinct nationalities,
#' # inherit related subState from the mother
#' fixInitStates <- 2 # give indices for attribute/subState that will be taken over
#'                    # from the mother, here: nat
#' varInitStates <- rbind(c("m","DE","0"), c("f","DE","0"),
#'                        c("m","AT","0"), c("f","AT","0"),
#'                        c("m","IT","0"), c("f","IT","0"))
#' initStatesProb <- c(0.515,0.485,
#'                     0.515,0.485,
#'                     0.515,0.485)
#' # Mind: depending on the inherited attribute nat="DE", nat="AT", or nat="IT"
#' # initials probabilities must sum to one
#'
#' # Definition of (possible) transition rates
#' # Fertility rates (Hadwiger mixture model)
#' fertRates <- function(age, calTime){
#'   b <- ifelse(calTime<=2020, 3.5, 3.0)
#'   c <- ifelse(calTime<=2020, 28, 29)
#'   rate <-  (b/c)*(c/age)^(3/2)*exp(-b^2*(c/age+age/c-2))
#'   rate[age<=15 | age>=45] <- 0
#'   return(rate)
#' }
#' # Mortality rates (Gompertz model)
#' mortRates <- function(age, calTime){
#'   a <- .00003
#'   b <- ifelse(calTime<=2020, 0.1, 0.097)
#'   rate <- a*exp(b*age)
#'   return(rate)
#' }
#'
#' fertTrMatrix <- cbind(c("f/DE/0->f/DE/1", "f/AT/0->f/AT/1", "f/IT/0->f/IT/1"),
#'                       c(rep("fertRates",3)))
#' allTransitions <- fertTrMatrix
#'
#' absTransitions <- cbind(c("f/DE/dead", "f/AT/dead", "f/IT/dead",
#'                           "m/DE/dead", "m/AT/dead", "m/IT/dead"),
#'                        c(rep("mortRates",6)))
#'
#' transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions,
#'                                           absTransitions=absTransitions,
#'                                           stateSpace=stateSpace)
#'
#' # Define transitions triggering a birth event
#' fertTr <- fertTrMatrix[,1]
#'
#' # Execute microsimulation
#' pop <- micSim(initPop=initPop,
#'                transitionMatrix=transitionMatrix, absStates=absStates,
#'                varInitStates=varInitStates, initStatesProb=initStatesProb,
#'                fixInitStates=fixInitStates,
#'                maxAge=maxAge, simHorizon=simHorizon,fertTr=fertTr)
#'
#' ######################################################################################
#' # 3. Complex example dealing with mortality, changes in the fertility and the marital
#' # status, in the educational attainment, as well as dealing with migration
#' ######################################################################################
#'
#' # Clean workspace
#' rm(list=ls())
#'
#' # Defining simulation horizon
#' startDate <- 20140101 # yyyymmdd
#' endDate   <- 20241231 # yyyymmdd
#' simHorizon <- c(startDate=startDate, endDate=endDate)
#'
#' # Seed for random number generator
#' set.seed(234)
#'
#' # Definition of maximal age
#' maxAge <- 100
#'
#' # Definition of non-absorbing and absorbing states
#' sex <- c("m","f")
#' fert <- c("0","1+")
#' marital <- c("NM","M","D","W")
#' edu <- c("no","low","med","high")
#' stateSpace <- expand.grid(sex=sex,fert=fert,marital=marital,edu=edu)
#' absStates <- c("dead","rest")
#'
#' # General month of enrollment to elementary school
#' monthSchoolEnrol <- 9
#'
#' # Definition of an initial population (for illustration purposes, create a random population)
#' N = 100
#' birthDates <- runif(N, min=getInDays(19500101), max=getInDays(20131231))
#' getRandInitState <- function(birthDate){
#'   age <- trunc((getInDays(simHorizon[1]) - birthDate)/365.25)
#'   s1 <- sample(sex,1)
#'   s2 <- ifelse(age<=18, fert[1], sample(fert,1))
#'   s3 <- ifelse(age<=18, marital[1], ifelse(age<=22, sample(marital[1:3],1),
#'                                            sample(marital,1)))
#'   s4 <- ifelse(age<=7, edu[1], ifelse(age<=18, edu[2], ifelse(age<=23, sample(edu[2:3],1),
#'                                                               sample(edu[-1],1))))
#'   initState <- paste(c(s1,s2,s3,s4),collapse="/")
#'   return(initState)
#' }
#' initPop <- data.frame(ID=1:N, birthDate=birthDates, initState=sapply(birthDates, getRandInitState))
#' initPop$birthDate <- getInDateFormat(initPop$birthDate)
#' range(initPop$birthDate)
#'
#' # Definition of immigrants entering the population (for illustration purposes, create immigrants
#' # randomly)
#' M = 20
#' immigrDates <- runif(M, min=getInDays(20140101), max=getInDays(20241231))
#' immigrAges <- runif(M, min=15*365.25, max=70*365.25)
#' immigrBirthDates <- immigrDates - immigrAges
#' IDmig <- max(as.numeric(initPop[,"ID"]))+(1:M)
#' immigrPop <- data.frame(ID = IDmig, immigrDate = immigrDates, birthDate=immigrBirthDates,
#'                         immigrInitState=sapply(immigrBirthDates, getRandInitState))
#' immigrPop$birthDate <- getInDateFormat(immigrPop$birthDate)
#' immigrPop$immigrDate <- getInDateFormat(immigrPop$immigrDate)
#'
#' # Definition of initial states for newborns
#' varInitStates <- rbind(c("m","0","NM","no"),c("f","0","NM","no"))
#' # Definition of related occurrence probabilities
#' initStatesProb <- c(0.515,0.485)
#'
#' # Definition of (possible) transition rates
#' # (1) Fertility rates (Hadwiger mixture model)
#' fert1Rates <- function(age, calTime, duration){  # parity 1
#'   b <- ifelse(calTime<=2020, 3.9, 3.3)
#'   c <- ifelse(calTime<=2020, 28, 29)
#'   rate <-  (b/c)*(c/age)^(3/2)*exp(-b^2*(c/age+age/c-2))
#'   rate[age<=15 | age>=45 | duration<0.75] <- 0
#'   return(rate)
#' }
#' fert2Rates <- function(age, calTime){  # partiy 2+
#'   b <- ifelse(calTime<=2020, 3.2, 2.8)
#'   c <- ifelse(calTime<=2020, 32, 33)
#'   rate <-  (b/c)*(c/age)^(3/2)*exp(-b^2*(c/age+age/c-2))
#'   rate[age<=15 | age>=45] <- 0
#'   return(rate)
#' }
#' # (2) Rates for first marriage (normal density)
#' marriage1Rates <- function(age, calTime){
#'   m <- ifelse(calTime<=2020, 25, 30)
#'   s <- ifelse(calTime<=2020, 3, 3)
#'   rate <- dnorm(age, mean=m, sd=s)
#'   rate[age<=16] <- 0
#'   return(rate)
#' }
#' # (3) Remariage rates (log-logistic model)
#' marriage2Rates <- function(age, calTime){
#'   b <- ifelse(calTime<=2020, 0.07, 0.10)
#'   p <- ifelse(calTime<=2020, 2.7,2.7)
#'   lambda <- ifelse(calTime<=1950, 0.04, 0.03)
#'   rate <- b*p*(lambda*age)^(p-1)/(1+(lambda*age)^p)
#'   rate[age<=18] <- 0
#'   return(rate)
#' }
#' # (4) Divorce rates (normal density)
#' divorceRates <- function(age, calTime){
#'   m <- 40
#'   s <- ifelse(calTime<=2020, 7, 6)
#'   rate <- dnorm(age,mean=m,sd=s)
#'   rate[age<=18] <- 0
#'   return(rate)
#' }
#' # (5) Widowhood rates (gamma cdf)
#' widowhoodRates <- function(age, calTime){
#'   rate <- ifelse(age<=30, 0, pgamma(age-30, shape=6, rate=0.06))
#'   return(rate)
#' }
#' # (6) Rates to change educational attainment
#' # Set rate to 'Inf' to make transition for age 7 deterministic.
#' noToLowEduRates <- function(age, calTime){
#'   rate <- ifelse(age==7,Inf,0)
#'   return(rate)
#' }
#' lowToMedEduRates <- function(age, calTime){
#'   rate <- dnorm(age,mean=16,sd=1)
#'   rate[age<=15 | age>=25] <- 0
#'   return(rate)
#' }
#' medToHighEduRates <- function(age, calTime){
#'   rate <- dnorm(age,mean=20,sd=3)
#'   rate[age<=18 | age>=35] <- 0
#'   return(rate)
#' }
#' # (7) Mortality rates (Gompertz model)
#' mortRates <- function(age, calTime){
#'   a <- .00003
#'   b <- ifelse(calTime<=2020, 0.1, 0.097)
#'   rate <- a*exp(b*age)
#'   return(rate)
#' }
#' # (8) Emigration rates
#' emigrRates <- function(age, calTime){
#'   rate <- ifelse(age<=18,0,0.0025)
#'   return(rate)
#' }
#'
#' # Transition pattern and assignment of functions specifying transition rates
#' fertTrMatrix <- cbind(c("0->1+","1+->1+"),
#'   c("fert1Rates", "fert2Rates"))
#' maritalTrMatrix <- cbind(c("NM->M","M->D","M->W","D->M","W->M"),
#'   c("marriage1Rates","divorceRates","widowhoodRates",
#'  "marriage2Rates","marriage2Rates"))
#' eduTrMatrix <- cbind(c("no->low","low->med","med->high"),
#'   c("noToLowEduRates","lowToMedEduRates","medToHighEduRates"))
#' allTransitions <- rbind(fertTrMatrix, maritalTrMatrix, eduTrMatrix)
#' absTransitions <- rbind(c("dead","mortRates"),c("rest","emigrRates"))
#' transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions,
#'   absTransitions=absTransitions, stateSpace=stateSpace)
#'
#' # Define transitions triggering a birth event
#' fertTr <- fertTrMatrix[,1]
#'
#' # Execute microsimulation
#' pop <- micSim(initPop=initPop, immigrPop=immigrPop,
#'               transitionMatrix=transitionMatrix,
#'               absStates=absStates,
#'               varInitStates=varInitStates,
#'               initStatesProb=initStatesProb,
#'               maxAge=maxAge,
#'               simHorizon=simHorizon,
#'               fertTr=fertTr,
#'               monthSchoolEnrol=monthSchoolEnrol)
#'
#' }
#' @export
#'
micSim <- function(initPop, immigrPop=NULL, transitionMatrix, absStates=NULL, fixInitStates = c(),
                   varInitStates=c(), initStatesProb=c(), maxAge=99, simHorizon, fertTr=c(), monthSchoolEnrol=c()) {
  
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # A. CHECK INPUT FOR CONSISTENCY
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  if(is.null(initPop))
    stop('No starting population has been defined.')
  if(!is.null(initPop)){
    if(paste(colnames(initPop),collapse='/')!='ID/birthDate/initState')
      stop('Matrix specifying the starting population has not been defined properly.')
  }
  if(!is.null(immigrPop)){
    if(paste(colnames(immigrPop),collapse='/')!='ID/immigrDate/birthDate/immigrInitState')
      stop('Matrix specifying immigrants has not been defined properly.')
  }
  if(is.null(transitionMatrix))
    stop('Matrix defining transition pattern und functions has not been defined properly.')
  if(maxAge<=0)
    stop('The maximal age until which individual life courses are simulated should exceed zero.')
  if(length(simHorizon)!=2)
    stop('The simulation horizon has not been defined properly.')
  if(is.null(absStates))
    absStates <- setdiff(colnames(transitionMatrix),rownames(transitionMatrix))
  if(length(fertTr)>0){
    if((is.null(varInitStates) & is.null(initStatesProb)))
      stop('For children potentially born during simulation no inital state(s) and/or corresponding occurrence probabilities have been defined.')
    if(length(fixInitStates)>0){
      for(i in 1:length(fixInitStates)){
        ssum <- sum(initStatesProb[apply(varInitStates,1, function(rr){varInitStates[,fixInitStates[i]][1] %in% rr})])
        if(ssum!=1)
          stop('The sum of the probabilities to assign initial states to newborns must equal 1.')
      }
    } else {
      if(sum(initStatesProb)!=1)
        stop('The sum of the probabilities to assign initial states to newborns must equal 1.')
    }
  }
  # check whether all functions delivering transition rates deliver vectors of rates as output (necessary for integration procedure later on)
  allTr <- unique(as.vector(transitionMatrix)[as.vector(transitionMatrix) !="0"])
  simStartInDays <- getInDays(simHorizon[1])
  simStopInDays  <- getInDays(simHorizon[2])
  ranYear <- c(getYear(simStartInDays), getYear(simStopInDays))
  if(length(fertTr)>0){
    minAge <- 0
  } else {
    minAge <- min(trunc(getAgeInDays(simHorizon[1], initPop$birthDate)/365.25))
  }
  ranAge <- c(minAge,maxAge)
  ran <- min(c(diff(ranYear), diff(ranAge)))
  depMatrix <- rate_cS(allTr)
  for (tr in 1:length(depMatrix)) {
    if(ncol(depMatrix) >= 3){
      tr_dur <- rownames(depMatrix)[depMatrix[,3] == 1]
      if(length(tr_dur) != 0){
        for(tr in 1:length(tr_dur)){ # check whether all functions delivering transition rates deliver vectors of rates as output (necessary for integration procedure later on)
          for(cal in c(getYear(simStartInDays): getYear(simStopInDays))){
            for(age in minAge:(maxAge-1)){
              for(dur in 0:ran){
                res <- eval(do.call(tr_dur[tr], args=list(age=age,calTime=cal,duration=dur)))
                if(anyNA(res)){
                  cat("The rates function for ", allTr[i], " does not deliver a vector of rates for an input vector of age, calendar time, and/or duration (all in years).\n")
                  cat("The missing rate occurs at year ",cal, " for age ", age, " and duration ", dur, "\n.")
                  cat("This is a requirement for the simulation procedure to run since it is based on integrated hazard rates.\n")
                  stop('Incorrect definition of input rates function!')
                }
              }
            }
          }
        }
      }
      else{
        tr_nd <- rownames(depMatrix)[depMatrix[,3] != 1]
        for(tr in 1:length(tr_nd)){ # check whether all functions delivering transition rates deliver vectors of rates as output (necessary for integration procedure later on)
          for(cal in c(getYear(simStartInDays): getYear(simStopInDays))){
            for(age in minAge:(maxAge-1)){
              res <- eval(do.call(tr_nd[tr], args=list(age=age,calTime=cal)))
              if(anyNA(res)){
                cat("The rates function for ", allTr[i], " does not deliver a vector of rates for an input vector of age and/ or calendar time (all in years).\n")
                cat("The missing rate occurs at year ",cal, " and for age ", age, "\n.")
                cat("This is a requirement for the simulation procedure to run since it is based on integrated hazard rates.\n")
                stop('Incorrect definition of input rates function!')
              }
            }
          }
        }
      }
    }
    if (ncol(depMatrix) == 2) {
      tr_nd <- rownames(depMatrix)
      for (tr in 1:length(tr_nd)) {
        for (cal in getYear(simStartInDays):getYear(simStopInDays)) {
          for (age in minAge:(maxAge - 1)) {
            res <- eval(do.call(tr_nd[tr], args = list(age = age, calTime = cal)))
            if (anyNA(res)) {
              cat("The rates function for ", tr_nd[tr], " does not deliver a vector of rates for an input vector of age and calendar time.\n")
              cat("The missing rate occurs at year ", cal, " for age ", age, "\n.")
              cat("This is a requirement for the simulation procedure to run since it is based on integrated hazard rates.\n")
              stop('Incorrect definition of input rates function!')
            }
          }
        }
      }
    }
  }
  if(length(monthSchoolEnrol)==0){
    schoolEnrol <- FALSE
  } else {
    schoolEnrol <- TRUE
  }
  
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # B. DEFINITION OF GLOBAL PARAMETERS
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  
  # Simulation horizon
  simStartInDays <- getInDays(simHorizon[1])
  simStopInDays  <- getInDays(simHorizon[2])
  # Assign to each string state and substate an unique numerical code
  codingScheme <- builtStatesCodes(transitionMatrix)
  nSubStates <- ceiling(max(floor(log10(abs(codingScheme[,2]))) + 1)/2) # number of substates (without absorbing states)
  # Put numerical codes to transition matrix as well
  transitionMatrixNum <- transitionMatrix
  indCodesR <- match(rownames(transitionMatrix), codingScheme[,1])
  rownames(transitionMatrixNum) <- as.numeric(codingScheme[indCodesR,2])
  indCodesC <- match(colnames(transitionMatrix), codingScheme[,1])
  colnames(transitionMatrixNum) <- as.numeric(codingScheme[indCodesC,2])
  absStatesNum <- -c(1:length(absStates))
  if(!is.null(varInitStates)){
    varInitStatesStr <- apply(varInitStates, 1, paste, collapse="/")
    varInitStatesNum <- codingScheme[codingScheme[,1] %in% varInitStatesStr, c(2,2+c(1:nSubStates))]
  }
  # Event queue
  queue <- matrix(NA,ncol=9,nrow=0) # columns: 'ID','currTime','currState','currAge','nextState','timeToNextState', 'birthtime', 'initState', "isMig"
  # Global time
  t.clock <- simStartInDays  # counts in days since 01-01-1970
  # Recording transitions performed
  transitions <- matrix(NA,ncol=5,nrow=0) # columns: ID, From, To, transitionTime, transitionAge
  # Record linkage of mothers to newborns via their IDs
  if(length(fertTr)>0){
    mothers <- matrix(NA,ncol=2,nrow=0) # 'motherID', 'childID'
  }
  # Maximal Id / counter for individuals in the simulation
  maxId <- 0
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # C. FUNCTIONS REQUIRED FOR SIMULATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  
  # Function building matrix indicating the transitions between states causing a newborn
  buildFertTrExpanded <- function(){
    allStates <- matrix(rownames(transitionMatrix), ncol=1)
    allStatesSplit <- apply(allStates,1,strsplit, split="/")
    fert <- do.call(rbind,(strsplit(fertTr,split='->')))
    fertTrExpandedStr <- NULL
    for(i in 1:nrow(allStates)){
      cS <- allStatesSplit[[i]][[1]]
      for(j in 1:nrow(allStates)){
        dS <- allStatesSplit[[j]][[1]]
        if(("f" %in% cS) & ("f" %in% dS)){
          for(k in 1:nrow(fert)){
            ff <- fert[k,]
            oS <- strsplit(ff[1],'/')[[1]]
            bS <- strsplit(ff[2],'/')[[1]]
            cond1 <- !(F %in% (oS %in% cS)) & !(F %in% (bS %in% dS))
            cond2 <- paste((cS[!(cS %in% oS)]),collapse="/") == paste((dS[!(dS %in% bS)]),collapse="/") # if there are a fertility event only one substate can change, namely that one belonging to the fertility attribute
            if(cond1 & cond2){
              fertTrExpandedStr <- rbind(fertTrExpandedStr, c(paste0(cS,collapse="/"), paste0(dS,collapse="/")))
            }
          }
        }
      }
    }
    indCodes1 <- match(fertTrExpandedStr[,1], codingScheme[,1]) # transform numerical state codes to strings according to codingScheme
    indCodes2 <- match(fertTrExpandedStr[,2], codingScheme[,1])
    fertTrExpanded <- cbind(codingScheme[indCodes1,2],codingScheme[indCodes2,2])
    return(fertTrExpanded)
  }
  
  # Function checks whether a transition causes a newborn.
  # (Demands 'fertTrExpanded': matrix indicating the transitions between states causing a newborn (defined by 'fertTr').)
  isBirthEvent <- function(currState, destState){
    oS <- which(fertTrExpanded[,1] %in% currState)
    dS <- which(fertTrExpanded[,2] %in% destState)
    if(length(intersect(oS,dS))>0)
      return(TRUE)
    return(FALSE)
  }
  
  # Function adds to simulation population a newborn (using combined sim. step with duration dep. trans. functions)
  addNewNewborn_dur <- function(birthTime=birthTime, motherID=motherID, motherState=motherState){
    if(length(fixInitStates)>0){
      motherStatePart <-codingScheme[codingScheme[,2] %in% motherState, 2+c(1:nSubStates)][fixInitStates]
      if(length(fixInitStates)>1){
        notSuitedStatesInd <- which(is.na(apply(varInitStatesNum[,-1,drop=F][,fixInitStates,drop=F],1,match, x=motherStatePart)), arr.ind=TRUE)[,2]
      } else{
        notSuitedStatesInd <- which(is.na(apply(varInitStatesNum[,-1,drop=F][,fixInitStates,drop=F],1,match, x=motherStatePart)), arr.ind=TRUE)
      }
      inSt <- setdiff(1:nrow(varInitStatesNum), notSuitedStatesInd)
      varInitStatesR <- varInitStatesNum[inSt,,drop=F]
      initStatesProbR <- initStatesProb[inSt]
      birthState <- varInitStatesR[sample(1:nrow(varInitStatesR),size=1,replace=T,prob=initStatesProbR),1]
    } else {
      birthState <- varInitStatesNum[sample(1:nrow(varInitStatesNum),size=1,replace=T,prob=initStatesProb),1]
    }
    maxId <<- maxId + 1
    newInd <- c(maxId,getInDateFormat(birthTime),codingScheme[codingScheme[,2] %in% birthState,1])
    #cat('NewBorn: ',newInd,'\n')
    initPop <<- rbind(initPop,newInd)
    mothers <<- rbind(mothers, c(motherID, maxId))
    nE <- getNextStep(c(maxId,birthState,0,birthTime,motherID))
    #cat('\n------------n')
  }
  
  # Function building matrix indicating the transitions between states causing a school enrollment
  buildEduTrExpanded <- function(){
    allStates <- matrix(rownames(transitionMatrix), ncol=1)
    eduTrExpandedStr <- NULL
    for(i in 1:nrow(allStates)){
      for(j in 1:nrow(allStates)){
        cond <- all( c(any('no' %in% strsplit(allStates[i],'/')[[1]]), any('low' %in% strsplit(allStates[j],'/')[[1]])))
        if(cond){
          eduTrExpandedStr <- rbind(eduTrExpandedStr,c(allStates[i], allStates[j]))
        }
      }
    }
    eduInd <- gsub(x=eduTrExpandedStr[,1], pattern="no", replacement="")==gsub(x=eduTrExpandedStr[,2], pattern="low", replacement="")
    eduTrExpandedStr <- eduTrExpandedStr[eduInd,]
    indCodes1 <- match(eduTrExpandedStr[,1], codingScheme[,1]) # transform numerical state codes to strings according to codingScheme
    indCodes2 <- match(eduTrExpandedStr[,2], codingScheme[,1])
    eduTrExpanded <- cbind(codingScheme[indCodes1,2],codingScheme[indCodes2,2])
    return(eduTrExpanded)
  }
  
  # Function checks whether a transition implies a school enrollment (in the year when child turns seven).
  # (If state 1 comprises value 'no' and state 2 comprises value 'low', the transition is marked as 'school enrollment'.)
  isSchoolEnrolment <- function(currState,destState){
    oS <- which(eduTrExpanded[,1] %in% currState)
    dS <- which(eduTrExpanded[,2] %in% destState)
    if(length(intersect(oS,dS))>0)
      return(TRUE)
    return(FALSE)
  }
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # D. SIMULATION STEP
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  
  # Function to compute the next transition state and time of an individual who at age 'currAge', at time 'calTime'
  # entered its current state 'currState'. At this, consider possible duration dependencies of transition rates.
  getNextStep <- function(inp, isIMInitEvent=F){
    
    # Extract input data
    id <- inp[1]
    currState <- inp[2] # current state in numerical code
    currAge <- inp[3] # age in days
    calTime <- inp[4] # calendar time in days since 01-01-1970
    birthTime <- inp[5] # birth time in days since 01-01-1970
    initState <- inp[6] # initial state at sim. start, NA for individuals not yet born at that time
    isMig <- inp[7] # is migrant
    
    # First event of an immigrant: he/she enters the population later than sim. starting time
    lagToWaitingTime <- ifelse(isIMInitEvent, (calTime - simStartInDays)/365.25,0) # in years
    #cat('\n-----\nID: ',id,'\n')
    #print(inp)
    ageInYears <- currAge/365.25
    #cat('Age: ',ageInYears,' - CalTime: ',getYear(calTime),'-',getMonth(calTime),'-',getDay(calTime),'\n')
    # Possible destination states
    possTr <- transitionMatrixNum[match(currState, rownames(transitionMatrixNum)),]
    possTr <- possTr[which(possTr !=0)]
    nextEventMatrix <- matrix(0, ncol=2, nrow=length(possTr))
    
    # How many years (along age scale) remain until 'maxAge'?
    ranMaxAge <- (maxAge-0.01)-ageInYears
    # How many years (along cal. time scale) remain until simulation end?
    ranMaxYear <-  (simStopInDays - calTime)/365.25
    # Identify the time range that should be considered.
    ran <- min(ranMaxYear,ranMaxAge)
    #cat('Ran: ',ran,' - ranMaxAge: ',ranMaxAge,' - ranMaxYear: ',ranMaxYear,'\n')
    #ranAge <- c(ageInYears,ageInYears+ranMaxAge) # age range in years
    #ranYear <- c(getYear(calTime), getYear(calTime)+ran) # year range in years
    #cat('RanAge: ',ranAge,' - ranYear: ',ranYear,'\n')
    
    if(ncol(depMatrix) >= 3){
      tr_dur <- rownames(depMatrix)[depMatrix[,3] == 1]
      # Extract transition history of individual until current cal. time.
      historiesInd <- transitions[transitions[,1] %in% id & transitions[,4] <= calTime,,drop=F]
      # Extract for each state the duration until transition (in days).
      # Here, we have to differ between states of which we do not know when they are entered (i.e., 'initial states' of members
      # of the starting population and the states of migrants when they entered the country), and
      # the states we know the 'entering date' as well as the 'leaving date' (if the state has been left).
      if(birthTime < simStartInDays | isMig %in% 1) {
        dur <- rbind(c(initState,NA),historiesInd[,c(3,4),drop=F]) # extract from historiesInd 'To' and 'transitionTime'
        dur <- cbind(dur,c(diff(dur[,2]),0)) # columns: TransitionTo, AtTime, durUntil
        dur[which(is.na(dur[,2])),3] <- NA
      } else {  # Individual is born during simulation.
        dur <- rbind(c(initState,birthTime),historiesInd[,c(3,4),drop=F])
        dur <- cbind(dur,c(diff(dur[,2]),0)) # columns: TransitionTo, AtTime, durUntil
      }
      # Compute for each possible destination state a waiting time.
      for(i in 1:length(possTr)){
        tr <- possTr[i]
        destState <-  as.numeric(names(tr))
        cS <- codingScheme[codingScheme[,2] %in% currState, 2+c(1:nSubStates)]
        dS <- codingScheme[codingScheme[,2] %in% destState, 2+c(1:nSubStates)]
        # To determine the duration (time elapsed since last transition) that applies for the considered destination state,
        # we have to determine the duration since the last change in the covariate affected.
        # For example, to specify the time being married, we have to determine the duration since (last) marriage.
        covToCh <- which((cS==dS)==F)
        durSinceLastCovCh <- Inf  # For the transition to 'dead' so far the time elapsed since the last transition does not play any role.
        if(length(covToCh)==1){
          covHist <- codingScheme[codingScheme[,2] %in% dur[,1], 2+c(1:nSubStates),drop=F][,covToCh]
          idd <- which(covHist==cS[covToCh])
          if(length(idd)>1){
            if(F %in% (diff(idd)==1)){
              y <- rev(idd)[c(-1,diff(rev(idd)))==-1]
              idd <- rev(y)[c(diff(rev(y)),1)==1]
            }
          }
          durSinceLastCovCh <- sum(dur[idd,3]) # If I do not know how long an individual already is in a state: This gives NA.
          if(is.na(durSinceLastCovCh))
            durSinceLastCovCh <- currAge # Then assume the individual is already for his/her whole life in the state.
        }
        if(length(covToCh)>1 & (!destState %in% absStatesNum)){
          cat('Recognized a possible transition implying a change of two or more covariates.',
              'Concerning the derivation of the time being elapsed since the last transition this feature is not yet implemented.',
              'Current State: ',currState,' -> Possible transition to ',destState,'\n')
        }
        tageInYears <- trunc(ageInYears)
        tCalTime <- trunc(1970.001+calTime/365.25)
        tdurSinceLastCovCh <- trunc(durSinceLastCovCh/365.25)
        if(tr %in% tr_dur) {
          indRateFctDET <- function(x){
            res <- eval(do.call(tr,
                                args=list(age=tageInYears+x,calTime=tCalTime+x,duration=tdurSinceLastCovCh+x)))
            return(res)
          }
          ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
          detE <- indRateFctDET(ranAccuracyInDays)
          daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
          if (Inf %in% detE) {
            timeToNext <- daysToTrInYears
          } else {
            u <- -log(1-runif(1))
            #cat('It: ',i,'--u: ',u,'\n')
            # Extract individual transition rate (depending on age, calendar time, and time elapsed)
            indRateFct <- function(x){
              ageIn <- ageInYears+x
              calIn <- 1970.001+calTime/365.25+x
              durIn <- durSinceLastCovCh/365.25+x
              res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn,duration=durIn)))
              if(TRUE %in% (res<0))
                stop('I have found negative rate value/s for transition: ',tr,'\n
                   This is implausible. Please check this. Simulation has been stopped.\n')
              #cat('x: ',x,' -- res', res,'\n')
              #cat('\n---\n')
              return(res)
            }
            
            if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
              intHaz <- 0
            } else {
              # Integrated hazard at max. value
              intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
              if(inherits(intHaz, 'try-error')){
                intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
              }
            }
            # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
            if(u<=intHaz){
              invHazFct <- function(x){
                #cat('x: ',x,'\n')
                try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
                #print(try.res)
                if(inherits(try.res, 'try-error')){
                  #cat('Seemingly, divergent intergral for ID ',id,
                  # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                  #  ' for random number: ',u,'\n')
                  try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
                }
                #cat('res: ',try.res,'\n-----\n')
                return(try.res)
              }
              # Find random waiting time.
              timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
            } else {
              timeToNext <- Inf
            }
          } }else {
            indRateFctDET <- function(x){
              res <- eval(do.call(tr,
                                  args=list(age=tageInYears+x,calTime=tCalTime+x)))
              return(res)
            }
            ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
            detE <- indRateFctDET(ranAccuracyInDays)
            daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
            if (Inf %in% detE) {
              timeToNext <- daysToTrInYears
            } else {
              u <- -log(1-runif(1))
              #cat('It: ',i,'--u: ',u,'\n')
              # Extract individual transition rate (depending on age and calendar time)
              indRateFct <- function(x){
                ageIn <- ageInYears+x
                calIn <- 1970.001+calTime/365.25+x
                res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn)))
                if(TRUE %in% (res<0))
                  stop('I have found negative rate value/s for transition: ',tr,'\n
                   This is implausible. Please check this. Simulation has been stopped.\n')
                #cat('x: ',x,' -- res', res,'\n')
                #cat('\n---\n')
                return(res)
              }
              if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
                intHaz <- 0
              } else {
                # Integrated hazard at max. value
                intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
                if(inherits(intHaz, 'try-error')){
                  intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
                }
              }
              # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
              if(u<=intHaz){
                invHazFct <- function(x){
                  #cat('x: ',x,'\n')
                  try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
                  #print(try.res)
                  if(inherits(try.res, 'try-error')){
                    #cat('Seemingly, divergent intergral for ID ',id,
                    # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                    #  ' for random number: ',u,'\n')
                    try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
                  }
                  #cat('res: ',try.res,'\n-----\n')
                  return(try.res)
                }
                # Find random waiting time.
                timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
              } else {
                timeToNext <- Inf
              }
            }
          }
        nextEventMatrix[i,1] <- destState
        nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
      }
      #print(nextEventMatrix)
      nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
      if(dim(nE)[1]>1)
        nE <- nE[1,,drop=F]
      if(nE[1,2]!=Inf){
        # Cal. time of next event of individual. (If there is one.)
        tt <- calTime + as.numeric(nE[1,2])
        #print(tt)
        #cat(nE[1,1],'---',tt,'\n')
        # Check whether next event implies school enrollment. If yes, adjust transition time to ensure that the individual
        # enters school at Sept. 1 in the year he/she turns seven.
        if(schoolEnrol){ # Is school enrollment considered in this simulation model? Yes, then continue; otherwise skip this part
          if(isSchoolEnrolment(currState,nE[1,1])){
            enYear <- getYear(tt)
            if(getMonth(tt) <= monthSchoolEnrol) {
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            } else {
              enYear <- enYear+1
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            }
            diffToEn <- as.numeric(enDate-tt)
            nE[1,2] <- as.numeric(nE[1,2]) + diffToEn
          }
        }
        # Enqueue new event (if there is one).
        queue <<- rbind(queue, c(id, t.clock, currState, currAge - lagToWaitingTime*365.25, nE[1,1], nE[1,2], birthTime, initState, isMig))
      }
      #cat('\n----------\n')
      return(nE)
    }
    if(ncol(depMatrix) == 2){
      # Compute for each possible destination state a waiting time.
      for(i in 1:length(possTr)){
        tr <- possTr[i]
        destState <-  as.numeric(names(tr))
        cS <- codingScheme[codingScheme[,2] %in% currState, 2+c(1:nSubStates)]
        dS <- codingScheme[codingScheme[,2] %in% destState, 2+c(1:nSubStates)]
        
        tageInYears <- trunc(ageInYears)
        tCalTime <- trunc(1970.001+calTime/365.25)
        
        indRateFctDET <- function(x){
          res <- eval(do.call(tr,
                              args=list(age=tageInYears+x,calTime=tCalTime+x)))
          return(res)
        }
        ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
        detE <- indRateFctDET(ranAccuracyInDays)
        daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
        if (Inf %in% detE) {
          timeToNext <- daysToTrInYears
        } else {
          u <- -log(1-runif(1))
          #cat('It: ',i,'--u: ',u,'\n')
          # Extract individual transition rate (depending on age and calendar time)
          indRateFct <- function(x){
            ageIn <- ageInYears+x
            calIn <- 1970.001+calTime/365.25+x
            res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn)))
            if(TRUE %in% (res<0))
              stop('I have found negative rate value/s for transition: ',tr,'\n
                   This is implausible. Please check this. Simulation has been stopped.\n')
            #cat('x: ',x,' -- res', res,'\n')
            #cat('\n---\n')
            return(res)
          }
          if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
            intHaz <- 0
          } else {
            # Integrated hazard at max. value
            intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
            if(inherits(intHaz, 'try-error')){
              intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
            }
          }
          # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
          if(u<=intHaz){
            invHazFct <- function(x){
              #cat('x: ',x,'\n')
              try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
              #print(try.res)
              if(inherits(try.res, 'try-error')){
                #cat('Seemingly, divergent intergral for ID ',id,
                # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                #  ' for random number: ',u,'\n')
                try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
              }
              #cat('res: ',try.res,'\n-----\n')
              return(try.res)
            }
            # Find random waiting time.
            timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
          } else {
            timeToNext <- Inf
          }
        }
        
        nextEventMatrix[i,1] <- destState
        nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
      }
      #print(nextEventMatrix)
      nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
      if(dim(nE)[1]>1)
        nE <- nE[1,,drop=F]
      if(nE[1,2]!=Inf){
        # Cal. time of next event of individual. (If there is one.)
        tt <- calTime + as.numeric(nE[1,2])
        #print(tt)
        #cat(nE[1,1],'---',tt,'\n')
        # Check whether next event implies school enrollment. If yes, adjust transition time to ensure that the individual
        # enters school at Sept. 1 in the year he/she turns seven.
        if(schoolEnrol){ # Is school enrollment considered in this simulation model? Yes, then continue; otherwise skip this part
          if(isSchoolEnrolment(currState,nE[1,1])){
            enYear <- getYear(tt)
            if(getMonth(tt) <= monthSchoolEnrol) {
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            } else {
              enYear <- enYear+1
              enDate <- getInDays_my(enYear, monthSchoolEnrol)
            }
            diffToEn <- as.numeric(enDate-tt)
            nE[1,2] <- as.numeric(nE[1,2]) + diffToEn
          }
        }
        # Enqueue new event (if there is one).
        queue <<- rbind(queue, c(id, t.clock, currState, currAge - lagToWaitingTime*365.25, nE[1,1], nE[1,2], birthTime, initState, isMig))
      }
      #cat('\n----------\n')
      return(nE)
    }
    else
      for(i in 1:length(possTr)){
        tr <- possTr[i]
        destState <-  as.numeric(names(tr))
        cS <- codingScheme[codingScheme[,2] %in% currState, 2+c(1:nSubStates)]
        dS <- codingScheme[codingScheme[,2] %in% destState, 2+c(1:nSubStates)]
        
        tageInYears <- trunc(ageInYears)
        tCalTime <- trunc(1970.001+calTime/365.25)
        indRateFctDET <- function(x){
          res <- eval(do.call(tr,
                              args=list(age=tageInYears+x,calTime=tCalTime+x)))
          return(res)
        }
        ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
        detE <- indRateFctDET(ranAccuracyInDays)
        daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
        if (Inf %in% detE) {
          timeToNext <- daysToTrInYears
        } else {
          u <- -log(1-runif(1))
          #cat('It: ',i,'--u: ',u,'\n')
          # Extract individual transition rate (depending on age and calendar time)
          indRateFct <- function(x){
            ageIn <- ageInYears+x
            calIn <- 1970.001+calTime/365.25+x
            res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn)))
            if(TRUE %in% (res<0))
              stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
            #cat('x: ',x,' -- res', res,'\n')
            #cat('\n---\n')
            return(res)
          }
          if(sum(indRateFct(0:ran))==0){ # Rate function contains only zeros.
            intHaz <- 0
          } else {
            # Integrated hazard at max. value
            intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
            if(inherits(intHaz, 'try-error')){
              intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
            }
          }
          # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.
          if(u<=intHaz){
            invHazFct <- function(x){
              #cat('x: ',x,'\n')
              try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
              #print(try.res)
              if(inherits(try.res, 'try-error')){
                #cat('Seemingly, divergent intergral for ID ',id,
                # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
                #  ' for random number: ',u,'\n')
                try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
              }
              #cat('res: ',try.res,'\n-----\n')
              return(try.res)
            }
            # Find random waiting time.
            timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root
          } else {
            timeToNext <- Inf
          }
        }
        nextEventMatrix[i,1] <- destState
        nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
      }
    #print(nextEventMatrix)
    nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
    if(dim(nE)[1]>1)
      nE <- nE[1,,drop=F]
    if(nE[1,2]!=Inf){
      # Cal. time of next event of individual. (If there is one.)
      tt <- calTime + as.numeric(nE[1,2])
      #print(tt)
      #cat(nE[1,1],'---',tt,'\n')
      # Check whether next event implies school enrollment. If yes, adjust transition time to ensure that the individual
      # enters school at Sept. 1 in the year he/she turns seven.
      #TODO check if correct! since noone is enrolling currently!
      if(schoolEnrol){ # Is school enrollment considered in this simulation model? Yes, then continue; otherwise skip this part
        if(isSchoolEnrolment(currState,nE[1,1])){
          enYear <- getYear(tt)
          if(getMonth(tt) <= monthSchoolEnrol) {
            enDate <- getInDays_my(enYear, monthSchoolEnrol)
          } else {
            enYear <- enYear+1
            enDate <- getInDays_my(enYear, monthSchoolEnrol)
          }
          diffToEn <- as.numeric(enDate-tt)
          nE[1,2] <- as.numeric(nE[1,2]) + diffToEn
        }
      }
      # Enqueue new event (if there is one).
      queue <<- rbind(queue, c(id, t.clock, currState, currAge - lagToWaitingTime*365.25, nE[1,1], nE[1,2], birthTime, initState, isMig))
    }
    #cat('\n----------\n')
    return(nE)
    
  }
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # E. INITIALIZATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Compute next events for members of starting population
  cat('Initialization ... \n')
  time_init_start = Sys.time()
  print(paste("Starting at: ", time_init_start))
  
  if(length(fertTr)>0){
    fertTrExpanded <- buildFertTrExpanded()
  }
  if(schoolEnrol){
    eduTrExpanded <- buildEduTrExpanded()
  }
  
  birthTimeInDays <- getInDays(initPop[,'birthDate'])
  IN <- matrix(c(initPop[,'ID'], # ID
                 rep(99,nrow(initPop)), # currState (= initState)
                 simStartInDays-birthTimeInDays, # age
                 rep(simStartInDays,nrow(initPop)), # calenderTime
                 birthTimeInDays, # birth time in days
                 rep(99,nrow(initPop)),  # initState
                 rep(0,nrow(initPop))), # is migrant (-> no)
               ncol=7, nrow=nrow(initPop))
  IN[,2] <- IN[,6] <- codingScheme[match(initPop[,'initState'], codingScheme[,1]),2]
  
  if(TRUE %in% (IN[,3]<0)) {
    cat("There are persons born later than simulation starting date in the initial population. Related IDs are: ")
    negAge <- IN[,1][IN[,3]<0]
    for(i in 1:length(negAge)){
      cat(negAge[i]," ")
    }
    stop("Error: Negative age in initial population.")
  }
  if(TRUE %in% (IN[,3]/365.25>maxAge)) {
    cat("In the initial population, there are persons older than maxAge at simulation starting date. Related IDs are: ")
    invalAge <- IN[,1][IN[,3]/365.25>maxAge]
    for(i in 1:length(invalAge)){
      cat(invalAge[i]," ")
    }
    stop("Error: Older than max. age in initial population.")
  }
  
  maxId <- max(IN[,1])
  
  init <- apply(IN, 1, getNextStep)
  
  
  # If immigrants enter the population, compute next events for them.
  if(!is.null(immigrPop)){
    #IM <- data.frame(ID=immigrPop[,'ID'], currState=immigrPop[,'immigrInitState'], age=getAgeInDays(immigrPop[,'immigrDate'],immigrPop[,'birthDate']),
    #                 calTime=getInDays(immigrPop[,'immigrDate']),stringsAsFactors=FALSE)
    
    # Check whether migrants are already born when they migrate
    if(TRUE %in% (immigrPop$immigrDate<immigrPop$birthDate)){
      cat("In the immigration population, there are persons who are not yet born when they immigrate. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate<immigrPop$birthDate]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }
      stop("Error: Not yet born at immigration date.")
    }
    # Check whether all migrants migrate after simulation starting date
    if(TRUE %in% (immigrPop$immigrDate<simHorizon[1])){
      cat("In the immigration population, there are persons who are specified to migrate into the virtual population before simulation starting time. That's against MicSim's concept of migration. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate<simHorizon[1]]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }
      stop("Error: Migrate before simulation starting date.")
    }
    # Check whether all migrants migrate before simulation stopping date
    if(TRUE %in% (immigrPop$immigrDate>simHorizon[2])){
      cat("In the immigration population, there are persons who are specified to migrate into the virtual population after simulation ending time. That's a bit meaningless. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate>simHorizon[2]]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }
      stop("Error: Migrate after simulation ending date.")
    }
    
    immigrTimeInDays <- getInDays(immigrPop$immigrDate)
    birthTimeInDays <- getInDays(immigrPop$birthDate)
    IM <- matrix(c(immigrPop[,'ID'], # ID
                   rep(99,nrow(immigrPop)), # currState (= initState)
                   immigrTimeInDays - birthTimeInDays, #age
                   immigrTimeInDays, # calenderTime
                   birthTimeInDays, # birth time in days
                   rep(99,nrow(immigrPop)),  # initState
                   rep(1,nrow(immigrPop))), # is migrant (-> yes)
                 ncol=7, nrow=nrow(immigrPop))
    IM[,2] <- IM[,6] <- codingScheme[match(immigrPop[,'immigrInitState'], codingScheme[,1]),2]
    
    maxId <- max(IN[,1], IM[,1])
    
    # Check whether all migrants are younger than maxAge when they migrate
    ageIm <- IM[,3]/365.25 # age at immigration in years
    if(TRUE %in% (ageIm>maxAge)){
      cat("In the immigration population, there are persons who older than 'maxAge' when they migrate. That's a bit meaningless. Related IDs are: ")
      invalImAge <- immigrPop$ID[ageIm>maxAge]
      for(i in 1:length(invalImAge)){
        cat(invalImAge[i]," ")
      }
      stop("Error: Migrants in the input data are older than 'maxAge'.")
    }
    
    imit <- apply(IM, 1, getNextStep, isIMInitEvent=T)
    
    immigrInitPop <- immigrPop[,c('ID','birthDate','immigrInitState')]
    colnames(immigrInitPop)[3] <- 'initState'
    initPop <- rbind(initPop, immigrInitPop)
  }
  
  time_init_end = Sys.time()
  print(paste("Ending at: ", time_init_end))
  #p_time = (time_init_end - time_init_start)
  #print(paste("#Time needed for initialization: ", p_time))
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # F. SIMULATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Run simulation either until queue is empty or until simulation horizon has been reached.
  cat('Simulation is running ... \n')
  currYear <- trunc(simHorizon[1]/10000)
  cat('Year: ',currYear,'\n')
  while(nrow(queue)>0 & t.clock <= simStopInDays){
    
    # Sort queue according to soonest event to happen.
    queue <- queue[order(queue[,2] + queue[,6]),,drop=F] # columns: 'ID','currTime','currState','currAge','nextState','timeToNextState','birthtime', 'initState', "isMig"; in queue currTime in days since 01-01-1970
    #print_t.clock <- paste(c(getDay(t.clock), getMonth(t.clock), getYear(t.clock)), collapse="/")
    #print(print_t.clock)
    #print(dim(queue)[1])
    # Enqueue individual who has the soonest event to happen.
    indS <- queue[1,]
    #print(indS)
    # Remove he/she from queue.
    queue <- queue[-1,,drop=F]
    # Set the global clock.
    t.clock <- indS[2] + indS[6] # in days since 01-01-1970
    cY <- getYear(t.clock) # transform days since 01-01-1970 back to years
    # If the global clock exceeds the simulation horizon, stop simulation.
    if(t.clock > simStopInDays)
      break
    if(cY>currYear){
      cat('Year: ',cY,'\n')
      currYear <- cY
    }
    # Age at current transition
    age <- indS[4] + indS[6] # in days
    # Register transition.
    transitions <- rbind(transitions, c(indS[c(1,3,5)], t.clock, age))
    # If current state is not an absorbent one, compute next event.
    if(!indS[5] %in% absStatesNum){
      # Current transition causes a newborn? If yes, add one to simulation population.
      if(length(fertTr)>0){
        if(isBirthEvent(indS[3],indS[5])){
          addNewNewborn_dur(birthTime=t.clock, motherID=indS[1], motherState=indS[5])
        }
      }
      res <- getNextStep(c(indS[c(1,5)], age, t.clock, indS[c(7,8,9)])) # ID, currState, age, calTime, birthtime, initState, isMig
      #print(res)
    }
    #cat('\n-----------\n')
  }
  transitions <- transitions[order(transitions[,1]),,drop=F] # columns: ID, From, To, transitionTime, transitionAge
  
  if (nrow(transitions) == 0){
    
    transitionsOut <- data.frame(ID=initPop[,'ID'], From= rep(NA,nrow(initPop)),
                                 To=rep(NA,nrow(initPop)), transitionTime = rep(NA,nrow(initPop)),
                                 transitionAge = rep(NA,nrow(initPop)), stringsAsFactors = FALSE)
    cat('Simulation has finished.\n')
    cat('Beware that along the simulation horizon the individual/s considered do/es not experience any transition/s.\n')
    cat('------------------\n')
    
  } else {
    
    cat('Simulation has finished.\n------------------\n')
    
    #     ----------------------------------------------------------------------------------------------------------------------
    #     ----------------------------------------------------------------------------------------------------------------------
    #     G. GENERATE OUTPUT
    #     ----------------------------------------------------------------------------------------------------------------------
    #     ----------------------------------------------------------------------------------------------------------------------
    
    indCodesTo <- match(transitions[,2], codingScheme[,2]) # transform numerical state codes to strings according to codingScheme
    transTo <- codingScheme[indCodesTo,1]
    indCodesFrom <- match(transitions[,3], codingScheme[,2])
    transFrom <- codingScheme[indCodesFrom,1]
    transitionsOut <- data.frame(ID=as.numeric(transitions[,1]), From=transTo , To=transFrom,
                                 transitionTime = getInDateFormat(transitions[,4]),
                                 transitionAge = round(transitions[,5]/365.25,2),
                                 stringsAsFactors = FALSE)
  }
  
  pop <- merge(initPop, transitionsOut, all=T, by='ID')
  pop <- pop[order(as.numeric(pop[,1]), as.numeric(pop[,7])),]
  
  if(length(fertTr)>0) {
    colnames(mothers) <- c("motherID", "ID")
    pop <- merge(pop, mothers, by="ID", all.x=TRUE)
    pop <- pop[order(as.numeric(pop[,1]), as.numeric(pop[,7])),]
  }
  
  return(pop)
}


