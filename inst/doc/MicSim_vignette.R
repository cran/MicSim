## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
# Load library
library(MicSim)

## -----------------------------------------------------------------------------
# Defining simulation horizon
startDate <- 20140101 # yyyymmdd
endDate   <- 20181231 # yyyymmdd
simHorizon <- c(startDate=startDate, endDate=endDate)

# Seed for random number generator
set.seed(12)

# Definition of maximal age (sharp, i.e. max age is 100.0)
maxAge <- 100

# Definition of state space
sex <- c("m","f")                     
country_R <- c("NL","ES","SE")
fert <- c("0","1","2","3","4+","999")           
stateSpace <- expand.grid(sex=sex,country_R=country_R,fert=fert)

# Definition of nonabsorbing and absorbing states
absStates <- c("dead","rest")   

## ----setup--------------------------------------------------------------------
# initial population included in the MicSim package
initPop <- initPopMigrExp
head(initPop)

# population of immigrants entering to virtual population over simulation time
immigrPop <- immigrPopMigrExp # included in the MicSim package
head(immigrPop)

# Definition of initial states for newborns
fixInitStates <- 2 # give indices for attribute/substate that will be taken over from the mother, here: nat (nationality)
varInitStates <- rbind(c("m","ES","0"), c("f","ES","0"), # to have possibility to define distinct sex ratios in distinct countries, 
                       c("m","NL","0"), c("f","NL","0"), # hold substate that are inherited by mother in the init state (i.e. nationality)
                       c("m","SE","0"), c("f","SE","0")) 
initStatesProb <- c(0.5151,0.4849, # probabilities for female / male newborns for each nationality separately
                    0.5124,0.4876,
                    0.5146,0.4854)

## -----------------------------------------------------------------------------
# Load rates from data included in the MicSim package (one column for one transition)
rates <- migrExpRates

# Define function to easily transform rates in function format required by MicSim 
require(glue)
for (i in 1:length(names(rates))) {
  script_var = names(rates[i])
  eval(parse(text = glue("{script_var} <- unlist(as.numeric(rates[,{i}]))")))
}

# --------------------------------------------------------------------------
# Fertility Rates
# --------------------------------------------------------------------------

fertRates_NL_0_1 <- function(age,calTime, duration){
  rate <- fert_NL_0_1[as.integer(age)+1]
  return(rate)  
}

fertRates_NL_1_2 <- function(age,calTime, duration){
  rate <- fert_NL_1_2[as.integer(age)+1]
  return(rate)  
}

fertRates_NL_2_3 <- function(age,calTime, duration){
  rate <- fert_NL_2_3[as.integer(age)+1]
  return(rate)  
}

fertRates_NL_3_4 <- function(age,calTime, duration){
  rate <- fert_NL_3_4[as.integer(age)-+1]
  return(rate)  
}

fertRates_ES_0_1 <- function(age,calTime, duration){
  rate <- fert_ES_0_1[as.integer(age)+1]
  return(rate)  
}

fertRates_ES_1_2 <- function(age,calTime, duration){
  rate <- fert_ES_1_2[as.integer(age)+1]
  return(rate)  
}

fertRates_ES_2_3 <- function(age,calTime, duration){
  rate <- fert_ES_2_3[as.integer(age)+1]
  return(rate)  
}

fertRates_ES_3_4 <- function(age,calTime, duration){
  rate <- fert_ES_3_4[as.integer(age)+1]
  return(rate)  
}

fertRates_SE_0_1 <- function(age,calTime, duration){
  rate <- fert_SE_0_1[as.integer(age)+1]
  return(rate)  
}

fertRates_SE_1_2 <- function(age,calTime, duration){
  rate <- fert_SE_1_2[as.integer(age)+1]
  return(rate)  
}

fertRates_SE_2_3 <- function(age,calTime, duration){
  rate <- fert_SE_2_3[as.integer(age)+1]
  return(rate)  
}

fertRates_SE_3_4 <- function(age,calTime, duration){
  rate <- fert_SE_3_4[as.integer(age)+1]
  return(rate)  
}

# --------------------------------------------------------------------------
# Internal Migration Rates
# --------------------------------------------------------------------------

`%!in%` <- Negate(`%in%`)

for(i in 1:length(country_R)) {
  other_provinces = country_R[which(country_R %!in% country_R[i])]
  for(k in 1:length(other_provinces)) {
    eq = paste(sprintf('%s_%s_rates', glue("{country_R[i]}"),  glue("{other_provinces[k]}")), 
               '<- function(age,calTime, duration)', '{',
               sprintf('rate <- rate_%s_%s[as.integer(age)+1]', glue("{country_R[i]}"), 
               glue("{other_provinces[k]}")), "\n ", 'return(rate)','}')
    eval(parse(text = eq))
  }
}

# --------------------------------------------------------------------------
# Mortality Rates
# --------------------------------------------------------------------------

# Female mortality
for(i in 1:length(country_R)) {
  eq = paste(sprintf('mortRates_f_%s', glue("{country_R[i]}")), '<- function(age,calTime, duration)',
             '{',
             sprintf('rate <- mort_f_%s[as.integer(age)+1]', glue("{country_R[i]}")),
             "\n ",
             'return(rate)',
             '}')
  eval(parse(text = eq))
}

# Male mortality
for(i in 1:length(country_R)) {
  eq = paste(sprintf('mortRates_m_%s', glue("{country_R[i]}")), '<- function(age,calTime, duration)',
             '{',
             sprintf('rate <- mort_m_%s[as.integer(age)+1]', glue("{country_R[i]}")),
             "\n ",
             'return(rate)',
             '}')
  eval(parse(text = eq))
}

# ---------------------------------------------------------------------------
# Emigration rates
# ---------------------------------------------------------------------------

# Emigration rates for females
for(i in 1:length(country_R)) {
  eq = paste(sprintf('emigrRates_f_%s', glue("{country_R[i]}")), '<- function(age,calTime, duration)',
             '{',
             sprintf('rate <- emig_f_%s[as.integer(age)+1]', glue("{country_R[i]}")),
             "\n ",
             'return(rate)',
             '}')
  eval(parse(text = eq))
}

# Emigration rates for males
for(i in 1:length(country_R)) {
  eq = paste(sprintf('emigrRates_m_%s', glue("{country_R[i]}")), '<- function(age,calTime, duration)',
             '{',
             sprintf('rate <- emig_m_%s[as.integer(age)+1]', glue("{country_R[i]}")),
             "\n ",
             'return(rate)',
             '}')
  eval(parse(text = eq))
}


## -----------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Transition matrix for fertility
# ---------------------------------------------------------------------------
fertTrMatrix <- cbind(c("f/ES/0->f/ES/1","f/ES/1->f/ES/2","f/ES/2->f/ES/3","f/ES/3->f/ES/4+",
                        "f/SE/0->f/SE/1","f/SE/1->f/SE/2","f/SE/2->f/SE/3","f/SE/3->f/SE/4+",
                        "f/NL/0->f/NL/1","f/NL/1->f/NL/2","f/NL/2->f/NL/3","f/NL/3->f/NL/4+"),
                      c("fertRates_ES_0_1", "fertRates_ES_1_2", "fertRates_ES_2_3", "fertRates_ES_3_4",
                        "fertRates_SE_0_1", "fertRates_SE_1_2", "fertRates_SE_2_3", "fertRates_SE_3_4",
                        "fertRates_NL_0_1", "fertRates_NL_1_2", "fertRates_NL_2_3", "fertRates_NL_3_4"))

# ---------------------------------------------------------------------------
# Transition matrix for migration
# ---------------------------------------------------------------------------
# First: make names for transition matrix, i.e. stateOfOrigin->stateOfDestination 
testo1 <- "intmigTrMatrix <- cbind(c("
for(i in 1:length(country_R)) {
  for(m in 1:length(country_R)) {
    if(m != i) {
      eq1 = paste(sprintf("'%s->%s'", glue("{country_R[i]}"), glue("{country_R[m]}")))
      if (i == length(country_R) & m == (i-1)) {
        testo1 = paste(testo1,eq1)
      } else {
        testo1 = paste(testo1 ,eq1, ",")
      }
    }
    if (m == i & m == length(country_R)){
      testo1 = paste(testo1, "),")
      
    }
  }
}

#Second: set names for transition functions
testo1 = paste (testo1,"c(")
for(i in 1:length(country_R)) {
  for(m in 1:length(country_R)) {
    if(m != i) {
      eq1 = paste(sprintf("'%s_%s_rates'", glue("{country_R[i]}"), glue("{country_R[m]}")))
      if (i == length(country_R) & m == (i-1)) {
        testo1 = paste(testo1,eq1)
      } else {
        testo1 = paste(testo1 ,eq1, ",")
      }
    }
    if (m == i & m == length(country_R)){
      testo1 = paste(testo1, "))")
    }
  }
}
eval(parse(text = testo1))

# ---------------------------------------------------------------------------
# Transition matrix for mortality and emigration
# ---------------------------------------------------------------------------

testo <- "absTransitions <- rbind("
for(i in 1:length(country_R)) {
  for(m in 1:length(sex)) {
    eq1 = paste(sprintf("c('%s/%s/dead', 'mortRates_%s_%s')", glue("{sex[m]}"), glue("{country_R[i]}") ,
                        glue("{sex[m]}"), glue("{country_R[i]}")),',',
                sprintf("c('%s/%s/rest', 'emigrRates_%s_%s')", glue("{sex[m]}"),glue("{country_R[i]}") ,
                        glue("{sex[m]}"), glue("{country_R[i]}")))
    if(i == length(country_R) & m == length(sex)) {
      testo = paste(testo, eq1, ")")
    } else {
      testo = paste(testo,eq1, ",")
    }
  }
}
eval(parse(text = testo))

# ---------------------------------------------------------------------------
# Combine all
# ---------------------------------------------------------------------------

allTransitions <- rbind(fertTrMatrix, intmigTrMatrix)
transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions,
                                          absTransitions=absTransitions, 
                                          stateSpace=stateSpace)

# ---------------------------------------------------------------------------
# Define transitions triggering a birth event
# ---------------------------------------------------------------------------

fertTr <- fertTrMatrix[,1]

## -----------------------------------------------------------------------------
pop <- micSim(initPop=initPop[1:500,], immigrPop=immigrPop[1:100,],
              transitionMatrix=transitionMatrix, absStates=absStates,
              varInitStates=varInitStates, initStatesProb=initStatesProb,
              fixInitStates=fixInitStates,
              maxAge=maxAge, simHorizon=simHorizon,fertTr=fertTr)
head(pop)


## -----------------------------------------------------------------------------
popLong <- convertToLongFormat(pop, migr=TRUE)
head(popLong)


## -----------------------------------------------------------------------------
popWide <- convertToWideFormat(pop)
head(popWide)

## -----------------------------------------------------------------------------

cores <- 3
seeds <- c(34,145,97)
pop <- micSimParallel(initPop=initPop, immigrPop=immigrPop,
                      transitionMatrix=transitionMatrix, absStates=absStates,
                      varInitStates=varInitStates, initStatesProb=initStatesProb,
                      fixInitStates=fixInitStates,
                      maxAge=maxAge, simHorizon=simHorizon,fertTr=fertTr,
                      cores=cores, seeds=seeds)
head(pop)

