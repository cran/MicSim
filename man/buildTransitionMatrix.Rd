\name{buildTransitionMatrix}
\alias{buildTransitionMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Determining transition pattern and transition functions 
}
\description{
 The function \code{buildTransitionMatrix} supports the constructing of the `transition matrix', which determines the transition pattern of the microsimulation model. The actual microsimulation is performed by \link{micSim} (sequentially) or by \link{micSimParallel} (parallel computing).
}
\usage{
buildTransitionMatrix(allTransitions, absTransitions, stateSpace)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{allTransitions}{
    A matrix comprising all possible transitions between values of state variables in the first column and in the second column the names of the functions defining the corresponding transition rates.
}
  \item{absTransitions}{
  A matrix comprising the names of the absorbing states which individuals are always exposed to (such as "dead"" and emigrated labeled as "rest") in the first column and in the second column the names of the functions defining the corresponding transition rates. 
}
  \item{stateSpace}{
 A matrix comprising all nonabsorbing states considered during simulation. 
}
}
\details{
 The function \code{buildTransitionMatrix} is an auxiliary function for building the transition matrix required to run the microsimulation using \link{micSim} or \link{micSimParallel}. 
 
In \code{stateSpace} all state variables considered during simulation including their values have to be defined. Values are always described using labels. For example, label "M" for being married. Each column of \code{stateSpace} refers to one state variable considered and each row refers to one state of the state space. Apart from "m" and "f" reserved for male and female (state variable: gender) and "no" and "low" reserved for no education and elementary school attended (state variable: educational attainment), labels can be set arbitrarily.

 Each element of the first column of \code{allTransitions} has to be of the form "A->B" with indicating "A" the starting value of a transition and "B" the arrival value. ("->" is the placeholder defined to mark a transition.) For example, "0" (childless) describes the starting value of the transition marking a first birth event and "1" (first child) its arrival value. All value labels used have to be identical to the value labels of the state variables specifying the simulation model.
  
All absorbing states listed in the first column of \code{absTransitions} have to be given as strings such as "dead" for being dead or "rest" for emigrated. Since dying is a competing risk all individuals are always exposed to, "dead" is a mandatory part of \code{absTransitions}.

All transitions can be defined to depend on several state variables. For example, a divorce rate depends on gender and on the fertility status. Therefore, the starting value and the arrival value of a transition have to be specified as a combination of the considered attributes, separated by a forward slash and in accordance with the ordering of the state variables in the state space. For example, "f/A->f/B" describes a female specific transition from "A" to "B" and "f/M/1 -> f/D/1" might describe a mother's (indicated by "1") transition from "M" (e.g., married) to "D" (e.g., divorced).
For absorbing states, a prefix indicates the attributes on which a transition is assumed to depend (also separated by forward slashs), e.g., "f/dead" and "m/dead" describe gender specific mortality transitions and "f/M/dead" and "m/M/dead" indicate gender specific mortality rates for married persons. 
}
\value{
  The \code{transitionMatrix} that is mandatory to perform a microsimulation run by \link{micSim} (sequentially) or by \link{micSimParallel} (parallel computing) is returned. The matrix has as many rows as the simulation model comprises nonabsorbing states and as many columns as the simulation model comprises absorbing and nonabsorbing states. The rows indicate starting states of transitions and the columns signify arrival states. At positions indicating impossible transitions, the matrix contains zeros. Otherwise the name of the function defining the respective transition rates is given.  
}
\author{
 Sabine Zinn
}
\examples{

#########################################################################################
# 1. Example: Transition rates are specified to depend on only one state variable
#########################################################################################
\donttest{
# Defintion of state space, i.e., nonabsorbing and absorbing states
sex <- c("m","f")                     
fert <- c("0","1","2","3+")           
marital <- c("NM","M","D","W")        
edu <- c("no","low","med","high")   
stateSpace <- expand.grid(sex=sex,fert=fert,marital=marital,edu=edu)

# Possible transitions indicating fertility behavior are "0->1", "1->2", "2->3+", 
# and "3+->3+". Here, "->" is the defined placeholder defining a transition.
# `fert1Rates' marks the name of the function defining the transition rates to 
# parity one and `fert2Rates' marks the name of the function defining the transition 
# rates to higher parities. 
# Note: The functions `fert1Rates' and `fert1Rates' are transition rate functions 
# defined by the user. Their naming depends on the user's choice.
fertTrMatrix <- cbind(c("0->1","1->2","2->3+","3+->3+"),                         
                c("fert1Rates", "fert2Rates", "fert2Rates","fert2Rates"))

# Possible transitions indicating changes in the marital status are "NM->M", "M->D", 
# "M->W", "D->M", and "W->M".
# `marriage1Rates' marks the name of the function defining the transition rates for first 
# marriage and `marriage2Rates' marks the name of the function defining the transition rates 
# for further marriages. `divorceRates' marks the name of the function defining divorce 
# rates and `widowhoodRates' marks the name of the function describing transition rates to 
# widowhood. 
# Note: The functions `marriage1Rates',`marriage2Rates', `divorceRates', and 
# `widowhoodRates' are transition rate functions defined by the user. 
# Their naming depends on the user's choice.
maritalTrMatrix <- cbind(c("NM->M","M->D","M->W","D->M","W->M"),              
                   c("marriage1Rates","divorceRates","widowhoodRates","marriage2Rates",
                     "marriage2Rates"))
                   
# Possible transitions indicating changes in the educational attainment are "no->low", 
# "low->med", and "med->high". 
# `noToLowEduRates' marks the name of the function defining transition rates for accessing 
# primary education, `noToLowEduRates' marks the name of the function defining transition 
# rates for graduating with a lower secondary education, and `medToHighEduRates' marks the 
# name of the function defining transition rates for graduating with a higher secondary 
# education.
# Note: The functions `noToLowEduRates',`noToLowEduRates', and `medToHighEduRates' are 
# transition rate functions defined by the user. Their naming depends on the user's 
# choice.                  
eduTrMatrix <- cbind(c("no->low","low->med","med->high"),
               c("noToLowEduRates","noToLowEduRates","medToHighEduRates")) 
               
# Combine all possible transitions and the related transition function into one matrix.          
allTransitions <- rbind(fertTrMatrix, maritalTrMatrix, eduTrMatrix)

# Possible absorbing states are `dead' and `rest'. (The latter indicates leaving the 
# population because of emigration). The accordant transition rate functions are named 
# `mortRates' and `emigrRates'. (Again, naming is up to the user.)
absTransitions <- rbind(c("dead","mortRates"),c("rest","emigrRates"))

# Construct `transition matrix'.
transitionMatrix <- buildTransitionMatrix(allTransitions,absTransitions,stateSpace)

#########################################################################################
# 2. Example: Transition rates are gender specific
#########################################################################################
# Defintion of nonabsorbing and absorbing states
sex <- c("m","f")                               
stateX <- c("H","P")        
stateSpace <- expand.grid(sex=sex,stateX=stateX)
absStates <- c("dead")   

# Transitions indicating changes in `stateX'. 
# We assume distinct transition rates for females and males.
# Note: The functions `ratesHP_f',`ratesHP_m', `ratesPH_f', and 
# `ratesPH_m' are transition rate functions defined by the user. 
trMatrix_f <- cbind(c("f/H->f/P","f/P->f/H"),c("ratesHP_f", "ratesPH_f"))
trMatrix_m <- cbind(c("m/H->m/P","m/P->m/H"),c("ratesHP_m", "ratesPH_m"))
allTransitions <- rbind(trMatrix_f,trMatrix_m)

# We assume gender specific mortality rates.
# Note: The naming and specification of the respective mortality rate functions 
# `mortRates_f' and `mortRates_m' depend on the user.
absTransitions <- rbind(c("f/dead","mortRates_f"), c("m/dead","mortRates_m"))

transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions,
                      absTransitions=absTransitions, stateSpace=stateSpace)

}

}

