\name{convertToLongFormat}
\alias{convertToLongFormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshaping microsimulation output into long format
}
\description{
The function reshapes the output given by \link{micSim} or by \link{micSimParallel} into long format. In long format, the data comprises for each episode which an individual experiences one row. 
}
\usage{
convertToLongFormat(pop,migr=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop}{
   The data frame \code{pop} contains the whole synthetic population considered during simulation including all events generated. For each individidual \code{pop} contains as many rows as the individual performed transitions during simulation.  
}
  \item{migr}{
   A logical variable indicating whether the simulation model considers immigration. The default setting is "no immigration considered": \code{migr=FALSE}. 
  }
}
\value{
  A data frame comprising the microsimulation output in long format.
  
  - \code{ID} is the unique numerical person identifier of an individual.
  
  - \code{birthDate} is the birth date of an individual.
  
  - The variables \code{Tstart} and \code{Tstop} mark the start und the ending dates of episodes.   
  
  - \code{statusEntry} specifies whether the entry into an episode has been observed. Value "1" marks an observed entry and "0" marks a left truncated episode.
 
  - \code{statusExit} specifies whether a transition between two states or right censoring completed an episode. Value "1" indicates a transition and "0" a censoring event. 
 
  - \code{OD} names the transition which completed an episode. Here, right censoring is marked by "cens".  
 
  - \code{ns} gives the number of episodes an individual has passed.
  
  - \code{Episode} enumerates the episodes an individual has passed.

  - The last columns of the data frame contain for each individual and episode the values of the state variables during that episode such as `sex', `education', etc. 
   
 - Birth and transition times are given as calendar dates in form of chron objects.
}
\details{
 \code{convertToLongFormat} uses information from the definition of the microsimulation model. In particular, it uses \code{stateSpace}, \code{absTransitions}, \code{allTransitions}, \code{simHorizon}, and optionally \code{immigrPop}. (For a description of these objects see \link{micSim}.) \code{stateSpace}, \code{absTransitions}, \code{allTransitions},  \code{simHorizon}, and \code{immigrPop} are globally defined, i.e., they are already part of the workspace. Thus, they do not have to be given to \code{convertToLongFormat} as extra input parameters.  
}
\author{
Sabine Zinn
}

\examples{

# Run microsimulation before, e.g., the complex example described on the 
# help page of the function "micSim".
\dontrun{
pop <- micSim(initPop, immigrPop, transitionMatrix, absStates, initStates, initStatesProb, 
        maxAge, simHorizon, fertTr)
popLong <- convertToLongFormat(pop,migr=TRUE)
}

}


