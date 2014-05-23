\name{convertToWideFormat}
\alias{convertToWideFormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reshaping microsimulation output into wide format
}
\description{
The function reshapes the output given by \link{micSim} or by \link{micSimParallel} into wide format. In wide format, the data comprises for each episode which an individual experiences additional column entries. 
}
\usage{
convertToWideFormat(pop)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pop}{
   The data frame \code{pop} contains the whole synthetic population considered during simulation including all events generated. For each individidual \code{pop} contains as many rows as the individual performed transitions during simulation.  
}
}
\value{
  A data frame comprising the microsimulation output in wide format.
  
  - \code{ID} is the unique numerical person identifier of an individual.
  
  - \code{birthDate} is the birth date of an individual.
  
  - \code{initState} is the state in which an individual initially entered the virtual population of the simulation.  
  
  - \code{ns} gives the number of (completed) episodes an individual has passed.
  
  - The variables \code{From.i} and \code{To.i} mark the start und the arrival state of the transition corresponding to episode \code{i}. The variables \code{transitionTime.i} and \code{transitionAge.i} give the corresponding transition time and age. The enumerator \code{i} ranges from 1 to the maximal number of transitions which an individual experienced during simulation. Only completed episodes are counted.  
}
\author{
Sabine Zinn
}

\examples{

# Run microsimulation before, e.g., the complex example described on the 
# help page of the function "micSim".
\dontrun{
pop <- micSim(initPop, immigrPop, transitionMatrix, absStates, initStates,    
    initStatesProb, maxAge, simHorizon, fertTr)
popWide <- convertToWideFormat(pop)
}

}


