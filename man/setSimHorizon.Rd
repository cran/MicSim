\name{setSimHorizon}
\alias{setSimHorizon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Set simulation horizon
}
\description{
 The function sets the simulation horizon of the microsimulation. The actual microsimulation is performed by \link{micSim} (sequentially) or by \link{micSimParallel} (parallel computing). 
}
\usage{
setSimHorizon(startDate, endDate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{startDate}{
    Starting date of simulation given as string of the format "dd/mm/yyyy". 
}
  \item{endDate}{
    End date of simulation given as string of the format "dd/mm/yyyy". 
}
}
\details{
The starting date has to precede the ending date.
}
\value{
 A vector of two chron objects indicating the simulation horizon of the simulation.
}
\author{
Sabine Zinn
}
\examples{
setSimHorizon(startDate="01/01/2000", endDate="31/12/2010")
}

