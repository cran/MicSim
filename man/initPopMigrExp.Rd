\name{initPopMigrExp}
\alias{initPopMigrExp}
\docType{data}
\title{
 One possible initial population for the MicSim package.
}
\description{
Initial population for the MicSim package with 72965 persons with birthdates from 08-03-1914 to 30-12-2013.
}
\usage{data("initPopMigrExp")}
\format{
  A data frame with information of birthdates, date of immigration and state at   immigration of 3758 migrants.
  \describe{
    \item{\code{ID}}{Personal identifier}
    \item{\code{birthDate}}{Birth date of entity}
    \item{\code{initState}}{Initial state}
  }
}
\details{
This is a example data set for the MicSim package. The initial population is already in the format that is required by the package. For more details on this see \link{micSim}.
The state space for the states and the related state domains are defined in the vignette of this package. The related application is also part of the vignette.
}
\source{
European Commission
}
\examples{
data(initPopMigrExp)
}
