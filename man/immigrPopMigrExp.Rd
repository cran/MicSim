\name{immigrPopMigrExp}
\alias{immigrPopMigrExp}
\docType{data}
\title{
  One possible population of migrants for the MicSim package.
}
\description{
Population of migrants for the MicSim package with 3758 migrants and migration dates between 01-01-2014 and 30-12-2018.
}
\usage{data("immigrPopMigrExp")}
\format{
  A data frame with information of birthdates, date of immigration and state at   immigration of 3758 migrants.
  \describe{
    \item{\code{ID}}{Personal identifier}
    \item{\code{immigrDate}}{Immigration date}
    \item{\code{birthDate}}{Birth date of migrants}
    \item{\code{immigrInitState}}{State at immigration}
  }
}
\details{
This is a example data set for the MicSim package. The population of migrants is already in the format that is required by the package. For more details on this see \link{micSim}. The state space for the states and the related state domains are defined in the vignette of this package. The related application is also part of the vignette.
}
\source{
European Commission
}
\examples{
data(immigrPopMigrExp)
}

