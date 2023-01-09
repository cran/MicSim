\name{getAgeInDays}
\alias{getAgeInDays}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Get from a given date the age in days
}
\description{
 Function computes for a given date the correct age in days.
}
\usage{
getAgeInDays(currDate, birthDate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{currDate}{
    Reference date given as string of the format "yyyymmdd". 
}
  \item{birthDate}{
    Birth date given as string of the format "yyyymmdd". 
}
}
\value{
 Correct age at the specific date \code{currDate} in days 
}
\author{
Sabine Zinn
}
\examples{
getAgeInDays("20200826", "19800605")
}

