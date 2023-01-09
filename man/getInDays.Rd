\name{getInDays}
\alias{getInDays}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Get from a date given in the numeric format yyyymmdd the number of days elapsed since 1970-01-01
}
\description{
 Function computes the days that have pasted since 1970-01-01 up to the currDate (in the numeric format yyyymmdd)
}
\usage{
getInDays(currDate)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{currDate}{
    Date given as string of the numeric format yyyymmdd. 
}
}
\value{
 Number of days elapsed since 1970-01-01.
}
\author{
Sabine Zinn
}
\examples{
getInDays(20200826)
}

