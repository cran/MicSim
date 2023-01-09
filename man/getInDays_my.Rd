\name{getInDays_my}
\alias{getInDays_my}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Get the number of days that have pasted from 1970-01-01 until 'yyyymm11'.
}
\description{
 Function computes the number of days that have pasted from 1970-01-01 until 'yyyymm11'.
}
\usage{
getInDays_my(year, month)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{year}{
    Year for which days elapsed should be computed, i.e., the yyyy in 'yyyymm11' 
}
  \item{month}{
    Month for which days elapsed should be computed, i.e., the mm in 'yyyymm11' 
}
}
\value{
 Number of days that have pasted from 1970-01-01 until 'yyyymm11' 
}
\author{
Sabine Zinn
}
\examples{
getInDays_my(2020, 12)
}

