\name{getMonth}
\alias{getMonth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Get the month in a year from days elapsed since 01-01-1970
}
\description{
 Function computes from days elapsed since 01-01-1970 the related month a year.
}
\usage{
getMonth(daysSince01011970)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{daysSince01011970}{
    Days elapsed since 1970-01-01
}
}
\value{
 Month in a year computed from days elapsed since 01-01-1970
}
\author{
Sabine Zinn
}
\examples{
getMonth(2561)
}

