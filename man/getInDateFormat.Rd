\name{getInDateFormat}
\alias{getInDateFormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Get date in the format 'yyyyddmm' from days elapsed since 01-01-1970
}
\description{
 Function generates from days elapsed since 01-01-1970 the date in the string format 'yyyyddmm'.
}
\usage{
getInDateFormat(daysSince01011970)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{daysSince01011970}{
    Days elapsed since 1970-01-01
}
}
\value{
 Date in string format 'yyyyddmm' from days elapsed since 01-01-1970
}
\author{
Sabine Zinn
}
\examples{
getInDateFormat(2561)
}

