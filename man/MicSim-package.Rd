\name{MicSim-package}
\alias{MicSim-package}
\alias{MicSim}
\docType{package}
\title{
MicSim: Continuous-time microsimulation for population projection
}
\description{
In life sciences, the central device of microsimulations is the life-course of an individual, which is defined by the sequence of states that the individual visits over time, and the waiting times between these state transitions. Modelling and simulating the life courses of a representative share of population members allows mapping population dynamics on a very detailed scale.

A standard approach to describe individual behavior is a continuous-time multi-state model. A multi-state model is a stochastic process that at any point in time occupies one out of a set of discrete states. These states summarize the demographically relevant categories an individual can belong to. Generally, the state space is determined by the problem to be studied, but commonly it will at least comprise the elementary demographic characteristics of sex and marital status. One element always present in the state space is "dead", a risk to which each individual is always exposed to.

In (demographic) microsimulations life-courses usually evolve along two time scales: individual age and calendar time. A possible third time scale is the time that an individual has already spent in his/her current state, e.g., the time that has elapsed since an individual's wedding. An event implies a change in the state of an individual. Age always runs parallel to the process time of a model. Therefore birthday, i.e., the completion of another year of life, is not an event in itself. 

A common way to characterize an individual life-course is via a trajectory of a stochastic process from the family of Markovian processes, where the process time maps the time span over which we "observe" an individual life-course. The MicSim package uses time-inhomogeneuous Markov models to describe individual life-courses. That way, transition intensities can vary at each point in time, i.e. are not assumed to be constant for predefined time intervals (such as whole years).

The transition intensities (also denoted as hazard rates or transition rates) of Markovian processes are their key quantities. Once they are known one can compute the distribution functions of sojourn times and thus simulate synthetic life-courses. That is, to run a microsimulation model, for all transitions and time scales considered transition rates have to be provided. A whole bunch of statistical estimation approaches exist to estimate transition rates from (e.g., register, survey, panel) data. Furthermore, also methods for approximating transition rates from probabilities are available, e.g. by assuming that they are constant in the interval covered by a probability, yielding the so called exponential model. More details on this are given in the description of the `micsim' function of this package, which is the actual workhorse of this toolkit. 
}
\details{
\tabular{ll}{
Package: \tab MicSim\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Date: \tab 2022-11-01\cr
License: \tab GPL-2\cr
}
}
\author{
Sabine Zinn

Maintainer: szinn@diw.de
}
\references{
S. Zinn (2014). The MicSim Package of R: An Entry-Level Toolkit for Continuous-Time Microsimulation. In International Journal of Microsimulation 7(3), 3-32.

Willekens, F., & Putter, H. (2014). Software for multistate analysis. Demographic Research, 31, 381-420. 

}
\keyword{ package }

