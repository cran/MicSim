\name{MicSim-package}
\alias{MicSim-package}
\alias{MicSim}
\docType{package}
\title{
MicSim: Continuous-time microsimulation for population projection
}
\description{
In demography, the central device of microsimulations is the life course of an individual, which is defined by the sequence of states that the individual visits over time, and the waiting times between these state transitions. Modelling and simulating the life courses of a representative share of population members allows mapping population dynamics on a very detailed scale.

A standard approach to describe individual behavior is a continuous-time multi-state model. A multi-state model is a stochastic process that at any point in time occupies one out of a set of discrete states. These states summarize the demographically relevant categories an individual can belong to. Generally, the state space is determined by the problem to be studied, but commonly it will at least comprise the elementary demographic characteristics of sex and marital status. One element always present in the state space is "dead", a risk to which each individual is always exposed to.

In demographic microsimulations life-courses usually evolve along two time scales: individual age and calendar time. A possible third time scale is the time that an individual has already spent in his/her current demographic state, e.g., the time that has elapsed since the individual's wedding.A demographic event implies a change in the state of an individual.It should be emphasized that age runs parallel to the process time in the model, and therefore birthdays, i.e., completion of another year of life, is not an event in itself. 

A common way to characterize an individual life-course is via a trajectory of a stochastic process from the family of Markovian processes, where the process time maps the time span over which we "observe" an individual life-course.If in the microsimulation model life-courses are specified only along age and calendar time, the accordant class used is the class of non-homogeneous continuous-time Markov chains. If in addition, the time that has elapsed since the last event is considered, non-homogeneous semi-Markov processes are used to describe individual life-courses.

The transition rates (also denoted as hazard rates or intensities) of Markovian processes are their key quantities. Once they are known one can compute the distribution functions of sojourn times and thus simulate synthetic life-courses. That is, to run a microsimulation model, for all transitions and time scales considered transition rates have to be provided.

}
\details{
\tabular{ll}{
Package: \tab MicSim\cr
Type: \tab Package\cr
Version: \tab 1.0.9\cr
Date: \tab 2015-07-23\cr
License: \tab GPL-2\cr
}
}
\author{
Sabine Zinn

Maintainer: sabine.zinn@lifbi.de
}
\references{
S. Zinn (2014). The MicSim Package of R: An Entry-Level Toolkit for Continuous-Time Microsimulation. In International Journal of Microsimulation 7(3), 3-32.

Willekens, F. and H. Putter eds. (2014) Multistate event history analysis. Special Collection of Demographic Research, forthcoming. 

Zinn, S., J.Himmelspach, J.Gampe, and A.M. Uhrmacher (2009). Mic-core: A tool for microsimulation. In M.D. Rossetti, R.R. Hill, B.Johansson, A.Dunkin, and R.G. Ingalls (Eds.), Proceedings of the 2009 Winter Simulation Conference.

Gampe, J. and S. Zinn (2007). Description of the microsimulation model. Technical report, MicMac project, MPIDR, Rostock.
}
\keyword{ package }

