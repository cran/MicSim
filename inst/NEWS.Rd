\name{NEWS}
\title{NEWS file for the MicSim package}
\section{Changes in version 1.1.0}{
 \itemize{
   \item Fix two bugs related to the creation of newborns during simulation: (1) Before that version all transitions of the type "1+"->"1+" led to newborns; now only transitions where no other attributes of the state change cause newborns, e.g. "f/M/+1" -> "f/D/1+" caused a newborn before this version. This is not any longer the case. Now only transitions of the form "f/M/1+" -> "f/M/+1" cause a newborn. (2) There was a bug in transforming the fertility matrix to detect birth events; instead of rows, columns were taken as transitions. This has been fixed now.
   \item Fix an inconsistency in the computation of transition ages (in years) related to leap years (now the related function is vectorized).
   \item A further functionality was added to the package: now it is possible to link individuals born during simulation to their mothers (see argument 'reportMothers' in the \code{micSim} or \code{micSimParallel} function). 
  }
}
\section{Changes in version 1.0.15}{
 \itemize{
   \item Fix bug in parallel computing function concerning case when there are migrants in the model (this bug      resulted from fix to version 1.0.14). 
  }
}
\section{Changes in version 1.0.14}{
 \itemize{
   \item Fix bug in parallel computing function concerning case when there are no migrants in the model. 
  }
}
\section{Changes in version 1.0.13}{
 \itemize{
   \item Update of chron package caused an error due to new date format. This problem has been fixed. 
  }
}
\section{Changes in version 1.0.12}{
 \itemize{
   \item Fix bug concerning the transition age of immigrants.
  }
}
\section{Changes in version 1.0.11}{
 \itemize{
   \item Up to version 1.0.10 deterministic events could only be defined to occur after whole numbered years. Now, they can also be defined to happen after whole numbered days.
  }
}
\section{Changes in version 1.0.10}{
 \itemize{
   \item Fix bug concerning the uniqueness of IDs of newborns when simulating natives and migrants on several cores, i.e., when using the function `micSimParallel'.
   \item Slight changes in the description of the (life-course) model used in MicSim.
  }
}

\section{Changes in version 1.0.9}{
 \itemize{
   \item Add error message arising if negative rate values appear.
   \item Besides not defining impossible transitions in the transition matrix, they can now also be specified by giving only zero rate values for the respective rate function/s. 
  }
}
\section{Changes in version 1.0.8}{
 \itemize{
   \item Missing package dependency has been added. MicSim requires `rlecuyer'.
  }
}
\section{Changes in version 1.0.7}{
 \itemize{
   \item Add a new reference.
   \item Modify function `micSim' that it can also deal with a single life course.
  }
}
\section{Changes in version 1.0.6}{
 \itemize{
   \item Fix bug in function `convertToLongFormat' leading to wrong ordering of the event sequences of individuals born during simulation.
  }
}
\section{Changes in version 1.0.5}{
 \itemize{
   \item Fix minor conversion problem in function `convertToWideFormat', related to request for NA values.
   \item Now function `buildTransitionMatrix' can handle covariate dependent mortality rates of every kind, e.g.,     mortality rates depending on health status and gender.
  }
}
\section{Changes in version 1.0.4}{
 \itemize{
   \item Add function `convertToWideFormat' to convert simulation output into wide format.
   \item Slight changes in the examples to make them more feasible. Concretely, in several rates definitions I adjusted the year in which dynamics are assumed to change.
   \item Fix a bug in the optional features of the `micSim' function: Children who are part of the base population had been enrolled one year to late. Now, enrollment dates are OK.
  }
}

\section{Changes in version 1.0.3}{
 \itemize{
   \item Refine funtionality of `micSim' function to handle deterministc events (marked by Inf rates).
   \item Fix bug in function `buildTransitionMatrix'.
  }
}
\section{Changes in version 1.0.2}{
 \itemize{
   \item Add News file.
  }
}
\section{Changes in version 1.0.1}{
 \itemize{
   \item Extend function `buildTransitionMatrix' to facilitate using transition rates which depend on covariates such as gender and marital status.
   \item Correct typing error in Details section of help page of function `buildTransitionMatrix'.
   \item Fix a problem in function `micSim' occuring if the simulation is asked not to produce newborns (concerns argument `fertTr=c()').
  }
}