\name{NEWS}
\title{NEWS file for the MicSim package}
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
   \item Extend function "buildTransitionMatrix" to facilitate using transition rates which depend on covariates such as gender and marital status.
   \item Correct typing error in Details section of help page of function "buildTransitionMatrix".
   \item Fix a problem in function "micSim" occuring if the simulation is asked not to produce newborns (concerns argument "fertTr=c()").
  }
}