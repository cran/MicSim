#####################################################################################
#####################################################################################
### FUNCTION TO CONVERT MICROSIMULATION OUTPUT INTO WIDE FORMAT                    ##
### SZ, April 2014                                                                 ##
#####################################################################################
#####################################################################################
#'
#' convertToWideFormat
#'
#' Reshaping microsimulation output into wide format
#'
#' @usage convertToWideFormat(pop)
#' @description The function reshapes the output given by \link{micSim} or by \link{micSimParallel} into wide format. In wide format, the data comprises for each episode which an individual experiences additional column entries.
#'
#' @param pop The data frame \code{pop} contains the whole synthetic population considered during simulation including all events generated. For each individual \code{pop} contains as many rows as the individual performed transitions during simulation.
#'
#' @returns  A data frame comprising the microsimulation output in wide format.
#' \itemize{
#'  \item \code{ID} is the unique numerical person identifier of an individual.
#'  \item \code{birthDate} is the birth date of an individual.
#'  \item \code{initState} is the state in which an individual initially entered the virtual population of the simulation.
#'  \item \code{ns} gives the number of (completed) episodes an individual has passed.
#'  \item  The variables \code{From.i} and \code{To.i} mark the start and the arrival state of the transition corresponding to episode \code{i}. The variables \code{transitionTime.i} and \code{transitionAge.i} give the corresponding transition time and age. The enumerator \code{i} ranges from 1 to the maximal number of transitions which an individual experienced during simulation. Only completed episodes are counted.
#' }
#' @examples
#'
#' \dontrun{
#' # Run microsimulation before, e.g., the complex example
#' # described on the help page of the function "micSim".
#'
#' pop <- micSim(initPop, immigrPop, transitionMatrix, absStates, initStates,
#'     initStatesProb, maxAge, simHorizon, fertTr)
#' popWide <- convertToWideFormat(pop)
#' }
#' @export
#'
#'
convertToWideFormat <- function(pop){
  
  giveSeq <- function(nu){
    return(1:nu)
  }
  
  if("motherID" %in% colnames(pop))
    pop <- pop[,!(colnames(pop) %in% "motherID")]
  
  popTemp <- pop
  ns <-  data.frame(table(popTemp$ID),stringsAsFactors=FALSE)
  ns <- ns[order(as.numeric(as.character(ns[,1]))),]
  colnames(ns) <- c('ID','ns')
  popTemp <- merge(popTemp, ns, by='ID')
  popTemp <- popTemp[order(as.numeric(popTemp[,c('ID')])),]
  nsU <- popTemp$ns[which(!duplicated(popTemp[,'ID']))]
  popTemp$Episode <- unlist(sapply(nsU,giveSeq))
  popTemp <- popTemp[,c('ID','birthDate','initState','ns','Episode','From','To','transitionTime','transitionAge')]
  popWide <- reshape(popTemp, timevar = 'Episode', idvar = 'ID', direction = 'wide',
                     v.names=c('From','To', 'transitionTime', 'transitionAge'))
  popWide$ns[which(is.na(popWide$transitionTime.1))] <- 0
  return(popWide)
}


