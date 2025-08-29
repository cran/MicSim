#' buildTransitionMatrix
#'
#' Determining transition pattern and transition functions
#' Construct matrix indicating transition pattern and naming the corresponding transition rate functions.
#'
#' @usage buildTransitionMatrix(allTransitions, absTransitions, stateSpace)
#' @description The function \code{buildTransitionMatrix} supports the constructing of the "transition matrix", which determines the transition pattern of the microsimulation model. The actual microsimulation is performed by \link{micSim} (sequentially) or by \link{micSimParallel} (parallel computing).
#'
#' @param allTransitions A matrix comprising all possible transitions between values of state variables in the first column and in the second column the names of the functions defining the corresponding transition rates.
#' @param absTransitions A matrix comprising the names of the absorbing states which individuals are always exposed to (such as "dead" and emigrated labeled as "rest") in the first column and in the second column the names of the functions defining the corresponding transition rates.
#' @param stateSpace A matrix comprising all non absorbing states considered during simulation.
#'
#'
#' @details The function \code{buildTransitionMatrix} is an auxiliary function for building the transition matrix required to run the microsimulation using \link{micSim} or \link{micSimParallel}.
#' \itemize{
#'  \item In \code{stateSpace} all state variables considered during simulation including their values have to be defined. Values are always described using labels.
#'  For example, label "M" for being married. Each column of \code{stateSpace} refers to one state variable considered and each row refers to one state of the state space.
#'  Apart from "m" and "f" reserved for male and female (state variable: gender) and "no" and "low" reserved for no education
#'  and elementary school attended (state variable: educational attainment), labels can be set arbitrarily.
#'  \item Each element of the first column of \code{allTransitions} has to be of the form "A->B" with indicating "A" the starting value of a transition
#'  and "B" the arrival value. ("->" is the placeholder defined to mark a transition).
#'  For example, "0" (childless) describes the starting value of the transition marking a first birth event and "1" (first child) its arrival value.
#'  All value labels used have to be identical to the value labels of the state variables specifying the simulation model.
#'  \item All absorbing states listed in the first column of \code{absTransitions} have to be given as strings such as "dead" for being dead or "rest" for emigrated.
#'  Since dying is a competing risk all individuals are always exposed to, "dead" is a mandatory part of \code{absTransitions}.
#'  \item All transitions can be defined to depend on several state variables.
#'  For example, a divorce rate depends on gender and on the fertility status.
#'  Therefore, the starting value and the arrival value of a transition have to be specified as a combination of the considered attributes,
#'  separated by a forward slash and in accordance with the ordering of the state variables in the state space.
#'  For example, "f/A->f/B" describes a female specific transition from "A" to "B" and "f/M/1 -> f/D/1" might describe a mothers (indicated by "1") transition from "M" (e.g., married) to "D" (e.g., divorced).
#'  \item For absorbing states, a prefix indicates the attributes on which a transition is assumed to depend (also separated by forward slashs),
#'  e.g., "f/dead" and "m/dead" describe gender specific mortality transitions and "f/M/dead" and "m/M/dead" indicate gender specific mortality rates for married persons.
#' }
#'
#' @returns A Transition Matrix. The \code{transitionMatrix} that is mandatory to perform a microsimulation run by \link{micSim} (sequentially) or by \link{micSimParallel} (parallel computing) is returned.
#' The matrix has as many rows as the simulation model comprises non-absorbing states and as many columns as the simulation model comprises absorbing and non-absorbing states.
#' The rows indicate starting states of transitions and the columns signify arrival states.
#' At positions indicating impossible transitions, the matrix contains zeros. Otherwise the name of the function defining the respective transition rates is given.
#'
#' @examples
#'
#' \dontrun{
#' # ########################################################################################
#' # 1. Example: Transition rates are specified to depend on only one state variable
#' # ########################################################################################
#' # Definition of state space, i.e., non-absorbing and absorbing states
#'
#' sex <- c("m","f")
#' fert <- c("0","1","2","3+")
#' marital <- c("NM","M","D","W")
#' edu <- c("no","low","med","high")
#' stateSpace <- expand.grid(sex=sex,fert=fert,marital=marital,edu=edu)
#'
#' # Possible transitions indicating fertility behavior are "0->1", "1->2", "2->3+", and "3+->3+".
#' # Here, "->" is the defined placeholder defining a transition.
#' # The "fert1Rates" marks the name of the function defining the transition rates to parity one and
#' # the "fert2Rates" marks the name of the function defining the transition rates to higher parities.
#' #
#' # Note:
#' # The functions "fert1Rates" and "fert2Rates" are transition rate functions defined by the user.
#' # Their naming depends on the users choice.
#'
#' fertTrMatrix <- cbind(c("0->1","1->2","2->3+","3+->3+"),
#'                c("fert1Rates", "fert2Rates", "fert2Rates","fert2Rates"))
#'
#' # Possible transitions indicating changes in the marital status are "NM->M", "M->D", "M->W", "D->M",
#' # and "W->M".
#' #
#' # The "marriage1Rates" marks the name of the function defining the transition rates
#' # for first marriage
#' # and "marriage2Rates" marks the name of the function defining the transition rates
#' # for further marriage
#' # The "divorceRates" marks the name of the function defining divorce rates
#' # and "widowhoodRates" marks the name of the function describing transition rates to widowhood.
#' #
#' # Note:
#' # The functions "marriage1Rates","marriage2Rates", "divorceRates", and "widowhoodRates" are
#' # transition rate functions defined by the user.
#' # Their naming depends on the users choice.
#'
#' maritalTrMatrix <- cbind(c("NM->M","M->D","M->W","D->M","W->M"),
#'                   c("marriage1Rates","divorceRates","widowhoodRates","marriage2Rates",
#'                     "marriage2Rates"))
#'
#' # Possible transitions indicating changes in the educational attainment
#' # are "no->low", "low->med", and "med->high".
#' # The "noToLowEduRates" marks transition rates for accessing primary education.
#' # The "noToLowEduRates" marks transition rates for graduating with a lower secondary education.
#' # And "medToHighEduRates" marks transition rates for graduating with a higher secondary education.
#' #
#' # Note:
#' # The functions "noToLowEduRates","noToLowEduRates", and "medToHighEduRates"
#' # are transition rate functions defined by the user.
#' # Their naming depends on the users choice.
#'
#' eduTrMatrix <- cbind(c("no->low","low->med","med->high"),
#'                 c("noToLowEduRates","noToLowEduRates","medToHighEduRates"))
#'
#' # Combine all possible transitions and the related transition function into one matrix.
#' # allTransitions <- rbind(fertTrMatrix, maritalTrMatrix, eduTrMatrix)
#' #
#' # Possible absorbing states are "dead" and "rest".
#' # (The latter indicates leaving the population because of emigration).
#' # The accordant transition rate functions are named "mortRates" and "emigrRates".
#' # (Again, naming is up to the user.)
#'
#' absTransitions <- rbind(c("dead","mortRates"),c("rest","emigrRates"))
#'
#' # Construct "transition matrix".
#'
#' transitionMatrix <- buildTransitionMatrix(allTransitions,absTransitions,stateSpace)
#'
#' # ########################################################################################
#' #  2. Example: Transition rates are gender specific
#' # ########################################################################################
#' # Definition of non-absorbing and absorbing states
#'
#' sex <- c("m","f")
#' stateX <- c("H","P")
#' stateSpace <- expand.grid(sex=sex,stateX=stateX)
#' absStates <- c("dead")
#'
#' # Transitions indicating changes in "stateX".
#' # We assume distinct transition rates for females and males.
#' #
#' # Note:
#' # The functions "ratesHP_f","ratesHP_m", "ratesPH_f", and "ratesPH_m"
#' # are transition rate functions defined by the user.
#'
#' trMatrix_f <- cbind(c("f/H->f/P","f/P->f/H"),c("ratesHP_f", "ratesPH_f"))
#' trMatrix_m <- cbind(c("m/H->m/P","m/P->m/H"),c("ratesHP_m", "ratesPH_m"))
#' allTransitions <- rbind(trMatrix_f,trMatrix_m)
#'
#' # We assume gender specific mortality rates.
#' # Note:
#' # The naming and specification of the respective mortality rate functions
#' # the "mortRates_f" and "mortRates_m" depend on the user.
#'
#' absTransitions <- rbind(c("f/dead","mortRates_f"), c("m/dead","mortRates_m"))
#'
#' transitionMatrix <- buildTransitionMatrix(allTransitions=allTransitions,
#'                      absTransitions=absTransitions, stateSpace=stateSpace)
#' }
#' @export
#'
#'
buildTransitionMatrix <- function(allTransitions,absTransitions,stateSpace){
  if(is.vector(allTransitions))
    allTransitions <- matrix(allTransitions, ncol=2, nrow=1)
  if(is.vector(absTransitions))
    absTransitions <- matrix(absTransitions, ncol=2, nrow=1)
  absStates <- absTransitions[,1]
  if(is.null(dim(stateSpace)))
    stateSpace <- matrix(stateSpace, ncol=1)
  absStNam <- c('dead')
  if('rest' %in% unlist(strsplit(absStates,'/')))
    absStNam <- c(absStNam, 'rest')
  transitionMatrix <- matrix(0,nrow=dim(stateSpace)[1], ncol=dim(stateSpace)[1]+length(absStNam))
  colnames(transitionMatrix) <- c(apply(stateSpace,1,paste,collapse='/'),absStNam)
  rownames(transitionMatrix) <- apply(stateSpace,1,paste,collapse='/')
  isInThisState <- function(ss,state){
    if(sum(ss %in% as.character(unlist(state)))==length(ss))
      return(TRUE)
    return(FALSE)
  }
  for(i in 1:length(absStates)){
    strAb <- unlist(strsplit(absStates[i],split='/'))
    if(length(strAb)==1){
      ia <- which(colnames(transitionMatrix)==absStates[i])
      transitionMatrix[,ia] <- absTransitions[i,2]
    } else {
      iAB <- which(strAb %in% c('dead','rest'))
      aS <- strAb[iAB]
      strAbCov <- strAb[-iAB]
      rA <- which(apply(stateSpace,1,isInThisState, ss=strAbCov)==TRUE)
      ia <- which(colnames(transitionMatrix)==aS)
      transitionMatrix[rA,ia] <- absTransitions[i,2]
    }
  }
  if(!is.null(allTransitions)){
    tr <- do.call(rbind,strsplit(allTransitions[,1],'->'))
    for(i in 1: dim(tr)[1]){
      trI <- tr[i,]
      oSPr <- unlist(strsplit(trI[1], split='/'))
      dSPr <- unlist(strsplit(trI[2], split='/'))
      idOS <- apply(stateSpace,1,isInThisState, ss=oSPr)
      idDS <- apply(stateSpace,1,isInThisState, ss=dSPr)
      stateSpaceOS <- stateSpace[idOS,,drop=F]
      stateSpaceDS <- stateSpace[idDS,,drop=F]
      for(j in 1:dim(stateSpaceOS)[1]){
        oS <- as.character(unlist(stateSpaceOS[j,]))
        for(k in 1:dim(stateSpaceDS)[1]){
          dS <- as.character(unlist(stateSpaceDS[k,]))
          c1 <- oS[!oS %in% oSPr]
          c2 <- dS[!dS %in% dSPr]
          if(sum(!(c1 %in% c2))==0 & sum(!(c2 %in% c1))==0){
            ir <- which(rownames(transitionMatrix)==paste(oS, collapse='/'))
            ic <- which(colnames(transitionMatrix)==paste(dS, collapse='/'))
            transitionMatrix[ir,ic] <- allTransitions[i,2]
          }
        }
      }
    }
  }
  return(transitionMatrix)
}

