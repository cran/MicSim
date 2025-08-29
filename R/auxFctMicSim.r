#####################################################################################
#####################################################################################
### AUXILIARY FUNCTIONS                                                            ##
### - TO DEFINE MICROSIMULATION INPUT                                              ##
### - FOR COMPUTATION WITH DATES                                                   ##
### SZ, June 2022                                                                  ##
#####################################################################################
#####################################################################################
#'
#' builtStatesCodes
#'
#' Assign to all states and subStates numerical codes
#' Build Transition Matrix
#'
#' @usage builtStatesCodes(transitionMatrix)
#' @description The buildStatesCodes function converts a state transition matrix into a numerical coding scheme.
#' Each state in the transition matrix is divided into sub-states (if separated by ???/???) and each sub-state is mapped to a unique numerical code.
#'
#' @param transitionMatrix Is the row and column names represent states.
#' States can be compound, separated by ???/???.
#'
#' @returns A data frame containing the coding scheme for the states.
#' @keywords internal
#' @noRd
#'
builtStatesCodes <- function(transitionMatrix){
  
  transitionMatrixNum <- transitionMatrix
  allStates <- rownames(transitionMatrix)
  allStatesMatrix <- do.call(cbind,sapply(allStates, strsplit, "/"))
  absStates <- setdiff(colnames(transitionMatrix), rownames(transitionMatrix))
  codeList <- vector(length=nrow(allStatesMatrix)+1, mode="list")
  for(j in 1:nrow(allStatesMatrix)){
    usj <- unique(allStatesMatrix[j,])
    codeList[[j]] <- cbind(usj, 1:length(usj))
    codeList[[j]][,2] <- ifelse(nchar(codeList[[j]][,2])%in% 1,paste("0", codeList[[j]][,2], sep=""),nchar(codeList[[j]][,2]))
  }
  codeList[[j+1]] <- cbind(absStates, -c(1:length(absStates)))
  
  allCodes <- c()
  allCodesSep <- NULL
  for(i in 1:length(allStates)){
    st <- unlist(strsplit(allStates[i], "/"))
    stNum <- c()
    for(k in 1:length(st)){
      subCode <- codeList[[k]][codeList[[k]][,1] %in% st[k],2]
      stNum <- c(stNum, subCode)
    }
    allCodes <- c(allCodes, paste(stNum, collapse = ""))
    allCodesSep <- rbind(allCodesSep, as.numeric(stNum))
  }
  codesSepAbs <- matrix(rep(NA, length(absStates)*ncol(allCodesSep)), ncol=ncol(allCodesSep), nrow=length(absStates))
  codesSepAbs[,1] <- -c(1:length(absStates))
  codingScheme <- data.frame(allStates = c(allStates,absStates),
                             allCodes= as.numeric(c(allCodes, -c(1:length(absStates)))),
                             allCodesSep=rbind(allCodesSep,codesSepAbs))
  return(codingScheme=codingScheme)
}
#'
#' rate_cS
#'
#' Construct matrix noting dependencies of functions
#' Build dependency matrix
#'
#' @usage rate_cS(allTr)
#' @description Helper function for MicSim. 
#' A function for building a matrix containing the arguments of the transition rate functions.
#'
#' @param allTr a matrix containing names of all transition rate functions.
#'
#' @returns A dependency matrix
#' @keywords internal
#' @noRd
#'
rate_cS <- function(allTr){
  rates <-  unique(allTr)
  form <- unique(unlist(lapply(rates, function(f) names(formals(f)))))
  depMatrix <- matrix(0, nrow = length(rates), ncol = length(form))
  colnames(depMatrix) <- form
  rownames(depMatrix) <- rates
  for(i in 1:nrow(depMatrix)){
    if(all(names(formals(rates[i])) %in% colnames(depMatrix))) {
      args <- names(formals(rates[i]))
      for(k in 1:length(args)) {
        depMatrix[i, match(args[k], colnames(depMatrix))] <- 1
      }
    }
  }
  depMatrix <- as.matrix(depMatrix)
  return(depMatrix)
}

