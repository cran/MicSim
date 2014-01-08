####################################################################################
####################################################################################
## AUXILIARY FUNCTIONS TO DEFINE MICROSIMULATION INPUT                            ##
## SZ, November 2013                                                              ##
####################################################################################
####################################################################################

# Set simulation horizon
setSimHorizon <- function(startDate, endDate){
  dts <- c(startDate,endDate)
  simHorizon <- chron(dts,format=c(dates="d/m/Y"),out.format=c(dates="d/m/year"))
  return(simHorizon)
}

# Construct matrix indicating transition pattern and naming the corresponding transition rate functions.
buildTransitionMatrix <- function(allTransitions,absTransitions,stateSpace){  
  if(is.vector(absTransitions))
    absTransitions <- matrix(absTransitions, ncol=2, nrow=1)  
  absStates <- absTransitions[,1]  
  if(is.null(dim(stateSpace)))
    stateSpace <- matrix(stateSpace, ncol=1)  
  transitionMatrix <- matrix(0,nrow=dim(stateSpace)[1], ncol=dim(stateSpace)[1]+length(absStates))
  colnames(transitionMatrix) <- c(apply(stateSpace,1,paste,collapse="/"),absStates)
  rownames(transitionMatrix) <- apply(stateSpace,1,paste,collapse="/")
  for(i in 1:length(absStates)){
    ia <- which(colnames(transitionMatrix)==absStates[i])
    transitionMatrix[,ia] <- absTransitions[i,2]
  }  
  if(!is.null(allTransitions)){
    tr <- do.call(rbind,strsplit(allTransitions[,1],"->"))
    for(i in 1: dim(stateSpace)[1]){
      nn <- stateSpace[i,]
      for(j in 1:length(nn)){
        inff <- which(tr[,1] %in% as.character(unlist(nn[j])))
        dSS <- tr[inff,2,drop=F]
        dFF <- allTransitions[inff,2,drop=F]
        if(dim(dSS)[1]>0){
          nc <- nn
          for(k in 1:dim(dSS)[1]){
            nc[j] <- dSS[k]
            d <- which(colnames(transitionMatrix) %in% paste(sapply(nc, as.character), collapse="/"))
            transitionMatrix[i,d] <- dFF[k]
            nc <- nn
          }
        }
      }
    }
  }
  return(transitionMatrix)
}