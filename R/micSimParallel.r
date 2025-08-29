# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# II. Execute microsimulation distributed (by executing as many single thread microsimulations in parallel as cores
#     are available)
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
#'
#' micSimParallel
#'
#' @usage micSimParallel(initPop=NULL, immigrPop=NULL, initPopList= c(),
#'                      immigrPopList= c(), transitionMatrix, absStates=NULL,
#'                      varInitStates=c(), initStatesProb=c(), fixInitStates = c(),
#'                      maxAge=99, simHorizon, fertTr=c(), monthSchoolEnrol=c(),
#'                      cores=1, seeds=1254)
#'
#' @description Performs a continuous-time microsimulation run (sequentially, i.e., using only one CPU core).
#'
#' @param initPop Data frame comprising the starting population of the simulation.
#' @param immigrPop Data frame comprising information about the immigrants entering the population across simulation time.
#' @param initPopList A list of matrices, where each matrix represents a subset of initial population.
#' @param immigrPopList A list of matrices, where each matrix represents a subset of the immigrant population.
#' @param transitionMatrix A matrix indicating the transition pattern and the names of the functions determining the respective transition rates (with rates to be returned as vectors,
#' i.e. for input age 0 to 10 eleven rate values have to be returned).
#' @param absStates A vector indicating the absorbing states of the model.
#' @param fixInitStates (Vector of) Indices of SubStates determining the attributes/subStates that a newborn will be taken over from the mother. If empty or not defined, no attributes will be inherited.
#' @param varInitStates (A vector comprising the) SubStates / attributes that are assigned to a newborn randomly according to the probabilities \code{initStatesProb}, i.e. that are not inherited from the mother.
#' @param initStatesProb A vector comprising the probabilities corresponding to \code{varInitStates}.
#' If \code{fixInitStates} are given (i.e. attributes from the mother are inherited), these probabilities have to sum to one conditioned on the inherited attributes,
#' i.e. for each (set of) inherited attribute(s) separately. Otherwise, the sum of \code{initStatesProb} has to be one.
#' @param maxAge A scalar indicating the exact maximal age (i.e., sharp 100.00 years) which an individual can reach during simulation. \code{maxAge} has to be greater than zero.
#' @param simHorizon A vector comprising the starting and ending date of the simulation. Both dates have to be given in the format 'yyyymmdd'. The starting date has to precede the ending date.
#' @param fertTr A vector indicating all transitions triggering a child birth event during simulation, that is, the creation of a new individual.
#' @param monthSchoolEnrol The month (as numeric value from 1 to 12) indicating the general enrollment month for elementary school, e.g., 9 for September.
#' If transition to elementary school is not defined (see below under 'details') and no such month is given school enrollment to elementary school is not modeled / simulated.
#' @param cores Number of cores as the user specifies. Note: The number is restricted by the maximum number of cores a computer (or a computer cluster) has.
#' @param seeds To ensure that the results are replicable and are therefore reasonable, the user should always assign a seed to each PRNG (pseudo random number generator) representation used.
#' @details
#' The \code{micSimParallel} function is designed to perform population simulations in parallel using multiple processing cores. 
#' This function is particularly useful for large-scale simulations where computational efficiency is crucial.
#' The function takes an initial population and an immigrant population, both of which can be provided as single matrices or as lists of matrices for parallel processing. 
#' If the populations are provided as lists, each element of the list is processed by a different core, allowing for efficient distribution of the computational load.
#'
#' @returns A matrix representing the combined population after running the simulation in parallel.
#' Each row corresponds to an individual, and columns represent different attributes or states of the individuals.
#' The function ensures that unique IDs are assigned to newborns and that the final output is ordered correctly by these IDs.
#' This matrix includes the results of the simulations for all subsets of the initial and immigrant populations processed across different cores.
#' 
#' @export
#'
#'
micSimParallel <- function(initPop=NULL, immigrPop=NULL, initPopList = c(), immigrPopList = c(),
                           transitionMatrix, absStates=NULL, varInitStates=c(), initStatesProb=c(),
                           fixInitStates = c(), maxAge=99, simHorizon, fertTr=c(), monthSchoolEnrol=c(),
                           cores=1, seeds=1254){
  
  cat('Starting at '); print(Sys.time())
  if(!is.null(initPop)) {
    N <- dim(initPop)[1]
  } else {
    if(length(initPopList)>0){
      N <- 0
      for(k in 1:cores){
        N <- N +nrow(initPopList[[k]])
      }
    } else {
      stop("No initial population for parallel computing has been given.\n")
    }
  }
  
  if(!is.null(immigrPop)) {
    M <- dim(immigrPop)[1]
  } else {
    M <- 0
    if(length(immigrPopList)>0){
      for(k in 1:cores){
        M <- M +nrow(immigrPopList[[k]])
      }
    }
  }
  
  # Split starting population and (if available) immigrant population according to available cores
  if(length(cores) %in% 0)
    stop("At least one core must be given.\n")
  
  condSplit <- ((length(initPopList) %in% cores) & is.null(immigrPop)) |
    ((length(initPopList) %in% cores) & (!is.null(immigrPop) & (length(immigrPopList) %in% cores)))
  
  if(condSplit){
    cat('\nAssign cases to distinct cores according to the split provided.\n')
    cat('Beware: It is not checked whether cases appear twice in the splits.\n')
    cat('If duplicates are in the different splits, this will result in duplicate life histories for the same entities.\n')
    cat('Thus, please check for duplicate IDs in advance.\n')
  }
  
  if(!condSplit) {
    
    if(length(initPopList)>0 & !(length(initPopList) %in% cores)){
      cat('\nSplit of initial population given for parallel computing does not match the number of cores determined.\n')
      cat('Therefore, MicSim makes an automated assignment of cases of the initial population to the distinct cores.\n')
      cat('At this, cases are distributed to the cores such that at each core approx. the same number of cases is simulated.\n')
    }
    if(!is.null(immigrPop) & (length(immigrPopList)>0 & !(length(immigrPopList) %in% cores))){
      cat('\nSplit of immigrant population given for parallel computing does not match the number of cores determined.\n')
      cat('Therefore, MicSim makes an automated assignment of cases of the immigrant population to the distinct cores.\n')
      cat('At this, cases are distributed to the cores such that at each core approx. the same number of immigrant cases is simulated.\n')
    }
    
    widthV <- max(trunc(N/cores), 10)
    widthW <- max(trunc(M/cores), 10)
    intV <- matrix(NA,ncol=2,nrow=cores)
    intW <- matrix(NA,ncol=2,nrow=cores)
    nI <- trunc(N/widthV)
    nIM <- trunc(M/widthW)
    ni <- 1
    for(i in 1:(nI-1)){
      intV[i,1] <- ni
      intV[i,2] <- ni+widthV-1
      ni <- ni+widthV
    }
    intV[nI,1] <- ni
    intV[nI,2] <- N
    ni <- 1
    if(nIM>1){
      for(i in 1:(nIM-1)){
        intW[i,1] <- ni
        intW[i,2] <- ni+widthW-1
        ni <- ni+widthW
      }
    }
    intW[nIM,1] <- ni
    intW[nIM,2] <- M
    initPopList <- list()
    immigrPopList <- list()
    for(core in 1:cores){
      if(!is.na(intV[core,1])){
        initPopList[[core]] <- initPop[intV[core,1]:intV[core,2],]
      } else {
        initPopList[[core]] <- NA
      }
      if(!is.na(intW[core,1])){
        immigrPopList[[core]] <- immigrPop[intW[core,1]:intW[core,2],]
      } else {
        immigrPopList[[core]] <- NA
      }
    }
  }
  
  nL <- cores - sum(unlist((lapply(initPopList, is.na))))
  mL <- cores - sum(unlist((lapply(immigrPopList, is.na))))
  
  sfInit(parallel=T,cpus=cores,slaveOutfile=NULL)
  sfExportAll(debug=FALSE)
#  sfLibrary(MicSim)
  sfClusterSetupRNGstream(seed=(rep(seeds,35)[1:length(cores)]))
  myPar <- function(itt){
    #cat('Starting thread: ',itt,'\n')
    if(itt<=mL){
      immigrPopL <- immigrPopList[[itt]]
    } else {
      immigrPopL <- NULL
    }
    if (itt<=nL) {
      initPopL <- initPopList[[itt]]
    } else {
      initPopL <- NULL
      stop("\nCompared to the number of migrants, the starting population is too small to justify running a distributed simulation on several cores.")
    }
    popIt <- micSim(initPop=initPopL, immigrPop=immigrPopL, transitionMatrix=transitionMatrix,
                    absStates=absStates, varInitStates=varInitStates, initStatesProb=initStatesProb,
                    fixInitStates=fixInitStates, maxAge=maxAge, simHorizon=simHorizon, fertTr=fertTr,
                    monthSchoolEnrol=monthSchoolEnrol)
    #cat('Thread: ',itt,' has stopped.\n')
    return(popIt)
  }
  pop <- sfLapply(1:max(nL,mL), myPar)
  # create unique IDs for newborns
  refID <- 0
  replaceID <- function(rr){
    pop[[i]][which(as.numeric(pop[[i]][,1])==rr[1]),1] <<- rr[2]
    return(NULL)
  }
  for(i in 1:length(pop)){
    if(!is.na(immigrPopList[[i]])[1]){
      allIDs <- c(initPopList[[i]]$ID, immigrPopList[[i]]$ID)
    } else {
      allIDs <- initPopList[[i]]$ID
    }
    exIDs <- unique(as.numeric(pop[[i]][,1]))
    repl <- setdiff(exIDs, allIDs)
    if(length(repl)>0) {
      newIDs <- cbind(repl,-(refID+(1:length(repl))))
      idch <- apply(newIDs,1,replaceID)
      refID <- max(abs(newIDs[,2]))
    }
  }
  pop <- do.call(rbind,pop)
  pop[as.numeric(pop[,1])<0,1]  <- abs(as.numeric(pop[as.numeric(pop[,1])<0,1]))+N+M
  pop <- pop[order(as.numeric(pop[,1])),]
  sfStop()
  cat('Stopped at '); print(Sys.time())
  return(pop)
}


