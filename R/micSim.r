####################################################################################
####################################################################################
## FUNCTION EXECUTING MICROSIMULATION                                             ##
## SZ, June 2022                                                                  ##
####################################################################################
####################################################################################
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# I. Execute microsimulation as single thread
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
micSim <- function(initPop, immigrPop=NULL, transitionMatrix, absStates=NULL, fixInitStates = c(), 
                   varInitStates=c(), initStatesProb=c(), maxAge=99, simHorizon, fertTr=c(), monthSchoolEnrol=c()) {
  
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # A. CHECK INPUT FOR CONSISTENCY
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------                   
  if(is.null(initPop))
    stop('No starting population has been defined.')
  if(!is.null(initPop)){
    if(paste(colnames(initPop),collapse='/')!='ID/birthDate/initState')
      stop('Matrix specifying the starting population has not been defined properly.')
  }
  if(!is.null(immigrPop)){
    if(paste(colnames(immigrPop),collapse='/')!='ID/immigrDate/birthDate/immigrInitState')
      stop('Matrix specifying immigrants has not been defined properly.')
  }
  if(is.null(transitionMatrix))
    stop('Matrix defining transition pattern und functions has not been defined properly.')
  if(maxAge<=0)
    stop('The maximal age until which individual life courses are simulated should exceed zero.')
  if(length(simHorizon)!=2)
    stop('The simulation horizon has not been defined properly.')
  if(is.null(absStates))
    absStates <- setdiff(colnames(transitionMatrix),rownames(transitionMatrix))  
  if(length(fertTr)>0){
    if((is.null(varInitStates) & is.null(initStatesProb)))
      stop('For children potentially born during simulation no inital state(s) and/or corresponding occurrence probabilities have been defined.')
    if(length(fixInitStates)>0){
      for(i in 1:length(fixInitStates)){
        ssum <- sum(initStatesProb[apply(varInitStates,1, function(rr){varInitStates[,fixInitStates[i]][1] %in% rr})])
        if(ssum!=1)
          stop('The sum of the probabilities to assign initial states to newborns must equal 1.')
      } 
    } else {
      if(sum(initStatesProb)!=1)
        stop('The sum of the probabilities to assign initial states to newborns must equal 1.')
    }
  }
  # check whether all functions delivering transition rates deliver vectors of rates as output (necessary for integration procedure later on)
  allTr <- unique(as.vector(transitionMatrix)[as.vector(transitionMatrix) !="0"])
  simStartInDays <- getInDays(simHorizon[1])
  simStopInDays  <- getInDays(simHorizon[2]) 
  ranYear <- c(getYear(simStartInDays), getYear(simStopInDays)) 
  if(length(fertTr)>0){
    minAge <- 0
  } else {
    minAge <- min(trunc(getAgeInDays(simHorizon[1], initPop$birthDate)/365.25))
  }
  ranAge <- c(minAge,maxAge)
  ran <- min(c(diff(ranYear), diff(ranAge)))
  durIn <- 0:ran
  for(tr in 1:length(allTr)){ # check whether all functions delivering transition rates deliver vectors of rates as output (necessary for integration procedure later on)
    for(cal in c(getYear(simStartInDays): getYear(simStopInDays))){
      for(age in minAge:(maxAge-1)){
        for(dur in 0:ran){
          res <- eval(do.call(allTr[tr], args=list(age=age,calTime= cal,duration=dur)))  
          if(anyNA(res)){
            cat("The rates function for ", allTr[i], " does not deliver a vector of rates for an input vector of age, calendar time, and/or duration (all in years).\n")
            cat("The missing rate occurs at year ",cal, " for age ", age, " and duration ", dur, "\n.")
            cat("This is a requirement for the simulation procedure to run since it is based on integrated hazard rates.\n")
            stop('Incorrect definition of input rates function!')
          }
        }
      }
    }
  }
  if(length(monthSchoolEnrol)==0){
    schoolEnrol <- FALSE
  } else {
    schoolEnrol <- TRUE 
  }
  
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------
  # B. DEFINITION OF GLOBAL PARAMETERS
  # --------------------------------------------------------------------------------------------------------------------
  # --------------------------------------------------------------------------------------------------------------------   
  
  # Simulation horizon
  simStartInDays <- getInDays(simHorizon[1])
  simStopInDays  <- getInDays(simHorizon[2])
  # Event queue
  queue <- matrix(NA,ncol=6,nrow=0)
  colnames(queue) <- c('ID','currTime','currState','currAge','nextState','timeToNextState')
  # Global time
  t.clock <- simStartInDays  # counts in days since 01-01-1970
  # Recording transitions performed
  transitions <- matrix(NA,ncol=5,nrow=0)
  colnames(transitions) <- c('ID', 'From', 'To', 'transitionTime', 'transitionAge') 
  # Record linkage of mothers to newborns via their IDs
  if(length(fertTr)>0){
    mothers <- matrix(NA,ncol=2,nrow=0) # 'motherID', 'childID' 
  }
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # C. FUNCTIONS REQUIRED FOR SIMULATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  
  # Function building matrix indicating the transitions between states causing a newborn
  buildFertTrExpanded <- function(){
    allStates <- matrix(rownames(transitionMatrix), ncol=1)
    allStatesSplit <- apply(allStates,1,strsplit, split="/")
    fert <- do.call(rbind,(strsplit(fertTr,split='->'))) 
    fertTrExpanded <- NULL
    for(i in 1:nrow(allStates)){
      cS <- allStatesSplit[[i]][[1]]
      for(j in 1:nrow(allStates)){
        dS <- allStatesSplit[[j]][[1]]
        if(("f" %in% cS) & ("f" %in% dS)){
          for(k in 1:nrow(fert)){
            ff <- fert[k,]
            oS <- strsplit(ff[1],'/')[[1]]
            bS <- strsplit(ff[2],'/')[[1]]  
            cond1 <- !(F %in% (oS %in% cS)) & !(F %in% (bS %in% dS))
            cond2 <- paste((cS[!(cS %in% oS)]),collapse="/") == paste((dS[!(dS %in% bS)]),collapse="/") # if there are a fertility event only one substate can change, namely that one belonging to the fertility attribute
            if(cond1 & cond2){
              fertTrExpanded <- rbind(fertTrExpanded, c(paste0(cS,collapse="/"), paste0(dS,collapse="/")))
            }
          }
        }
      }
    }
    return(fertTrExpanded)
  }
  
  # Function checks whether a transition causes a newborn. 
  # (Demands `fertTrExpanded': matrix indicating the transitions between states causing a newborn (defined by `fertTr').)
  isBirthEvent <- function(currState, destState){
    oS <- which(fertTrExpanded[,1] %in% currState)  
    dS <- which(fertTrExpanded[,2] %in% destState)    
    if(length(intersect(oS,dS))>0)
      return(TRUE)
    return(FALSE)
  }
  
  # Function adds to simulation population a newborn.
  addNewNewborn <- function(birthTime=birthTime, motherID=motherID, motherState=motherState){    
     
    birthTime=t.clock
    motherID=indS['ID']
    motherState=indS['nextState']
    
    if(length(fixInitStates)>0){
      motherStatePart <- strsplit(motherState, split="/")[[1]][fixInitStates]
      inSt <- which(apply(varInitStates[,fixInitStates, drop=F],1,paste0, collapse=",") %in% paste(motherStatePart, collapse=",")) 
      varInitStatesR <- varInitStates[inSt,,drop=F]
      initStatesProbR <- initStatesProb[inSt]
      probStatePart <- varInitStatesR[sample(1:nrow(varInitStatesR),size=1,replace=T,prob=initStatesProbR),]
    } else {
      probStatePart <- varInitStates[sample(1:nrow(varInitStates),size=1,replace=T,prob=initStatesProb),]
    }
    birthState <- paste0(probStatePart,collapse="/")
    if(is.null(immigrPop)){
      id <- as.numeric(max(as.numeric(initPop[,'ID'])))+1
    } else {
      id <- as.numeric(max(c(as.numeric(immigrPop[,'ID']),as.numeric(initPop[,'ID']))))+1
    }
    newInd <- c(id,getInDateFormat(birthTime),birthState)
    #cat('NewBorn: ',newInd,'\n')
    initPop <<- rbind(initPop,newInd)
    mothers <<- rbind(mothers, c(motherID, id))
    nE <- getNextStep(c(id,birthState,0,birthTime))
    #cat('\n------------n')
  }
  
  # Function checks whether a transition implies a school enrollment (in the year when child turns seven).
  # (If state 1 comprises value `no' and state 2 comprises value `low', the transition is marked as `school enrollment'.)
  isSchoolEnrolment <- function(currState,destState){ # TODO: speed up
    enrol <- all( c(any('no' %in% strsplit(currState,'/')[[1]]), any('low' %in% strsplit(destState,'/')[[1]])) ) 
    return(enrol)
  }
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # D. SIMULATION STEP
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Function to compute the next transition state and time of an individual who at age `currAge' and at time `calTime' 
  # entered its current state `currState'.
  getNextStep <- function(inp, isIMInitEvent=F){     
    # Transform input
    id <- as.numeric(unlist(inp[1], use.names = FALSE)) 
    currState <- as.character(unlist(inp[2], use.names = FALSE)) 
    currAge <- as.numeric(unlist(inp[3], use.names = FALSE)) # age in days 
    calTime <- as.numeric(unlist(inp[4], use.names = FALSE)) # calendar time in days since 01-01-1970
    # first event of an immigrant: he/she enters the population later than sim. starting time
    lagToWaitingTime <- ifelse(isIMInitEvent, (calTime - simStartInDays)/365.25,0) # in years
     #cat('\n-----\nID: ',id,'\n')
     #print(inp)
    ageInYears <- currAge/365.25 
     #cat('Age: ',ageInYears,' - CalTime: ',getYear(calTime),'-',getMonth(calTime),'-',getDay(calTime),'\n')
    # Possible destination states
    possTr <- transitionMatrix[match(currState, rownames(transitionMatrix)),]    
    possTr <- possTr[which(possTr !=0)]
    nextEventMatrix <- matrix(0, ncol=2, nrow=length(possTr))   

    # How many years (along age scale) remain until `maxAge'? 
    ranMaxAge <- (maxAge-0.01)-ageInYears     
    # How many years (along cal. time scale) remain until simulation end?        
    ranMaxYear <-  (simStopInDays - calTime)/365.25  
    # Identify the time range that should be considered. 
    ran <- min(ranMaxYear,ranMaxAge)   
    #cat('Ran: ',ran,' - ranMaxAge: ',ranMaxAge,' - ranMaxYear: ',ranMaxYear,'\n')
    #ranAge <- c(ageInYears,ageInYears+ranMaxAge) # age range in years
    #ranYear <- c(getYear(calTime), getYear(calTime)+ran) # year range in years
    #cat('RanAge: ',ranAge,' - ranYear: ',ranYear,'\n')
    
    # Extract transition history of individual until current cal. time.
    historiesInd <- transitions[as.numeric(transitions[,'ID']) %in% id & transitions[,'transitionTime'] <= calTime,,drop=F]
    iPIDout <- as.numeric(initPop[,'ID']) %in% id 
    initPopInd <- initPop[iPIDout,]
    birthTime <-  initPopInd['birthDate']
    initState <-  as.character(unlist(initPopInd['initState']))
    # Extract for each state the duration until transition (in days). 
    # Here, we have to differ between states of which we do not know when they are entered (i.e., `initial states' of members
    # of the starting population and the states of migrants when they entered the country), and 
    # the states we know the `entering date' as well as the `leaving date' (if the state has been left). 
    if(getInDays(birthTime) < simStartInDays | id %in% as.numeric(immigrPop[,'ID'])) {  
      dur <- rbind(c(initState,NA),historiesInd[,c('To','transitionTime')]) 
      dur <- cbind(dur,c(diff(as.numeric(dur[,2])),0))
      colnames(dur) <- c('TransitionTo','AtTime','durUntil')
      dur[which(is.na(dur[,'AtTime'])),'durUntil'] <- NA
    } else {  # Individual is born during simulation.
      birthTime <- getInDays(initPop[as.numeric(initPop[,'ID']) %in% id, 'birthDate']) # in days since 01-01-1970 
      dur <- rbind(c(initState,birthTime),historiesInd[,c('To','transitionTime')]) 
      dur <- cbind(dur,c(diff(as.numeric(dur[,2])),0))
      colnames(dur) <- c('TransitionTo','AtTime','durUntil')
    }
    # Compute for each possible destination state a waiting time.  
    for(i in 1:length(possTr)){
      tr <- possTr[i]
      destState <-  names(tr)
      cS <- strsplit(currState,'/')[[1]] # TODO: to speed up, this operation has to be improved/changed
      dS <- strsplit(destState, '/')[[1]] # TODO: to speed up, this operation has to be improved/changed
      # To determine the duration (time elapsed since last transition) that applies for the considered destination state,
      # we have to determine the duration since the last change in the covariate affected. 
      # For example, to specify the time being married, we have to determine the duration since (last) marriage. 
      covToCh <- which((cS==dS)==F)
      durSinceLastCovCh <- Inf  # For the transition to `dead' so far the time elapsed since the last transition does not play any role.  
      if(length(covToCh)==1){
        covHist <- do.call(rbind,sapply(dur[,'TransitionTo'],strsplit,split='/'))[,covToCh]
        idd <- which(covHist==cS[covToCh])
        if(length(idd)>1){
          if(F %in% (diff(idd)==1)){
            y <- rev(idd)[c(-1,diff(rev(idd)))==-1]
            idd <- rev(y)[c(diff(rev(y)),1)==1]
          }
        }
        durSinceLastCovCh <- sum(as.numeric(dur[idd,'durUntil'])) # If I do not know how long an individual already is in a state: This gives NA.
        if(is.na(durSinceLastCovCh))
          durSinceLastCovCh <- 365.25 # Then assume the individual is already for one year (=365.26 days) in that state.
      }  
      if(length(covToCh)>1 & (!destState %in% absStates)){
        cat('Recognized a possible transition implying a change of two or more covariates.',
            'Concerning the derivation of the time being elapsed since the last transition this feature is not yet implemented.', 
            'Current State: ',currState,' -> Possible transition to ',destState,'\n') 
      }
      tageInYears <- trunc(ageInYears)            
      tCalTime <- trunc(1970.001+calTime/365.25)  
      tdurSinceLastCovCh <- trunc(durSinceLastCovCh/365.25)
      indRateFctDET <- function(x){               
        res <- eval(do.call(tr,                   
                            args=list(age=tageInYears+x,calTime=tCalTime+x,duration=tdurSinceLastCovCh+x)))
        return(res)                               
      }                                                 
      ranAccuracyInDays <- (0:(trunc(ran*365.25)+0.99))/365.25
      detE <- indRateFctDET(ranAccuracyInDays)
      daysToTrInYears <- (which(detE == Inf)[1] - 1)/365.25
      if (Inf %in% detE) {
        timeToNext <- daysToTrInYears
      } else {   
        u <- -log(1-runif(1)) 
        #cat('It: ',i,'--u: ',u,'\n')
        # Extract individual transition rate (depending on age, calendar time, and time elapsed)  
        indRateFct <- function(x){
          ageIn <- ageInYears+x
          calIn <- 1970.001+calTime/365.25+x 
          durIn <- durSinceLastCovCh/365.25+x
          res <- eval(do.call(tr, args=list(age=ageIn,calTime= calIn,duration=durIn)))   
          if(TRUE %in% (res<0))
            stop('I have found negative rate value/s for transition: ',tr,'\n
                 This is implausible. Please check this. Simulation has been stopped.\n')
          #cat('x: ',x,' -- res', res,'\n')
          #cat('\n---\n')
          return(res)
        }
        if(sum(indRateFct(0:ran))==0){ # Rate funtion contains only zeros.
          intHaz <- 0
        } else {        
          # Integrated hazard at max. value
          intHaz <- try(integrate(indRateFct, lower=0, upper=ran)$value, silent=TRUE)
          if(inherits(intHaz, 'try-error')){          
            intHaz <- integrate(indRateFct, lower=0, upper=ran, stop.on.error = FALSE, rel.tol = 0.01)$value
          }
        }
        # If transformed random variate exceeds max. value of integr. hazard, we will not find a finite random waiting time.      
        if(u<=intHaz){
          invHazFct <- function(x){
            #cat('x: ',x,'\n')
            try.res <- try(integrate(indRateFct, lower=0, upper=x)$value-u, silent=TRUE)
            #print(try.res)
            if(inherits(try.res, 'try-error')){  
              #cat('Seemingly, divergent intergral for ID ',id,
              # ' in state ',currState,' at age ',currAge,' at time ',calTime, ' to state ',destState,
              #  ' for random number: ',u,'\n')  
              try.res <- integrate(indRateFct, lower=0, upper=x, stop.on.error = FALSE, rel.tol = 0.01)$value-u
            } 
            #cat('res: ',try.res,'\n-----\n')
            return(try.res)
          }  
          # Find random waiting time. 
          timeToNext <- uniroot(invHazFct,interval=c(0,ran))$root      
        } else {
          timeToNext <- Inf
        }
      }
      nextEventMatrix[i,1] <- destState    
      nextEventMatrix[i,2] <- (timeToNext+lagToWaitingTime)*365.25    # time to next event in days
    }
    #print(nextEventMatrix)
    nE <- nextEventMatrix[which(nextEventMatrix[,2]==min(as.numeric(nextEventMatrix[,2]))),,drop=F]
    if(dim(nE)[1]>1)
      nE <- nE[1,,drop=F]
    if(nE[1,2]!=Inf){
      # Cal. time of next event of individual. (If there is one.)
      tt <- calTime + as.numeric(nE[1,2]) 
      #print(tt)
      #cat(nE[1,1],'---',tt,'\n')  
      # Check whether next event implies school enrollment. If yes, adjust transition time to ensure that the individual
      # enters school at Sept. 1 in the year he/she turns seven.
      if(schoolEnrol){ # Is school enrollment considered in this simulation model? Yes, then continue; otherwise skip this part
        if(isSchoolEnrolment(currState,nE[1,1])){
          enYear <- getYear(tt)
          if(getMonth(tt) <= monthSchoolEnrol) {
            enDate <- getInDays_my(enYear, monthSchoolEnrol) 
          } else {
            enYear <- enYear+1
            enDate <- getInDays_my(enYear, monthSchoolEnrol) 
          }      
          diffToEn <- as.numeric(enDate-tt)
          nE[1,2] <- as.numeric(nE[1,2]) + diffToEn 
        }
      }
      # Enqueue new event (if there is one).
      queue <<- rbind(queue, c(id, t.clock, currState, currAge - lagToWaitingTime*365.25, nE[1,1], nE[1,2]))
    }  
    #cat('\n----------\n')      
    return(nE)
  }
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # E. INITIALIZATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Compute next events for members of starting population
  cat('Initialization ... \n')
  time_init_start = Sys.time()
  print(paste("Starting at: ", time_init_start))   

  IN <- data.frame(ID=initPop[,'ID'],currState=initPop[,'initState'],age= simStartInDays-getInDays(initPop[,'birthDate']),
                   calTime=rep(simStartInDays,dim(initPop)[1]),stringsAsFactors=FALSE) 
  if(TRUE %in% (IN$age<0)) {
    cat("There are persons born later than simulation starting date in the initial population. Related IDs are: ")
    negAge <- IN$ID[IN$age<0]
    for(i in 1:length(negAge)){
      cat(negAge[i]," ")
    }
    stop("Error: Negative age in initial population.")
  }
  if(TRUE %in% (IN$age/365.25>maxAge)) {
    cat("In the initial population, there are persons older than maxAge at simulation starting date. Related IDs are: ")
    invalAge <- IN$ID[IN$age/365.25>maxAge]
    for(i in 1:length(invalAge)){
      cat(invalAge[i]," ")
    }
    stop("Error: Older than max. age in initial population.")
  }
  init <- apply(IN, 1, getNextStep)

  # If immigrants enter the population, compute next events for them.
  if(!is.null(immigrPop)){
    IM <- data.frame(ID=immigrPop[,'ID'], currState=immigrPop[,'immigrInitState'], age=getAgeInDays(immigrPop[,'immigrDate'],immigrPop[,'birthDate']),
                     calTime=getInDays(immigrPop[,'immigrDate']),stringsAsFactors=FALSE) 
    # Check whether migrants are already born when they migrate
    if(TRUE %in% (immigrPop$immigrDate<immigrPop$birthDate)){
      cat("In the immigration population, there are persons who are not yet born when they immigrate. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate<immigrPop$birthDate]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }
      stop("Error: Not yet born at immigration date.") 
    }
    # Check whether all migrants migrate after simulation starting date
    if(TRUE %in% (immigrPop$immigrDate<simHorizon[1])){
      cat("In the immigration population, there are persons who are specified to migrate into the virtual population before simulation starting time. That's against MicSim's concept of migration. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate<simHorizon[1]]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }
      stop("Error: Migrate before simulation starting date.") 
    }
    # Check whether all migrants migrate before simulation stopping date  
    if(TRUE %in% (immigrPop$immigrDate>simHorizon[2])){
      cat("In the immigration population, there are persons who are specified to migrate into the virtual population after simulation ending time. That's a bit meaningless. Related IDs are: ")
      invalImDate <- immigrPop$ID[immigrPop$immigrDate>simHorizon[2]]
      for(i in 1:length(invalImDate)){
        cat(invalImDate[i]," ")
      }      
      stop("Error: Migrate after simulation ending date.") 
    }
    # Check whether all migrants are younger than maxAge when they migrate  
    ageIm <- (getInDays(immigrPop$immigrDate)-getInDays(immigrPop$birthDate))/365.25 # age at immigration in years
    if(TRUE %in% (ageIm>maxAge)){
      cat("In the immigration population, there are persons who older than `maxAge' when they migrate. That's a bit meaningless. Related IDs are: ")
      invalImAge <- immigrPop$ID[ageIm>maxAge]
      for(i in 1:length(invalImAge)){
        cat(invalImAge[i]," ")
      }  
      stop("Error: Migrants in the input data are older than `maxAge'.") 
    }    
    immigrInitPop <- immigrPop[,c('ID','birthDate','immigrInitState')]
    colnames(immigrInitPop)[3] <- 'initState'
    initPop <- rbind(initPop, immigrInitPop)
    imit <- apply(IM, 1, getNextStep, isIMInitEvent=T)
  } 
  if(length(fertTr)>0){
    fertTrExpanded <- buildFertTrExpanded()
  }
  
  time_init_end = Sys.time()
  print(paste("Ending at: ", time_init_end)) 
  #p_time = (time_init_end - time_init_start)
  #print(paste("#Time needed for initialization: ", p_time)) 
  
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # F. SIMULATION
  # ----------------------------------------------------------------------------------------------------------------------
  # ----------------------------------------------------------------------------------------------------------------------
  # Run simulation either until queue is empty or until simulation horizon has been reached.
  cat('Simulation is running ... \n')
  currYear <- trunc(simHorizon[1]/10000)
  cat('Year: ',currYear,'\n')
  while(dim(queue)[1]>0 & t.clock <= simStopInDays){ 
    
    # Sort queue according to soonest event to happen. 
    queue <- queue[order(as.numeric(queue[,'currTime']) + as.numeric(queue[,'timeToNextState'])),,drop=F] # In queue currTime in days since 01-01-1970
    #print_t.clock <- paste(c(getDay(t.clock), getMonth(t.clock), getYear(t.clock)), collapse="/")
    #print(print_t.clock)
    #print(dim(queue)[1])
    # Enqueue individual who has the soonest event to happen.
    indS <- queue[1,]  
    #print(indS)
    # Remove he/she from queue.
    queue <- queue[-1,,drop=F]
    # Set the global clock.
    t.clock <- as.numeric(indS['currTime']) + as.numeric(indS['timeToNextState']) # in days since 01-01-1970
    cY <- getYear(t.clock) # transform days since 01-01-1970 back to years
    # If the global clock exceeds the simulation horizon, stop simulation.
    if(t.clock > simStopInDays)
      break 
    if(cY>currYear){
      cat('Year: ',cY,'\n')
      currYear <- cY
    }   
    # Age at current transition  
    age <- as.numeric(indS['currAge']) +  as.numeric(indS['timeToNextState']) # in days
    # Register transition.
    transitions <- rbind(transitions, c(indS[c('ID','currState','nextState')], t.clock, age))
    # If current state is not an absorbent one, compute next event.
    if(!indS['nextState'] %in% absStates){
      # Current transition causes a newborn? If yes, add one to simulation population.
      if(length(fertTr)>0){ 
        if(isBirthEvent(indS['currState'],indS['nextState'])){
          addNewNewborn(birthTime=t.clock, motherID=indS['ID'], motherState=indS['nextState'])
        } 
      }
      res <- getNextStep(c(indS[c('ID','nextState')], age, t.clock))
      #print(res)
    }
    #cat('\n-----------\n')
  }
  transitions <- transitions[order(as.numeric(transitions[,1])),,drop=F]
  
  if (nrow(transitions) == 0){
    
    transitionsOut <- data.frame(ID=initPop[,'ID'], From= rep(NA,nrow(initPop)), 
                                 To=rep(NA,nrow(initPop)), transitionTime = rep(NA,nrow(initPop)), 
                                 transitionAge = rep(NA,nrow(initPop)), stringsAsFactors = FALSE)
    cat('Simulation has finished.\n')
    cat('Beware that along the simulation horizon the individual/s considered do/es not experience any transition/s.\n')
    cat('------------------\n')
    
  } else {
    
    cat('Simulation has finished.\n------------------\n')
    
    # ----------------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    # G. GENERATE OUTPUT 
    # ----------------------------------------------------------------------------------------------------------------------
    # ----------------------------------------------------------------------------------------------------------------------
    transitionsOut <- data.frame(ID=transitions[,'ID'], From=transitions[,'From'], To=transitions[,'To'],
                                 transitionTime = getInDateFormat(as.numeric(transitions[,'transitionTime'])), 
                                 transitionAge = round(as.numeric(transitions[,'transitionAge'])/365.25,2),
                                 stringsAsFactors = FALSE)
  }
  
  pop <- merge(initPop, transitionsOut, all=T, by='ID')
  pop <- pop[order(as.numeric(pop[,1])),] 
  
  if(length(fertTr)>0) {
   colnames(mothers) <- c("motherID", "ID")
   pop <- merge(pop, mothers, by="ID", all.x=TRUE) 
   pop <- pop[order(as.numeric(pop[,1])),] 
  }
  
  return(pop)
}

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# II. Execute microsimulation distributed (by executing as many single thread microsimulations in parallel as cores 
#     are available)
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
micSimParallel <- function(initPop=NULL, immigrPop=NULL, initPopList = c(), immigrPopList = c(),
                           transitionMatrix, absStates=NULL, varInitStates=c(), initStatesProb=c(), 
                           fixInitStates = c(), maxAge=99, simHorizon, fertTr=c(), monthSchoolEnrol=c(), cores=1, seeds=1254){
  
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
    
  sfInit(parallel=T,cpus=cores,slaveOutfile="output.txt")        
  sfExportAll(debug=FALSE)  
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


