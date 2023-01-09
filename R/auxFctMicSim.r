####################################################################################
####################################################################################
## AUXILIARY FUNCTIONS                                                            ## 
## - TO DEFINE MICROSIMULATION INPUT                                              ##
## - FOR COMPUTATION WITH DATES                                                   ##
## SZ, June 2022                                                                  ##
####################################################################################
####################################################################################

# Function computes the days that have pasted since 1970-01-01 up to the currDate (in format 'yyyymmdd')
getInDays <- function(currDate) {
  currDate <- as.numeric(currDate)
  cD_year  <- trunc(currDate/10000)
  cD_month <- trunc((round(currDate/10000,2)-cD_year)*100) 
  cD_days  <- trunc(currDate - (trunc(currDate/100)*100))
  cD_fullYearsInDays <- (cD_year - 1970)*365.25 
  daysInCurrYear <- cD_month*30.42 # approx. days in a month over the year
  cD_fracYearsInDays <- daysInCurrYear + (cD_days - 1)
  cD_daysSince01011970 <- cD_fullYearsInDays + cD_fracYearsInDays    
  return(cD_daysSince01011970)
}

# Function computes the correct age in days; arguments are the birth date and the current date in numeric format 'yyyymmdd'
getAgeInDays <- function(currDate, birthDate) {
  return(getInDays(currDate) - getInDays(birthDate))
}

# Get the number of days that have pasted from 1970-01-01 till 'yyyymm11'.
getInDays_my <- function(year, month){
  return((year - 1970)*365.25 + month*30.42 - 30.41)
}

# Function computes from days since 01-01-1970: year
getYear <- function(daysSince01011970){
  return(trunc(1970+daysSince01011970/365.25))
}

# Function computes from days since 01-01-1970: month (approx.)
getMonth <- function(daysSince01011970){
  y <- getYear(daysSince01011970) 
  fracInDays <- ((1970+daysSince01011970/365.25) - y)*365.25
  return(trunc(fracInDays/30.42)+1)
}

# Function computes from days since 01-01-1970: day (approx.)
getDay <- function(daysSince01011970){
  y <- getYear(daysSince01011970) 
  fracInDays <- ((1970+daysSince01011970/365.25) - y)*365.25
  month_b <- trunc(fracInDays/30.42)
  return((fracInDays- month_b*30.42)+1)
}

getInDateFormat <- function(daysSince01011970){
  y <- getYear(daysSince01011970) 
  m <- getMonth(daysSince01011970)
  d <- trunc(getDay(daysSince01011970))
  y <- ifelse(m %in% 13, y+1, y) # check for 13th month (occurs rarely due to rounding, simply replace by correct date)
  d <- ifelse(m %in% 13, 1, d)
  m <- ifelse(m %in% 13, 1, m)
  falseDates <- c("231", "431", "631", "931", "1131") # check for 31th in month
  conFD <- paste(m,d, sep="") %in% falseDates
  m <- ifelse(conFD, m+1, m)  
  d <- ifelse(conFD, 1, d)   
  cond0230 <- paste(m,d, sep="") %in% "230"  # check for 30th in Feb (irrespective of leap years)  
  m <- ifelse(cond0230, m+1, m)  
  d <- ifelse(cond0230, 1, d)     
  leapYear <- (y%%100!=0 & y%%4==0) | (y%%400==0) # check for 29th in Feb in years that are not leap years
  cond0229 <- paste(m,d, sep="") %in% "229"  
  m <- ifelse(!leapYear & cond0229, m+1,m)
  d <- ifelse(!leapYear & cond0229, 1,d)
  m <- ifelse(nchar(m)<2,as.character(paste(0,m, sep="")),m)
  d <- ifelse(nchar(d)<2,as.character(paste(0,d, sep="")),d)    
  return(apply(cbind.data.frame(cbind.data.frame(y,m),d), 1, paste, collapse ="")        )
}

# Construct matrix indicating transition pattern and naming the corresponding transition rate functions.
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
  # Function to identify whether a set of attributes (`substates') is part of a state space state
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