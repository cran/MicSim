#' getInDays
#'
#' Get from a date given in the numeric format yyyymmdd the number of days elapsed since 1970-01-01
#'
#' @usage getInDays(currDate)
#' @description Function computes the days that have pasted since 1970-01-01 up to the currDate (in the numeric format yyyymmdd)
#'
#' @param currDate current date - Date given as string of the numeric format yyyymmdd.
#'
#' @rdname getInDays
#' @returns cD_daysSince01011970 - Number of days elapsed since 1970-01-01.
#' @examples getInDays(20200826)
#' @export
#'
#'
getInDays <- function(currDate) {
  currDate <- as.numeric(currDate)
  cD_year  <- currDate %/% 10000
  cD_month <- (currDate %/% 100) %% 100
  cD_days  <- currDate %% 100
  cD_fullYearsInDays <- (cD_year - 1970)*365.25
  daysInCurrYear <- (cD_month - 1)*30.42 # approx. days in a month over the year
  cD_fracYearsInDays <- daysInCurrYear + (cD_days - 1)
  cD_daysSince01011970 <- cD_fullYearsInDays + cD_fracYearsInDays
  return(cD_daysSince01011970)
}
#'
#' getAgeInDays
#'
#' Get from a given date the age in days
#'
#' @usage getAgeInDays(currDate, birthDate)
#' @description Function computes for a given date the correct age in days.
#'
#' @param currDate current date - Reference date given as string of the format "yyyymmdd"
#' @param birthDate Birth date given as string of the format "yyyymmdd".
#'
#' @rdname getAgeInDays
#' @returns Age in Days - Correct age at the specific date \code{currDate} in days
#' @examples getAgeInDays("20200826", "19800605")
#' @export
#'
#'
getAgeInDays <- function(currDate, birthDate) {
  return(getInDays(currDate) - getInDays(birthDate))
}
#'
#' getInDays_my
#'
#' Get the number of days that have passed from 1970-01-01 till 'yyyymm11'.
#'
#' @usage getInDays_my(year, month)
#' @description Function computes the number of days that have pasted from 1970-01-01 until 'yyyymm11'.
#'
#' @param year Year for which days elapsed should be computed, i.e., the yyyy in 'yyyymm11'
#' @param month Month for which days elapsed should be computed, i.e., the mm in 'yyyymm11'
#'
#' @rdname getInDays_my
#' @returns Number of days that have pasted from 1970-01-01 until 'yyyymm11'
#' @examples getInDays_my(2020, 12)
#' @export
#'
#'
getInDays_my <- function(year, month){
  return((year - 1970)*365.25 + month*30.42 - 30.41)
}
#'
#' getYear
#'
#' Get the calendar year from days elapsed since 01-01-1970
#'
#' @usage getYear(daysSince01011970)
#' @description Function computes from days elapsed since 01-01-1970 the related calendar year.
#'
#' @param daysSince01011970 Days elapsed since 1970-01-01
#'
#' @rdname getYear
#' @returns Calendar year from days elapsed since 01-01-1970
#' @examples getYear(2561)
#' @export
#'
#'
getYear <- function(daysSince01011970){
  return(trunc(1970+daysSince01011970/365.25))
}
#'
#' getMonth
#'
#' Get the month in a year from days elapsed since 01-01-1970
#'
#' @usage getMonth(daysSince01011970)
#' @description Function computes from days elapsed since 01-01-1970 the related month a year.
#'
#' @param daysSince01011970 Days elapsed since 1970-01-01
#'
#' @returns Month in a year computed from days elapsed since 01-01-1970
#' @examples getMonth(2561)
#' @rdname getMonth
#' @export
#'
#'
getMonth <- function(daysSince01011970){
  y <- getYear(daysSince01011970)
  fracInDays <- ((1970+daysSince01011970/365.25) - y)*365.25
  return(trunc(fracInDays/30.42)+1)
}
#'
#' getDay
#'
#' Get the day in a month (in a year) from days elapsed since 01-01-1970
#'
#' @usage getDay(daysSince01011970)
#' @description Function computes from days elapsed since 01-01-1970 the day in a month (in a year).
#'
#' @param daysSince01011970 Days elapsed since 1970-01-01
#'
#' @returns Day in a month (in a year) computed from days elapsed since 01-01-1970
#' @examples getDay(2561)
#' @rdname getDay
#' @export
#'
#'
getDay <- function(daysSince01011970){
  y <- getYear(daysSince01011970)
  fracInDays <- ((1970+daysSince01011970/365.25) - y)*365.25
  month_b <- trunc(fracInDays/30.42)
  return((fracInDays- month_b*30.42)+1)
}
#'
#' getInDateFormat
#'
#' Get date in the format 'yyyyddmm' from days elapsed since 01-01-1970
#'
#' @usage getInDateFormat(daysSince01011970)
#' @description Function generates from days elapsed since 01-01-1970 the date in the string format 'yyyyddmm'.
#'
#' @param daysSince01011970 Days elapsed since 1970-01-01
#'
#' @returns Date in string format 'yyyyddmm' from days elapsed since 01-01-1970
#' @examples getInDateFormat(2561)
#' @rdname getInDateFormat
#' @export
#'
#'
getInDateFormat <- function(daysSince01011970){
  y <- getYear(daysSince01011970)
  m <- getMonth(daysSince01011970)
  d <- trunc(getDay(daysSince01011970))
  return(sprintf("%4d%02d%02d", y, m, d))
}


