#' Initial Population Data for Examples on Simulating Migration Flows
#'
#' Data used to construct initial population for Migration Flow Vignette
#'
#'
#' @format ## `immigrPopMigrExp`
#' A data frame with 3758 rows and 4 columns:
#' \describe{
#'   \item{ID}{ID of Individual}
#'   \item{immigrDate}{Immigration Date}
#'   \item{birthDate}{Birth Date}
#'   \item{immigrInitState}{initial State of Individual as per stateSpace}
#' }

"immigrPopMigrExp"

#' Immigrant Population Data for Examples on Simulating Migration Flows
#'
#' Data used to construct initial immigrant population for Migration Flow Vignette
#'
#'
#' @format ## `initPopMigrExp`
#' A data frame with 72965 rows and 3 columns:
#' \describe{
#'  \item{ID}{Id of Individual}
#'  \item{birthDate}{Birth Date}
#'  \item{immigrInitState}{initial State of Individual as per stateSpace}
#' }

"initPopMigrExp"

#' Rate Data for Examples on Simulating Migration Flows
#'
#' Data used to construct transition rates for Migration Flow Vignette
#'
#'
#' @format ## `migrExpRates`
#' A data frame with 100 rows and 30 columns:
#' \describe{
#'   \item{mort_f_ES}{Mortality Rate for females from Spain}
#'   \item{mort_f_NL}{Mortality Rate for females from the Netherlands}
#'   \item{mort_f_SE}{Mortality Rate for females from Swedem}
#'   ...
#' }

"migrExpRates"

#' Rate Data for Simulation of MDD example
#'
#' Data used to construct transition rates from @lepe2024
#'
#'
#' @format ## `rates_mdd`
#' A data frame with 49 rows and 17 columns:
#' \describe{
#'   \item{age}{Age of individual}
#'   \item{incidence_lo_m}{Incidence Rates low education men}
#'   \item{cf_incidence_lo_f}{Counterfactual Incidence Rate low female}
#'   ...
#' }
#' @source https://doi.org/10.1093/eurpub/ckae066
"rates_mdd"