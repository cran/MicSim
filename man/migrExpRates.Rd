\name{migrExpRates}
\alias{migrExpRates}
\docType{data}
\title{
 Transition rates for migration example of MicSim package.
}
\description{
Transition rates for fertility (up to parity 4), migration between Spain / Netherlands / Sweden, mortality rates and emigration rates (leaving the Spain, the Netherlands, Sweden to some other country than these three). Rates have been estimated by the Unit ``Migration, Demography and Governance Unit'' of the European Commission.
}
\usage{data("migrExpRates")}
\format{
  A data frame with transition rates for 30 states and ages from 0 to 99.
  \describe{
    \item{\code{mort_f_ES}}{mortality rates for females in Spain}
    \item{\code{mort_f_NL}}{mortality rates for females in the Netherlands}
    \item{\code{mort_f_SE}}{mortality rates for females in Sweden}
    \item{\code{mort_m_ES}}{mortality rates for males in Spain}
    \item{\code{mort_m_NL}}{mortality rates for males in the Netherlands}
    \item{\code{mort_m_SE}}{mortality rates for males in Sweden}
    \item{\code{fert_ES_0_1}}{fertility rates for partity 1 for Spain}
    \item{\code{fert_ES_1_2}}{fertility rates for partity 2 for Spain}
    \item{\code{fert_ES_2_3}}{fertility rates for partity 3 for Spain}
    \item{\code{fert_ES_3_4}}{fertility rates for partity 4 for Spain}    
    \item{\code{fert_NL_0_1}}{fertility rates for partity 1 for the Netherlands}
    \item{\code{fert_NL_1_2}}{fertility rates for partity 2 for the Netherlands}
    \item{\code{fert_NL_2_3}}{fertility rates for partity 3 for the Netherlands}
    \item{\code{fert_NL_3_4}}{fertility rates for partity 4 for the Netherlands}   
    \item{\code{fert_SE_0_1}}{fertility rates for partity 1 for Sweden}
    \item{\code{fert_SE_1_2}}{fertility rates for partity 2 for Sweden}
    \item{\code{fert_SE_2_3}}{fertility rates for partity 3 for Sweden}
    \item{\code{fert_SE_3_4}}{fertility rates for partity 4 for Sweden}  
    \item{\code{rate_ES_NL}}{migration rates for Spain to the Netherlands}
    \item{\code{rate_ES_SE}}{migration rates for Spain to the Sweden}
    \item{\code{rate_NL_ES}}{migration rates for the Netherlands to Spain}
    \item{\code{rate_NL_SE}}{migration rates for the Netherlands to Sweden}  
    \item{\code{rate_SE_ES}}{migration rates for Sweden to Spain}  
    \item{\code{rate_SE_NL}}{migration rates for Sweden to the Netherlands}      
    \item{\code{emig_f_ES}}{emigration rates for females in Spain}
    \item{\code{emig_f_NL}}{emigration rates for females in the Netherlands}
    \item{\code{emig_f_SE}}{emigration rates for females in Sweden}
    \item{\code{emig_m_ES}}{emigration rates for females in Spain}  
    \item{\code{emig_m_NL}}{emigration rates for females in the Netherlands} 
    \item{\code{emig_m_SE}}{emigration rates for females in Sweden}   
  }
}
\source{
 European Commission
}
\examples{
data(migrExpRates)
}

