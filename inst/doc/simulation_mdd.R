## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----libraries, warning=FALSE, message=FALSE----------------------------------
# loading libraries
library(MicSim)
library(knitr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Clearing the environment
rm(list = ls())
invisible(gc())

## ----sim frame----------------------------------------------------------------
startDate <- 20180101
endDate   <- 20651231
simHorizon <- c(startDate = startDate, endDate = endDate)

N <- 1000

minage <- 18
maxage <- 65

## ----sim frame state space----------------------------------------------------
health    <- c("NoMDD", "MDD")
education <- c("Low", "High")
sex   <- c("Male", "Female")
stateSpace <- expand.grid(health = health, education = education, sex = sex)

absStates <- "dead"

## ----init states--------------------------------------------------------------
set.seed(1234)

birthDates <- runif(N, min = getInDays(20000101), max = getInDays(20001231))

prev_lo_m <- 0.0122
prev_lo_f <- 0.0191
prev_hi_m <- 0.0014
prev_hi_f <- 0.0022

m_lo_mdd <- sample(
  x = c("MDD/Low/Male", "NoMDD/Low/Male"),
  prob = c(prev_lo_m, 1 - prev_lo_m),
  size = N / 4, replace = TRUE
)

f_lo_mdd <- sample(
  x = c("MDD/Low/Female", "NoMDD/Low/Female"),
  prob = c(prev_lo_f, 1 - prev_lo_f),
  size = N / 4, replace = TRUE
)

m_hi_mdd <- sample(
  x = c("MDD/High/Male", "NoMDD/High/Male"),
  prob = c(prev_hi_m, 1 - prev_hi_m),
  size = N / 4, replace = TRUE
)

f_hi_mdd <- sample(
  x = c("MDD/High/Female", "NoMDD/High/Female"),
  prob = c(prev_hi_f, 1 - prev_hi_f),
  size = N / 4, replace = TRUE
)

initStates <- c(m_lo_mdd, f_lo_mdd, m_hi_mdd, f_hi_mdd)

## ----init pop-----------------------------------------------------------------

initPop <- data.frame(ID = 1:N, birthDate = birthDates, initState = initStates)

initPop$birthDate <- getInDateFormat(initPop$birthDate)

initPop %>% 
  head() %>% 
  kable("pipe", caption = "**Example of what initPop should look like**")


## ----transition rates data----------------------------------------------------
mortRates <- function(age, calTime) {
  return(0)
}

rates <- rates_mdd

rates %>%
  head(n = c(6,8)) %>%
  kable("html", caption = "**Subset of transition rates per age**")

## ----transition rates---------------------------------------------------------

create_rate_function <- function(vec) {
  # Ensures the correct column of transition rates called on by the function
  force(vec)
  
  function(age, calTime) {
    rate <- vec[as.integer(age) - 16]
    return(rate)
  }
}

# Loop over relevant columns to create rate functions.
for (col in names(rates)[-1]) {
  # Extracting relevant transition rates
  vec <- as.numeric(rates[[col]])

  # Creating a function for each transition rate
  f <- create_rate_function(vec)

  # Assign the function to the global environment
  assign(col, f, envir = .GlobalEnv)
}

## ----transition definition----------------------------------------------------
absTransitions <- c("dead", "mortRates")

allTransitions <- cbind(
  c("NoMDD/Low/Male->MDD/Low/Male", "NoMDD/Low/Female->MDD/Low/Female",
    "NoMDD/High/Male->MDD/High/Male", "NoMDD/High/Female->MDD/High/Female",
    "MDD/Low/Male->NoMDD/Low/Male", "MDD/Low/Female->NoMDD/Low/Female",
    "MDD/High/Male->NoMDD/High/Male", "MDD/High/Female->NoMDD/High/Female"),
  c("incidence_lo_m", "incidence_lo_f",
    "incidence_hi_m", "incidence_hi_f",
    "remission_lo_m", "remission_lo_f",
    "remission_hi_m", "remission_hi_f")
)

transitionMatrix <- buildTransitionMatrix(
  allTransitions = allTransitions,
  stateSpace = stateSpace,
  absTransitions = absTransitions
)

## ----transition matrix--------------------------------------------------------
cf_allTransitions <- allTransitions
cf_allTransitions[, 2] <- paste0("cf_", cf_allTransitions[, 2])

cf_transitionMatrix <- buildTransitionMatrix(
  allTransitions = cf_allTransitions,
  stateSpace = stateSpace,
  absTransitions = absTransitions
)

## ----run sims-----------------------------------------------------------------
set.seed(1234)

simpop <- micSim(initPop=initPop, transitionMatrix=transitionMatrix,
                 absStates=absStates, maxAge=maxage, simHorizon=simHorizon)

cf_simpop <- micSim(initPop=initPop, transitionMatrix=cf_transitionMatrix,
                    absStates=absStates, maxAge=maxage, simHorizon=simHorizon)

## ----processing output--------------------------------------------------------
simpop_long <- simpop %>%
  convertToLongFormat() %>%
  mutate(
    agestart = getAgeInDays(Tstart, birthDate) / 365.25,
    agestart = case_when(
      agestart < minage ~ minage,
      agestart > maxage ~ maxage,
      TRUE ~ agestart
    ),
    agestop = getAgeInDays(Tstop, birthDate) / 365.25,
    agestop = case_when(
      agestop < minage ~ minage,
      agestop > maxage ~ maxage,
      TRUE ~ agestop
    ),
    sim = "Observed"
  )

cf_simpop_long <- cf_simpop %>%
  convertToLongFormat() %>%
  mutate(
    agestart = getAgeInDays(Tstart, birthDate) / 365.25,
    agestart = case_when(
      agestart < minage ~ minage,
      agestart > maxage ~ maxage,
      TRUE ~ agestart
    ),
    agestop = getAgeInDays(Tstop, birthDate) / 365.25,
    agestop = case_when(
      agestop < minage ~ minage,
      agestop > maxage ~ maxage,
      TRUE ~ agestop
    ),
    sim = "Counterfactual"
  )

mydata <- simpop_long %>%
  bind_rows(cf_simpop_long) %>%
  mutate(
    sim = factor(sim, levels = c("Observed", "Counterfactual")),
    table_order = case_when(
      sim == "Observed" & education == "High" ~ 0,
      sim == "Observed" & education == "Low" ~ 1,
      sim == "Counterfactual" & education == "Low" ~ 2
    )
  )

## ----age-specific prev, fig.width=8,fig.height=8------------------------------
prev <- list()
for (i in minage:maxage){
  prev[[i]] <- simpop_long %>%
    filter((health == "MDD" & agestart <= i & agestop >= i)) %>%
    group_by(sex, education) %>%
    summarise(prev = n()/(N*.25) * 100, .groups = "keep") %>%
    mutate(age = i) %>%
    ungroup()
}

prev <- bind_rows(prev)

prev %>%
  ggplot() +
  geom_line(
    aes(age, prev, color = sex, linetype = education), 
    linewidth = 1
  ) +
  scale_y_continuous(breaks = c(2, 4, 6)) +
  labs(
    y = "Prevalence of MDD %",
    x = "Age (years)",
    colour = "Sex",
    linetype = "Education"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.margin = margin(0, 0, 0, 0)
  )

## ----life course prev---------------------------------------------------------
lc_prev <- mydata %>%
  filter(Episode <= 2, health == "MDD", !is.na(table_order)) %>%
  group_by(sim, education, sex) %>%
  summarise(proportion_MDD = round(n() / (N * .25) * 100, 1), .groups = "drop") %>%
  pivot_wider(
    names_from = c(sex), names_glue = "{sex} (n%)",
    values_from = proportion_MDD
  ) %>%
  rename(Scenario = sim, Education = education)


kable(lc_prev)

