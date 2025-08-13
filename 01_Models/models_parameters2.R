# -------------------------------------------------------------- #
# Script Name: ~/MResProjectTwo/Parameters/model_parameters2.R
# Script Purpose: Functions to set parameters for the models
# Date Created: 2025-04-11
#
# Author: Kate Turpie
# Email: k.turpie24@imperial.ac.uk
# -------------------------------------------------------------- #

#### Input Parameters for the Dynamic Transmission Model ####

## Infection and Prophylaxis Parameters ##

set_ip_params <- function(R0 = 3,
                          contacts = contact_params$contacts,
                          gamma = 7, 
                          delta = 275, 
                          chi=0, 
                          omega=0, 
                          rho=c(0.9,0,0,0,0,0,0,0)) {
  
  #Days to Rates
  if (delta != 0) {delta = 1/delta}
  if (gamma != 0) {gamma = 1/gamma}
  if (omega != 0) {omega = 1/omega}
  
  #Beta Calculation using NGM & Eigenvalues
  math_D <- 1/gamma
  
  contacts <- t(contacts) # c(0.9, 0.75, 0.5, 0.5, 0.5, 0.5, 0.5,0.5)
  NGM <- contacts*math_D
  
  eig <- eigen(NGM, only.values=TRUE)$values
  max_eig <- max(Re(eig))
  
  beta = R0/max_eig
  
  
  return(list(beta=beta, #infection rate
              gamma=gamma, #recovery rate
              delta=delta, #resusceptibility
              chi=chi, #efficacy
              omega=omega, # wane rate
              rho=rho)) #protection rate
}

## Population Parameters ##

set_init_compartments <- function(population=68350000, demog=contact_params$demog, set.seed_true=TRUE, seed=99, start_P_seed=TRUE) {
  
  #divide up susceptibles by age group
  start_S = demog*population
  
  #remove older adults
  start_R = round(start_S*c(0.1, 0.25, 0.5, 0.5, 0.5, 0.5, 0.5,0.5))
  start_S = start_S - start_R
  
  # seed infection
  if (set.seed_true == TRUE) {
    set.seed(seed)
  }
  
  proportion_Iu <- (sample(1:2,1))/100
  start_Iu <- round(start_S * proportion_Iu)
  start_S <- start_S - start_Iu
  
  
  # seed prophylaxis (if required)
  if (start_P_seed == FALSE) {
    start_P = rep(0, length(start_S))
  }
  else {
    start_P <- round(start_S * c(0.9, 0, 0, 0, 0, 0, 0,0))
    start_S <- start_S - start_P
  }
  
  # initialize other compartments to 0 across all age groups
  
  start_Ip <- rep(0, length(start_S))
  
  return(list(start_S = start_S, 
              start_Iu = start_Iu, 
              start_Ip = start_Ip, 
              start_P = start_P, 
              start_R = start_R))
}

set_contact_params <- function(country = "United Kingdom", 
                               age_limits = c(0,1,3,5,20,40,60,75), 
                               age_groups = c("0-1","1-2","3–4","5–19","20-39","40-59", "60-74","75+")){
  data(polymod)
  contacts <- contact_matrix(polymod, countries = country, age.limits = age_limits, symmetric = TRUE, return.demography=TRUE)
  age_number = length(age_groups)
  return(list(contacts = contacts$matrix, demog = contacts$demography$proportion, age_number = age_number, age_groups = age_groups))
}

## Combine all parameter lists for entry to Odin

set_all_params <- function(ip_params, init_compartments, contact_params){
  all_params <- c(ip_params, init_compartments, contact_params)
  return(all_params)
}

#### Input Parameters for Trial Model ####

set_trial_parameters <- function(
    beta = 0.3159, gamma = 1/7,
    efficacy = 1.0, wane_rate = 600, 
    n_exposed = 994, n_participants = 1490){
  
  wane_rate = 1/wane_rate
  
  params <- list(
    beta = beta, 
    gamma = gamma, 
    efficacy = efficacy, 
    wane_rate = wane_rate,
    n_exposed = n_exposed,
    n_participants = n_participants)
}

#### Model Settings Parameters ####

set_duration_and_runs <- function(trial_duration = 0:183, n_runs = 100) {
  
  duration_and_runs <- list(trial_duration = trial_duration, n_runs = n_runs)
}


