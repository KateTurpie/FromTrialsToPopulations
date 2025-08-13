# -------------------------------------------------------------- #
# Script Name: ~/MResProjectTwo/Models/models_run.R
# Script Purpose: initalizes and runs the odin models
# Date Created: 2025-04-17
#
# Author: Kate Turpie
# Email: k.turpie24@imperial.ac.uk
# -------------------------------------------------------------- #

#### Run the Model Once ####

run_models <- function(model, params, duration_and_runs, stochastic=TRUE) {
  
  #### stochastic run ####
  
  if (stochastic == TRUE) {
    
    ## initialize system
    stochastic_model <- dust2::dust_system_create(model, 
                                           params, 
                                           n_particles=duration_and_runs$n_runs)
    dust2::dust_system_set_state_initial(stochastic_model)
    
    ## run
    stochastic_run <- dust2::dust_system_simulate(stochastic_model, 
                                         duration_and_runs$trial_duration)
    
    ## tidy
    rownames(stochastic_run) <- names(unlist(dust2::dust_unpack_index(stochastic_model)))
    stochastic_run <- tidyup(stochastic_run)
    
    return(stochastic_run)
  } 
  
  #### deterministic run ####
  
  if (stochastic == FALSE){
    
    ## initialize system
    deterministic_model <- dust2::dust_system_create(model, 
                                               params, 
                                               n_particles=duration_and_runs$n_runs, 
                                               deterministic = TRUE)
    dust2::dust_system_set_state_initial(deterministic_model)
    
    ## run
    deterministic_run <- dust2::dust_system_simulate(deterministic_model, 
                                             duration_and_runs$trial_duration)
    
    ## tidy
    rownames(deterministic_run) <- names(unlist(dust2::dust_unpack_index(deterministic_model)))
    deterministic_run <- tidyup(deterministic_run)
    
    return(deterministic_run)
  }
  
}

#### Run the Model Multiple Times Varying A Single Parameter ####

run_param_variation <- function(param_name, 
                                param_values, 
                                param_list=params, 
                                model, 
                                dar = duration_and_runs,
                                stochastic = FALSE){
  
  lapply(param_values, function(current){
    
    ## cycle through parameter values
    p <- param_list
    p[param_name] <- current
    
    ## run model with current parameter value
    run_models(model, p, dar,stochastic)
    
  })
}

#### Run the Model Multiple Times Varying Two Parameters ####

run_params_variations <- function(param_one, 
                                  one_values,
                                  param_two, 
                                  two_values,
                                  all_combos = TRUE,
                                  param_list = params, 
                                  model, 
                                  dar = duration_and_runs,
                                  stochastic = FALSE) {
  
  if (all_combos == TRUE) {
  param_grid <- expand_grid(
    one = one_values, 
    two = two_values
  )}
  
  else {param_grid <- data.frame(one = one_values, two = two_values)}
  
  result <- pmap(param_grid, function(one, two) {
    p <- param_list
    p[param_one] <- one
    p[param_two] <- 1/two 
    run_models(model, p, dar, stochastic)
  })
  
  param_grid$result <- result
  final_results <- param_grid %>% unnest(result)
  
  return(final_results)
}

#### Calculate All Combos of Initial Efficacy & Duration of Protection for a VE ####
## technically not a model but it runs ##

get_de_combos <- function(VE_targets = c(0.25, 0.5, 0.75, 0.9),
                          beta = 0.00034,
                          trial_len = 90,
                          omegas = seq(1/1000, 1/30, length.out = 500), 
                          end_t = 365) {
  
  valid_combinations <- list()  # empty list, will store valid combos per VE
  
  for (VE in VE_targets) {
    valid_omegas <- c()
    valid_chis <- c()
    actual_VE <- c()
    
    for (test_omega in omegas) {
      denom <- beta * (1 - exp(-test_omega * trial_len))
      log_arg <- 1 - VE * (1 - exp(beta * trial_len))
      
      if (log_arg > 0) {
        chi <- (test_omega * log(log_arg)) / denom
        
        if (chi > 0 && chi <= 1) {
          valid_omegas <- c(valid_omegas, test_omega)
          valid_chis <- c(valid_chis, chi)
          actual_VE <- c(actual_VE, VE)
        }}}
    
    all_valid_combos <- data.frame(
      VE = actual_VE,
      omega = valid_omegas,
      chi = valid_chis
    )
    
    valid_combinations[[as.character(VE)]] <- all_valid_combos
  }
  
  valid_combinations_all_VE <- do.call(rbind, valid_combinations)
  t <- seq(0, end_t, by = 1)
  
  efficacy_curves <- valid_combinations_all_VE %>%
    mutate(id = row_number()) %>% #assign ID number to all combinations
    crossing(time = t) %>%
    mutate(true_efficacy = chi*exp(-omega*t)) %>%
    mutate(VE = as.numeric(VE))
  
  return(list(valid_combos = valid_combinations_all_VE, efficacy_curves = efficacy_curves))
}
