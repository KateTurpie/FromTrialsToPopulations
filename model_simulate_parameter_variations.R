# -------------------------------------------------------------- #
# Script Name: ~/MResProjectTwo/02_Run/model_simulate_parameter_variations.R
# Script Purpose: run the models with different parameter variations
# Date Created: 2025-06-23
#
# Author: Kate Turpie
# Email: k.turpie24@imperial.ac.uk
# -------------------------------------------------------------- #

#### Set Up, Packages & Sourcing ####

source("../MResProjectTwo/01_Models/models_set_up.R")

#### Clinical Trial Model ####

## Parameters ##
params <- set_trial_parameters(beta=0.00034)
duration_and_runs <- set_duration_and_runs(0:365,1)

## Parameter Variations ##
wane_rate_values <- seq(30,1080,30)
efficacy_values <- seq(0.2,1, 0.05)

## Run Model, varying one parameter ##
efficacy_1_0.5 <- run_param_variation("efficacy", efficacy_values, model = trial_simulation, stochastic = TRUE)
wane_ranges <- run_param_variation("wane_rate", wane_rate_values, model = trial_simulation)

## Run Model, varying two parameters ##
trial_runs <- run_params_variations("efficacy", efficacy_values, "wane_rate", wane_rate_values, model = trial_simulation, stochastic = FALSE)


#### DTM Model ####

## Parameters ##
contact_params <- set_contact_params()
ip_params <- set_ip_params(chi = 0.9997, omega = 1/73, rho = c(0.0,0,0,0,0,0,0,0))
init_compartments <- set_init_compartments()

params <- set_all_params(ip_params, contact_params, init_compartments)
duration_and_runs <- set_duration_and_runs(0:365,1)

# Parameters from the Equation #
VE_targets <- c(0.9, 0.75, 0.5, 0.25)
omegas <- 1/(1000:30)
trial_lengths <- c(30, 90, 150, 180, 270, 360)

all_sir_runs <- list()

for (trial_len in trial_lengths) {
  efficacies <- get_de_combos(trial_len = trial_len, omegas = omegas)
  
  sir_runs_for_trial <- list()
  
  for (ve in VE_targets) {
    VE_eff <- efficacies$valid_combos$VE == ve
    efficacy_values <- efficacies$valid_combos$chi[VE_eff]
    wane_rate_values <- 1 / efficacies$valid_combos$omega[VE_eff]
    
    sir_runs <- run_params_variations(
      "chi", efficacy_values,
      "omega", wane_rate_values,
      all_combos = FALSE,
      model = sirs_model,
      stochastic = FALSE
    )
    
    # Name keys like ve90, ve75, etc.
    key_name <- paste0("ve", ve * 100)
    sir_runs_for_trial[[key_name]] <- sir_runs
  }
  
  all_sir_runs[[paste0("sir_runs_", trial_len)]] <- sir_runs_for_trial
}

# Incorrect (but purposeful) Parameters # 
efficacy_values <- c(0.9,0.75,0.5,0.25,0.9, 0.75, 0.5,0.25, 0.9, 0.75, 0.5,0.25, 0.9, 0.75, 0.5,0.25)
wane_rate_values <- c(1000000000000000000000,1000000000000000000000,1000000000000000000000,1000000000000000000000,1000, 1000, 1000, 1000, 500, 500, 500, 500, 270,270,270,270)

sir_runs <- run_params_variations("chi", efficacy_values, "omega", wane_rate_values, all_combos=FALSE, model = sirs_model, stochastic = FALSE)


## Baseline ##
ip_params <- set_ip_params(rho = c(0.0,0,0,0,0,0,0,0))
init_compartments <- set_init_compartments(start_P_seed = FALSE)
params <- set_all_params(ip_params, contact_params, init_compartments)
sir_baseline <- run_models(sirs_model, 
                           params, 
                           duration_and_runs, 
                           FALSE)  # Toggle: True for stochastic, False for deterministic
