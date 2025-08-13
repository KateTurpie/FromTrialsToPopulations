# -------------------------------------------------------------- #
# Script Name: ~/MResProjectTwo/Run/model_simulation.R
# Script Purpose: Run the models 
# Date Created: 2025-06-02
#
# Author: Kate Turpie
# Email: k.turpie24@imperial.ac.uk
# -------------------------------------------------------------- #

#### Set Up, Packages & Sourcing ####

source("../FromTrialsToPopulations/01_Models/models_set_up.R")

#### Clinical Trial Model ####

## Parameters ##

foi <- retrieve_foi(0,40)
params <- set_trial_parameters(beta=0.00034, wane_rate = 1/600, efficacy = 0.8477305)
duration_and_runs <- set_duration_and_runs(0:150,)

## Run Clinical Trial ##

trial_run <- run_models(trial_simulation, 
                       params,
                       duration_and_runs,
                       FALSE) # Toggle: True for stochastic, False for deterministic

## Save Run ##

save_run(trial_run)

#### SIRS DTM Model ####

## Set Parameters ##

contact_params <- set_contact_params()
ip_params <- set_ip_params(chi = 0.9997, omega = 1/73, rho = c(0.0,0,0,0,0,0,0,0,0))
init_compartments <- set_init_compartments()

params <- set_all_params(ip_params, contact_params, init_compartments)
duration_and_runs <- set_duration_and_runs(0:365,2)

## Run DTM Trial ##

sir_run <- run_models(sirs_model, 
                      params, 
                      duration_and_runs, 
                      FALSE)  # Toggle: True for stochastic, False for deterministic

## Save Run ##

save_run(sir_run)

