# -------------------------------------------------------------- #
# Script Name: ~/MResProjectTwo/Models/clinical_trial_odin.R
# Script Purpose: Odin model for the clinical trial simulation
# Date Created: 2025-04-17
#
# Author: Kate Turpie
# Email: k.turpie24@imperial.ac.uk
# -------------------------------------------------------------- #

#### Odin Model ####

trial_simulation <- odin2::odin({
  
  ## Define Parameters ##
  
  beta <- parameter() #infection rate
  gamma <- parameter() #recovery rate
  
  efficacy <- parameter() #initial efficacy of prophylactic
  wane_rate <- parameter() #waning rate, function of time
  
  n_participants <- parameter() #population size
  n_control <- n_participants - n_exposed
  n_exposed <- parameter()
 
  ## Set Array Dimensions ##
  
  dim(S, I, R, p_SI, n_SI, n_IR) <- 2
  
  ## Define State Transitions ##
  
  p_SI[1] <- 1 - exp(-beta*dt) #constant FoI, no drug protection [control]
  p_SI[2] <- 1 - exp(-beta * (1-(efficacy*exp(-wane_rate*time))) * dt) #constant FoI with drug [exposed]
  p_IR <- 1 - exp(-gamma * dt)
  
  ## Set Initial Values ##
  
  initial(S[1]) <- n_control
  initial(S[2]) <- n_exposed
  initial(I[1]) <- 0
  initial(I[2]) <- 0
  initial(R[1]) <- 0
  initial(R[2]) <- 0
  
  ## Update States ##
  
  n_SI[] <- Binomial(S[i], p_SI[i])
  n_IR[] <- Binomial(I[i], p_IR)
  
  update(S[]) <- S[i] - n_SI[i]
  update(I[]) <- I[i] + n_SI[i] - n_IR[i]
  update(R[]) <- R[i] + n_IR[i]
  
  ## Other Outputs ##
  
  # Total Cases Over Time #
  dim(total_cases) <- 2
  initial(total_cases[]) <- 0
  update(total_cases[]) <- total_cases[i] + n_SI[i]
  
  # Drug Efficacy Over Time #
  initial(efficacy_true) <- efficacy
  update(efficacy_true) <- (efficacy*exp(-wane_rate*time))
  
  initial(efficacy_measured) <- 0
  update(efficacy_measured) <-((total_cases[1]/n_control)-(total_cases[2]/n_exposed))/(total_cases[1]/n_control)
  
})

