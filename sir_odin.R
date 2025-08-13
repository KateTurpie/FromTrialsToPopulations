# -------------------------------------------------------------- #
# Script Name: ~/MResProjectTwo/sir_odin.R
# Script Purpose: Dynamic transmission model in odin
# Date Created: 2025-04-08
#
# Author: Kate Turpie
# Email: k.turpie24@imperial.ac.uk
# -------------------------------------------------------------- #

#### Odin Model ####

sirs_model <- odin2::odin({ 
  
  #### define parameters ####
  
  ## population parameters
  start_S <- parameter()
  start_Iu <- parameter()
  start_Ip <- parameter()
  start_R <- parameter()
  start_P <- parameter()
  start_I[] <- start_Ip[i] + start_Iu[i]
  n_age[] <- start_S[i] + start_I[i] + start_P[i] + start_R[i]
  initial(n_pop) <- sum(n_age)
  
  contacts <- parameter() #contact rates [matrix]
  contacts_ij[, ] <- contacts[i, j] * (I[j]/n_age[j]) #overall contacts by age group
  
  ## drug parameters
  chi <- parameter() #drug efficacy (initial)
  omega <- parameter() #duration of protection of drug
  rho <- parameter() #prophylactic coverage rate
  
  ## infection parameters 
  beta <- parameter() #infection rate (natural)
  beta_p <- beta*(1-chi) #infection rate (protected)
  
  lambda[] <- beta * sum(contacts_ij[i, ])
  lambda_p[] <- beta_p * sum(contacts_ij[i,])
  
  gamma <- parameter() #recovery rate
  delta <- parameter() #re-susceptibility rate
  
  ## model structure parameters
  age_number <- parameter()
  
  ## set array dimensions
  dim(S, I, Iu, Ip, R, P, n_age, age_groups) <- age_number #for age stratification
  dim(start_S, start_I, start_Iu, start_Ip, start_R, start_P,) <- age_number #for age stratification
  dim(p_SIu, p_SP,p_PS, p_PIp, p_IR, p_RS) <- age_number
  dim(lambda, lambda_p, tau, rho) <- age_number
  dim(n_SIu, n_SP, n_PS, n_PIp, n_IuR, n_IpR, n_RS) <- age_number
  dim(contacts, contacts_ij) <- c(age_number, age_number)
  
  #### Set Initial Conditions ####
  
  ## initial compartment values
  initial(S[]) <- start_S[i]
  initial(Iu[]) <- start_Iu[i]
  initial(Ip[]) <- start_Ip[i]
  initial(R[]) <- start_R[i]
  initial(P[]) <- start_P[i]
  initial(I[]) <- start_I[i]
  
  ## initial state transition rates 
  p_SP[] <- 1 - exp(-rho[i] * dt)
  p_PS[] <- 1 - exp(-omega * dt)
  
  p_SIu[] <- 1 - exp(-lambda[i] * dt)
  p_PIp[] <- 1 - exp(-lambda_p[i] * dt)
  
  p_IR[] <- 1 - exp(-gamma * dt)
  p_RS[] <- 1 - exp(-delta * dt)
  
  #### update number of individuals in each SIR state ####
  
  ## numbers flowing through
  n_SP[] <- Binomial(S[i], p_SP[i])
  n_PS[] <- Binomial(P[i], p_PS[i])
  
  n_SIu[] <- Binomial(S[i], p_SIu[i])
  n_PIp[] <- Binomial(P[i], p_PIp[i])
  
  n_IuR[] <- Binomial(Iu[i], p_IR[i])
  n_IpR[] <- Binomial(Ip[i], p_IR[i])
  
  n_RS[] <- Binomial(R[i], p_RS[i])
  
  ## update states
  update(S[]) <- S[i] - n_SIu[i] - n_SP[i] + n_RS[i] + n_PS[i]
  
  update(Iu[]) <- Iu[i] + n_SIu[i] - n_IuR[i]
  update(Ip[]) <- Ip[i] + n_PIp[i] - n_IpR[i]
  update(I[]) <- Iu[i] + Ip[i]
  
  update(R[]) <- R[i] + n_IuR[i] + n_IpR[i] - n_RS[i]
  update(P[]) <- P[i] + n_SP[i] - n_PS[i] - n_PIp[i]
  
  update(n_pop) <- sum(S) + sum(Iu) + sum(R) + sum(P) + sum(Ip)
  
  
  ## extras 
  
  dim(lambda_t, total_infections) <- age_number
  initial(lambda_t[]) <- 0
  update(lambda_t[]) <- lambda[i]
  initial(total_infections[]) <- start_I[i]
  update(total_infections[]) <- total_infections[i] + n_SIu[i] + n_PIp[i]
  
})
