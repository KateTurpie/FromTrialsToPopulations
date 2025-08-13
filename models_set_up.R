# -------------------------------------------------------------- #
# Script Name: ~/MResProjectTwo/Run/models_set_up.R
# Script Purpose: packages & source files for models
# Date Created: 2025-04-17
#
# Author: Kate Turpie
# Email: k.turpie24@imperial.ac.uk
# -------------------------------------------------------------- #

#### Load the Packages ####

library(odin2)
library(dust2)
library(monty)
library(tidyverse)
library(socialmixr)
library(ggthemes)
library(MoMAColors)
library(scales)

#### Source all Relevant Files & Functions ####

## Files required to run the models ##
source("../MResProjectTwo/01_Models/clinical_trial_odin.R")
source("../MResProjectTwo/01_Models/sir_odin.R")
source("../MResProjectTwo/01_Models/models_parameters2.R")
source("../MResProjectTwo/01_Models/models_run.R")

## Data tidying and such ##
source("../MResProjectTwo/04_Data/tidying.R")
source("../MResProjectTwo/04_Data/save_files.R")