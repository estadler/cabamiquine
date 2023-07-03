# -------------------------------------------------------------------------
#' This file is used to set up everything that's needed across the project
#' It loads libraries, creates functions, sets themes and defaults
#' 
# -------------------------------------------------------------------------


# load packages -----------------------------------------------------------

library("readxl")
library("writexl")
library("reshape2")
library("ggplot2")
library("viridis")
library("matrixStats")

# load helper functions ----------------------------------------------------

source("helper_functions/estimates-from-data-helpers.R")
source("helper_functions/deterministic-model-functions.R")
source("helper_functions/stochastic-model-functions.R")
source("helper_functions/general-helper-functions.R")

# set defaults ------------------------------------------------------------

# for estimates from the data
alpha_conf <- 0.95 # confidence level, e.g. for confidence intervals
hRBC_vol <- 90e-12 # human RBC volume: 90e-12 mL

# defaults for the discrete model:
pmr <- 4 # PMR of parasites sensitive to cabamiquine
fitn_cost <- 1-7/100 # reduction factor of the PMR of resistant parasites, i.e. the PMR for resistant parasites is pmr*lambda, fitness cost of 7% reduction of PMR
mut_rate <- 1.05e-9 # mutation rate in base pair substitutions per base pair per generation
n_res_mut <- 11 # number of mutations to give resistance
initial_freq_guess <- 2e-8 # initial guess for the frequency of resistant mutants for the data estimates
init_res_par <- 0 # initial number of resistant parasites

# for plots
plot_width <- 8
plot_height <- 6
