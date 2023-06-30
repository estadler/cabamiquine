# -------------------------------------------------------------------------
#' This file is used to run the project from start to finish
#' It can be run using the shortcut `cmd`/`ctrl` + `shift` + `S`
#' 
# -------------------------------------------------------------------------


# choose analysis to be done ----------------------------------------------

data_analysis <- TRUE # do the data analysis (estimates from the data, comparison of estimates)
model_estimates <- TRUE # plot deterministic and stochastic model simulations

# setup -------------------------------------------------------------------

source("setup.R")

# processing --------------------------------------------------------------

if(data_analysis){
  source("processing/01_make-data-table.R") # load data and clean data
}

# analysis ----------------------------------------------------------------

# run analysis and create outputs
if(data_analysis){
  source("analysis/01_estimates-from-data.R")
  source("analysis/02_data-estimates-comparison.R")
}
if(model_estimates){
  source("analysis/03_deterministic-model.R")
  source("analysis/04_stochastic-model.R")
  source("analysis/05_stochastic-sim-experiments.R")
}
if(data_analysis & model_estimates){
  source("analysis/06_visualisations.R")
}

