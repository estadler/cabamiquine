
# central estimate of the frequency of resistant mutants ------------------

# assume at least one resistant mutant is present at treatment if there is a recrudescence
# negative loglikelihood function for estimating the log-transformed fraction of resistant mutants
nllh <- function(log_f,data){ 
  # log_f: log-transformed fraction of resistant mutants
  # data: data used for the estimation, with column names inoculum, n_wells and n_pos_wells
  -sum(log(dbinom(data$n_pos_wells,data$n_wells,1-dbinom(0,data$inoculum,exp(log_f)))))
}


# negative loglikelihood function for the comparison of estimates ---------

# include a difference between two data sets as a variable of the function
nllh_d <- function(log_f,d,data){ 
  # log_f: log-transformed fraction of resistant mutants
  # d: addition/difference of the log-transformed fraction of resistant mutants
  # data: data used for the estimation, with column names inoculum, n_wells and n_pos_wells
  -sum(log(dbinom(data$n_pos_wells,data$n_wells,1-dbinom(0,data$inoculum,d+exp(log_f)))))
}


# estimate assuming x% of positive wells are false negatives --------------

# assume that (100-x)% of wells with at least 1 resistant parasite are positive wells but x% are false negatives
nllh_false_neg <- function(log_f,data,false_neg){ 
  # log_f: log-transformed fraction of resistant mutants
  # data: data used for the estimation, with column names inoculum, n_wells and n_pos_wells
  # false_neg: fraction of wells with at least one resistant parasite that show up as false negatives
  -sum(log(dbinom(data$n_pos_wells,data$n_wells,(1-false_neg)*(1-dbinom(0,data$inoculum,exp(log_f))))))
}
