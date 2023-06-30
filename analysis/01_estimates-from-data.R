# -------------------------------------------------------------------------
#' Estimate the frequency of resistant parasites from the data
#' 
# -------------------------------------------------------------------------


# Estimated fraction of resistant mutants ---------------------------------

# estimate the fraction of resistant mutants for each group in the data:
data.estimates <- data.frame(group=sort(unique(data.all$group)),estimate=c(rep(NA,length(unique(data.all$group)))),
                             CI_lower=c(rep(NA,length(unique(data.all$group)))),CI_upper=c(rep(NA,length(unique(data.all$group)))))

for(i in c(1:nrow(data.estimates))){
  # negative log-likelihood function for the current data subset (group):
  nllh_tmp <- function(log_f){nllh(log_f,data.all[data.all$group==data.estimates$group[i],])}
  
  # find the MLE by minimizing the negative log-likelihood function:
  fit_tmp <- nlm(nllh_tmp,log(initial_freq_guess),hessian=TRUE)
  
  # add estimate and Confidence Interval (CI) bounds to the data:
  data.estimates$estimate[i] <- exp(fit_tmp$estimate)
  data.estimates$CI_lower[i] <- exp(fit_tmp$estimate-qnorm((1+alpha_conf)/2)*sqrt(diag(solve(fit_tmp$hessian))))
  data.estimates$CI_upper[i] <- exp(fit_tmp$estimate+qnorm((1+alpha_conf)/2)*sqrt(diag(solve(fit_tmp$hessian))))
}


# add frequency as 1 resistant mutant per x parasites ---------------------

data.estimates <- cbind(data.estimates,freq_per_par=1/data.estimates$estimate,
                        freq_CI_lower=1/data.estimates$CI_upper,
                        freq_CI_upper=1/data.estimates$CI_lower)


# save output -------------------------------------------------------------

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(data.estimates,glue::glue("output/Data-estimates_{Sys.Date()}_{cur_time}.xlsx"))


# Estimate from all data combined -----------------------------------------

# use all data that are used in the analysis to estimate the overall fraction of resistant parasites:
nllh_tmp <- function(log_f){nllh(log_f,data.all[data.all$group%in%c("MIR_3D7","MIR_Dd2","NSG","VIS","field","ld_3D7_x52"),])}

# find the MLE by minimizing the negative log-likelihood function:
fit_tmp <- nlm(nllh_tmp,log(initial_freq_guess),hessian=TRUE)

# make data with the estimate and CI:
data_est_all <- exp(fit_tmp$estimate)
data_est_all_CI_lower <- exp(fit_tmp$estimate-qnorm((1+alpha_conf)/2)*sqrt(diag(solve(fit_tmp$hessian))))
data_est_all_CI_upper <- exp(fit_tmp$estimate+qnorm((1+alpha_conf)/2)*sqrt(diag(solve(fit_tmp$hessian))))


# Estimates assuming some wells are false negatives -----------------------

# parameters:
false_neg <- 0.05 # 5% of wells with at least one resistant parasite are false negatives

# estimates:
data_est_false_neg <- data.frame(group=sort(unique(data.all$group)),estimate=c(rep(NA,length(unique(data.all$group)))),
                             CI_lower=c(rep(NA,length(unique(data.all$group)))),CI_upper=c(rep(NA,length(unique(data.all$group)))))

for(i in c(1:nrow(data_est_false_neg))){
  # negative log-likelihood function for the current data subset (group):
  nllh_tmp <- function(log_f){nllh_false_neg(log_f,data.all[data.all$group==data.estimates$group[i],],false_neg)}
  
  # find the MLE by minimizing the negative log-likelihood function:
  fit_tmp <- nlm(nllh_tmp,log(initial_freq_guess),hessian=TRUE)
  
  # add estimate to the data:
  data_est_false_neg$estimate[i] <- exp(fit_tmp$estimate)
  data_est_false_neg$CI_lower[i] <- exp(fit_tmp$estimate-qnorm((1+alpha_conf)/2)*sqrt(diag(solve(fit_tmp$hessian))))
  data_est_false_neg$CI_upper[i] <- exp(fit_tmp$estimate+qnorm((1+alpha_conf)/2)*sqrt(diag(solve(fit_tmp$hessian))))
}

# add frequency as 1 resistant mutant per x parasites:
data_est_false_neg <- cbind(data_est_false_neg,freq_per_par=1/data_est_false_neg$estimate,
                        freq_CI_lower=1/data_est_false_neg$CI_upper,
                        freq_CI_upper=1/data_est_false_neg$CI_lower)

# save output:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(data_est_false_neg,glue::glue("output/Data-estimates-false-neg_{Sys.Date()}_{cur_time}.xlsx"))


# comparison of data estimates with and without false negatives -----------

data_estimates_all <- data.frame(data=data.estimates$group,estimate=data.estimates$estimate,CI_lower=data.estimates$CI_lower,
                                 CI_upper=data.estimates$CI_lower,false_neg_estimate=data_est_false_neg$estimate,
                                 false_neg_CI_lower=data_est_false_neg$CI_lower,false_neg_CI_upper=data_est_false_neg$CI_upper)
data_estimates_all <- dplyr::mutate(data_estimates_all,ratio=false_neg_estimate/estimate,
                                    ratio_CI_lower=false_neg_CI_lower/CI_upper,ratio_CI_upper=false_neg_CI_upper/CI_lower)

# save output:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(data_estimates_all,glue::glue("output/Data-estimates-ratio_{Sys.Date()}_{cur_time}.xlsx"))


# clean up ----------------------------------------------------------------

rm(i,nllh_tmp,fit_tmp,cur_time,false_neg)
