# -------------------------------------------------------------------------
#' Model estimates for the frequency of resistant parasites
#' stochastic, synchronous generations model
#' Simulation of the experiments and estimation of the data in the same way
#' as in the original data to compare the data with the stochastic model
#' 
#' Includes sensitivity analysis for the stochastic simulations
#' 
# -------------------------------------------------------------------------

# specify simulations to do:
sim_do <- c("all") # all simulations: c("VIS","NSG","MIR_Dd2","MIR_3D7","LD_field","LD_3D7") or c("all")
sens_analysis_do <- TRUE # do a sensitivity analysis for all the experiments specified in sim_do


# parameters for the model ------------------------------------------------

# default paramters initialized in setup.R, and:
n_sim <- 1e5 # number of simulations of experiments
eq <- det_mod_equil_par(n_res_mut,mut_rate,fitn_cost,pmr) # equilibrium distribution of resistant parasites with standard parameters

# load the data-specific parameters: 
data_settings <- read_excel("raw-data/Data-information.xlsx",col_names = TRUE)


# parameters for the sensitivity analysis ---------------------------------

# list of parameters to use for the sensitivity analysis:
sens_analysis_list <- list(mut_rate=c(1e-10,3.2e-10,1.05e-9,5.6e-9,3.0e-8),
                           n_res_mut=c(5,11,20,30,40),
                           fitn_cost=1-c(3,5,7,9,11)/100,
                           pmr=c(1.2,2,4,5,6))
# dataframe with the standard parameters:
standard_param <- data.frame(mut_rate,n_res_mut,fitn_cost,pmr)


# simulation of VIS data --------------------------------------------------

if("VIS"%in%sim_do | any(sim_do=="all")){
  # parameter values for the VIS experiment:
  VIS_inoculum <- as.numeric(data_settings$inoculum[data_settings$data=="VIS"])
  VIS_PMR <- as.numeric(data_settings$PMR[data_settings$data=="VIS"])
  VIS_n_gen <- as.numeric(data_settings$generations[data_settings$data=="VIS"])
  VIS_n <- data_settings$n_cultures[data_settings$data=="VIS"]
  
  # make simulation:
  f_VIS <- stoch_sim_exps(eq,inoc=VIS_inoculum,n_inoc=VIS_n,p=VIS_PMR,n_gen=VIS_n_gen,n_sim)
  
  # make summary data from the simulation:
  f_VIS_summary <- data.frame(median=median(f_VIS$frac_res),CI_lower=quantile(f_VIS$frac_res,probs=(1-alpha_conf)/2),
                              CI_upper=quantile(f_VIS$frac_res,probs=(1+alpha_conf)/2),mean=mean(f_VIS$frac_res),
                              geom_mean=exp(mean(log(f_VIS$frac_res))))
  
  # save the simulated data: 
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  save(f_VIS,f_VIS_summary,file=glue::glue("output/Simulations-VIS-experiment_{Sys.Date()}_{cur_time}.RData"))
  
  # sensitivity analysis:
  if(sens_analysis_do){
    # sensitivity analysis
    tmp <- stoch_sim_exps_sens_analysis(exp_inoc=VIS_inoculum,exp_n_inoc=VIS_n,exp_pmr=VIS_PMR,exp_n_gen=VIS_n_gen,n_sim,
                                        param_standard=standard_param,sens_list=sens_analysis_list)
    
    # save sensitiviy analysis results:
    VIS_sens_analysis_all <- tmp$all
    VIS_sens_analysis_summary <- tmp$summary
    VIS_sens_analysis_pos <- tmp$frac_pos
    
    cur_time <- format(Sys.time(),"h%H_m%M_s%S")
    save(VIS_sens_analysis_all,VIS_sens_analysis_summary,VIS_sens_analysis_pos,
         file=glue::glue("output/VIS-exp-sensitivity-analysis_{Sys.Date()}_{cur_time}.RData"))
  }
  
  # cleanup:
  rm(VIS_inoculum,VIS_PMR,VIS_n_gen,VIS_n,f_VIS,tmp,VIS_sens_analysis_all,cur_time)
}


# simulation of NSG data --------------------------------------------------

if("NSG"%in%sim_do | any(sim_do=="all")){
  # parameter values for the NSG experiment:
  NSG_inoculum <- as.numeric(data_settings$inoculum[data_settings$data=="NSG"])
  NSG_PMR <- as.numeric(data_settings$PMR[data_settings$data=="NSG"])
  NSG_n_gen <- as.numeric(data_settings$generations[data_settings$data=="NSG"])
  NSG_n_mice <- data_settings$n_cultures[data_settings$data=="NSG"] # number of mice in the data
  
  # make simulation:
  f_NSG <- stoch_sim_exps(eq,inoc=NSG_inoculum,n_inoc=NSG_n_mice,p=NSG_PMR,n_gen=NSG_n_gen,n_sim)
  
  # make summary data from the simulation:
  f_NSG_summary <- data.frame(median=median(f_NSG$frac_res),CI_lower=quantile(f_NSG$frac_res,probs=(1-alpha_conf)/2),
                              CI_upper=quantile(f_NSG$frac_res,probs=(1+alpha_conf)/2),mean=mean(f_NSG$frac_res),
                              geom_mean=exp(mean(log(f_NSG$frac_res))))
  
  # save the simulated data: 
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  save(f_NSG,f_NSG_summary,file=glue::glue("output/Simulations-NSG-experiment_{Sys.Date()}_{cur_time}.RData"))
  
  # sensitivity analysis:
  if(sens_analysis_do){
    # sensitivity analysis
    tmp <- stoch_sim_exps_sens_analysis(exp_inoc=NSG_inoculum,exp_n_inoc=NSG_n_mice,exp_pmr=NSG_PMR,exp_n_gen=NSG_n_gen,n_sim,
                                        param_standard=standard_param,sens_list=sens_analysis_list)
    
    # save sensitiviy analysis results:
    NSG_sens_analysis_all <- tmp$all
    NSG_sens_analysis_summary <- tmp$summary
    NSG_sens_analysis_pos <- tmp$frac_pos
    
    cur_time <- format(Sys.time(),"h%H_m%M_s%S")
    save(NSG_sens_analysis_all,NSG_sens_analysis_summary,NSG_sens_analysis_pos,
         file=glue::glue("output/NSG-exp-sensitivity-analysis_{Sys.Date()}_{cur_time}.RData"))
  }
  
  # cleanup:
  rm(NSG_inoculum,NSG_PMR,NSG_n_gen,NSG_n_mice,tmp,NSG_sens_analysis_all,f_NSG,cur_time)
}


# simulation of MIR Dd2 data ----------------------------------------------

if("MIR_Dd2"%in%sim_do | any(sim_do=="all")){#DONE
  # parameter values for the MIR_Dd2 experiment:
  MIR_Dd2_inocula <- data_settings$inoculum[data_settings$data=="MIR in vitro Dd2"]
  MIR_Dd2_n_inocula <- data_settings$n_cultures[data_settings$data=="MIR in vitro Dd2"]
  MIR_Dd2_PMR <- as.numeric(data_settings$PMR[data_settings$data=="MIR in vitro Dd2"])
  MIR_Dd2_n_gen <- as.numeric(data_settings$generations[data_settings$data=="MIR in vitro Dd2"])
  
  # make simulation:
  f_MIR_Dd2 <- stoch_sim_exps(eq,MIR_Dd2_inocula,MIR_Dd2_n_inocula,p=MIR_Dd2_PMR,n_gen=MIR_Dd2_n_gen,n_sim)
  
  # make summary data from the simulation:
  f_MIR_Dd2_summary <- data.frame(median=median(f_MIR_Dd2$frac_res),CI_lower=quantile(f_MIR_Dd2$frac_res,probs=(1-alpha_conf)/2),
                                CI_upper=quantile(f_MIR_Dd2$frac_res,probs=(1+alpha_conf)/2),mean=mean(f_MIR_Dd2$frac_res),
                                geom_mean=exp(mean(log(f_MIR_Dd2$frac_res))))
  
  # save the simulated data: 
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  save(f_MIR_Dd2,f_MIR_Dd2_summary,file=glue::glue("output/Simulations-MIR-Dd2-experiment_{Sys.Date()}_{cur_time}.RData"))
  
  # sensitivity analysis:
  if(sens_analysis_do){
    # sensitivity analysis
    tmp <- stoch_sim_exps_sens_analysis(exp_inoc=MIR_Dd2_inocula,exp_n_inoc=MIR_Dd2_n_inocula,exp_pmr=MIR_Dd2_PMR,exp_n_gen=MIR_Dd2_n_gen,n_sim,
                                        param_standard=standard_param,sens_list=sens_analysis_list)
    
    # save sensitiviy analysis results:
    MIR_Dd2_sens_analysis_all <- tmp$all
    MIR_Dd2_sens_analysis_summary <- tmp$summary
    MIR_Dd2_sens_analysis_pos <- tmp$frac_pos
    
    cur_time <- format(Sys.time(),"h%H_m%M_s%S")
    save(MIR_Dd2_sens_analysis_all,MIR_Dd2_sens_analysis_summary,MIR_Dd2_sens_analysis_pos,
         file=glue::glue("output/MIR-Dd2-exp-sensitivity-analysis_{Sys.Date()}_{cur_time}.RData"))
  }
  
  # cleanup:
  rm(MIR_Dd2_inocula,MIR_Dd2_n_inocula,MIR_Dd2_PMR,MIR_Dd2_n_gen,f_MIR_Dd2,tmp,MIR_Dd2_sens_analysis_all,cur_time)
}


# simulation of MIR 3D7 data ----------------------------------------------

if("MIR_3D7"%in%sim_do | any(sim_do=="all")){
  # parameter values for the MIR_3D7 experiment:
  MIR_3D7_inocula <- data_settings$inoculum[data_settings$data=="MIR in vitro 3D7"]
  MIR_3D7_n_inocula <- data_settings$n_cultures[data_settings$data=="MIR in vitro 3D7"]
  MIR_3D7_PMR <- as.numeric(data_settings$PMR[data_settings$data=="MIR in vitro 3D7"])
  MIR_3D7_n_gen <- as.numeric(data_settings$generations[data_settings$data=="MIR in vitro 3D7"])
  
  # make simulation:
  f_MIR_3D7 <- stoch_sim_exps(eq,MIR_3D7_inocula,MIR_3D7_n_inocula,p=MIR_3D7_PMR,n_gen=MIR_3D7_n_gen,n_sim)
  
  # make summary data from the simulation:
  f_MIR_3D7_summary <- data.frame(median=median(f_MIR_3D7$frac_res),CI_lower=quantile(f_MIR_3D7$frac_res,probs=(1-alpha_conf)/2),
                                  CI_upper=quantile(f_MIR_3D7$frac_res,probs=(1+alpha_conf)/2),mean=mean(f_MIR_3D7$frac_res),
                                  geom_mean=exp(mean(log(f_MIR_3D7$frac_res))))
  
  # save the simulated data: 
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  save(f_MIR_3D7,f_MIR_3D7_summary,file=glue::glue("output/Simulations-MIR-3D7-experiment_{Sys.Date()}_{cur_time}.RData"))
  
  # sensitivity analysis:
  if(sens_analysis_do){
    # sensitivity analysis
    tmp <- stoch_sim_exps_sens_analysis(exp_inoc=MIR_3D7_inocula,exp_n_inoc=MIR_3D7_n_inocula,exp_pmr=MIR_3D7_PMR,exp_n_gen=MIR_3D7_n_gen,n_sim,
                                        param_standard=standard_param,sens_list=sens_analysis_list)
    
    # save sensitiviy analysis results:
    MIR_3D7_sens_analysis_all <- tmp$all
    MIR_3D7_sens_analysis_summary <- tmp$summary
    MIR_3D7_sens_analysis_pos <- tmp$frac_pos
    
    cur_time <- format(Sys.time(),"h%H_m%M_s%S")
    save(MIR_3D7_sens_analysis_all,MIR_3D7_sens_analysis_summary,MIR_3D7_sens_analysis_pos,
         file=glue::glue("output/MIR-3D7-exp-sensitivity-analysis_{Sys.Date()}_{cur_time}.RData"))
  }
  
  # cleanup:
  rm(MIR_3D7_inocula,MIR_3D7_n_inocula,MIR_3D7_PMR,MIR_3D7_n_gen,f_MIR_3D7,tmp,MIR_3D7_sens_analysis_all,cur_time)
}


# simulation of field data ------------------------------------------------

if("LD_field"%in%sim_do | any(sim_do=="all")){
  # parameter values for the LD_field experiment:
  LD_field_inocula <- data_settings$inoculum[data_settings$data=="Laurent field"]
  LD_field_distr <- data_settings$n_cultures[data_settings$data=="Laurent field"]
  LD_field_n_gen <- data_settings$generations[data_settings$data=="Laurent field"]
  LD_field_PMR <- as.numeric(data_settings$PMR[data_settings$data=="Laurent field"])
  
  # make simulation:
  f_LD_field <- stoch_sim_exps(eq,LD_field_inocula,LD_field_distr,p=LD_field_PMR,n_gen=LD_field_n_gen,n_sim)
  
  # make summary data from the simulation:
  f_LD_field_summary <- data.frame(median=median(f_LD_field$frac_res),CI_lower=quantile(f_LD_field$frac_res,probs=(1-alpha_conf)/2),
                                   CI_upper=quantile(f_LD_field$frac_res,probs=(1+alpha_conf)/2),mean=mean(f_LD_field$frac_res),
                                   geom_mean=exp(mean(log(f_LD_field$frac_res))))
  
  # save the simulated data: 
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  save(f_LD_field,f_LD_field_summary,file=glue::glue("output/Simulations-LD-field-experiment_{Sys.Date()}_{cur_time}.RData"))
  
  # sensitivity analysis:
  if(sens_analysis_do){
    # sensitivity analysis
    tmp <- stoch_sim_exps_sens_analysis(exp_inoc=LD_field_inocula,exp_n_inoc=LD_field_distr,exp_pmr=LD_field_PMR,exp_n_gen=LD_field_n_gen,n_sim,
                                        param_standard=standard_param,sens_list=sens_analysis_list)
    
    # save sensitiviy analysis results:
    LD_field_sens_analysis_all <- tmp$all
    LD_field_sens_analysis_summary <- tmp$summary
    LD_field_sens_analysis_pos <- tmp$frac_pos
    
    cur_time <- format(Sys.time(),"h%H_m%M_s%S")
    save(LD_field_sens_analysis_all,LD_field_sens_analysis_summary,LD_field_sens_analysis_pos,
         file=glue::glue("output/LD-field-exp-sensitivity-analysis_{Sys.Date()}_{cur_time}.RData"))
  }
  
  # cleanup:
  rm(LD_field_inocula,LD_field_distr,LD_field_PMR,LD_field_n_gen,f_LD_field,tmp,LD_field_sens_analysis_all,cur_time)
}


# simulation of Laurent's 3D7 data ----------------------------------------

if("LD_3D7"%in%sim_do | any(sim_do=="all")){
  # parameter values for the LD_3D7 experiment:
  LD_3D7_inoculum <- as.numeric(data_settings$inoculum[data_settings$data=="Laurent 3D7"])
  LD_3D7_PMR <- as.numeric(data_settings$PMR[data_settings$data=="Laurent 3D7"])
  LD_3D7_n_gen <- as.numeric(data_settings$generations[data_settings$data=="Laurent 3D7"])
  LD_3D7_distr <- data_settings$n_cultures[data_settings$data=="Laurent 3D7"] # 52 cultures, all with the mean inoculation size from the data
  
  # make simulation:
  f_LD_3D7 <- stoch_sim_exps(eq,LD_3D7_inoculum,n_inoc=LD_3D7_distr,p=LD_3D7_PMR,n_gen=LD_3D7_n_gen,n_sim)
  
  # make summary data from the simulation:
  f_LD_3D7_summary <- data.frame(median=median(f_LD_3D7$frac_res),CI_lower=quantile(f_LD_3D7$frac_res,probs=(1-alpha_conf)/2),
                                CI_upper=quantile(f_LD_3D7$frac_res,probs=(1+alpha_conf)/2),mean=mean(f_LD_3D7$frac_res),
                                geom_mean=exp(mean(log(f_LD_3D7$frac_res))))
  
  # save the simulated data: 
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  save(f_LD_3D7,f_LD_3D7_summary,file=glue::glue("output/Simulations-LD-3D7-experiment_{Sys.Date()}_{cur_time}.RData"))
  
  # sensitivity analysis:
  if(sens_analysis_do){
    # sensitivity analysis
    tmp <- stoch_sim_exps_sens_analysis(exp_inoc=LD_3D7_inoculum,exp_n_inoc=LD_3D7_distr,exp_pmr=LD_3D7_PMR,exp_n_gen=LD_3D7_n_gen,n_sim,
                                        param_standard=standard_param,sens_list=sens_analysis_list)
    
    # save sensitiviy analysis results:
    LD_3D7_sens_analysis_all <- tmp$all
    LD_3D7_sens_analysis_summary <- tmp$summary
    LD_3D7_sens_analysis_pos <- tmp$frac_pos
    
    cur_time <- format(Sys.time(),"h%H_m%M_s%S")
    save(LD_3D7_sens_analysis_all,LD_3D7_sens_analysis_summary,LD_3D7_sens_analysis_pos,
         file=glue::glue("output/LD-3D7-exp-sensitivity-analysis_{Sys.Date()}_{cur_time}.RData"))
  }
  # cleanup:
  rm(LD_3D7_inoculum,LD_3D7_distr,LD_3D7_PMR,LD_3D7_n_gen,f_LD_3D7,tmp,LD_3D7_sens_analysis_all,cur_time)
}


# clean up ----------------------------------------------------------------

rm(sim_do,n_sim,eq,data_settings)

