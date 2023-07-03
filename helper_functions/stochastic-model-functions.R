
# stochastic model --------------------------------------------------------

# simulate next generation of the stochastic model:
stoch_mod_next_gen <- function(n_sus_par,n_res_par,p=pmr,fitness=fitn_cost,mut=mut_rate,n_res=n_res_mut){
  # using that the sum of independent Poisson distributions is a Poisson distribution 
  # (parameter=sum of the parameters of the independent Poisson distributions):
  offspring_sus <- as.numeric(rpois(n=1,lambda=n_sus_par*p))
  offspring_res <- as.numeric(rpois(n=1,lambda=n_res_par*p*fitness))
  
  # mutation: 
  mut_sus_to_res <- as.numeric(rbinom(n=1,size=n_sus_par,prob=mut*n_res))
  
  c(offspring_sus-mut_sus_to_res,offspring_res+mut_sus_to_res)
}

# simulate the fraction of resistant parasites from generation 0 to n_gen:
stoch_mod_frac_res_sim <- function(n_gen){
  tmp_res <- init_res_par
  tmp_sus <- init_sus_par
  tmp_frac <- c(tmp_res/(tmp_res+tmp_sus),rep(NA,n_gen))
  for(i in 2:(1+n_gen)){
    tmp <- stoch_mod_next_gen(tmp_sus,tmp_res)
    tmp_sus <- tmp[1]
    tmp_res <- tmp[2]
    tmp_frac[i] <- tmp_res/(tmp_res+tmp_sus)
  }
  
  data.frame(generation=c(0:n_gen),frac_res=tmp_frac)
}

# simulate the fraction of resistant parasites n_sim times for n_gen generations:
stoch_mod_rep_sim <- function(n_gen,n_sim){
  tmp <- stoch_mod_frac_res_sim(n_gen)
  
  for(i in 2:n_sim){
    tmp2 <- stoch_mod_frac_res_sim(n_gen)
    tmp <- cbind(tmp,tmp2$frac_res)
  }
  
  tmp
}


# Stochastic simulations of the experiments -------------------------------

# simulate one experiment and calculate the estimated frequency of resistant mutants in this experiment:
stoch_sim_1_exp <- function(eq=4.125e-08,inoc,n_inoc,p=pmr,n_gen,
                            fitness=fitn_cost,mut=mut_rate,n_res=n_res_mut){
  # input parameters:
  # eq          equilibrium distribution of resistant parasites in the stock
  # inoc        inocula (vector of the different inocula in the experiment)
  # n_inoc      numbers of different inocula (same length as inoc)
  # p           PMR (same length as inocula)
  # n_gen       number of generations from inoculation to treatment (same length as inocula)
  # fitness     fitness cost for resistance
  # mut         mutation rate
  # n_res       number of resistance mutations
  
  # sample from the parasite stock:
  tmp_res <- unlist(sapply(c(1:length(inoc)),function(x){rbinom(n=n_inoc[x],size=inoc[x],prob=eq)}))
  tmp_sus <- rep(inoc,n_inoc)-tmp_res
  tmp_p <- rep(p,n_inoc)
  
  # simulate each experiment for n_gen generations:
  if(length(n_gen)==1){
    if(n_gen>0){
      for(i in 1:n_gen){
        for(j in 1:length(tmp_res)){
          # next generation:
          offspring_sus <- as.numeric(rpois(n=1,lambda=tmp_sus[j]*tmp_p[j]))
          offspring_res <- as.numeric(rpois(n=1,lambda=tmp_res[j]*tmp_p[j]*fitness))
          
          # mutation: 
          mut_sus_to_res <- as.numeric(rbinom(n=1,size=tmp_sus[j],prob=mut*n_res))
          
          tmp_sus[j] <- offspring_sus-mut_sus_to_res
          tmp_res[j] <- offspring_res+mut_sus_to_res
        }
      }
    }
  }else{
    tmp_n_gen <- rep(n_gen,n_inoc)
    for(j in 1:length(tmp_res)){
      if(tmp_n_gen[j]>0){
        for(i in 1:tmp_n_gen[j]){
          # next generation:
          offspring_sus <- as.numeric(rpois(n=1,lambda=tmp_sus[j]*tmp_p[j]))
          offspring_res <- as.numeric(rpois(n=1,lambda=tmp_res[j]*tmp_p[j]*fitness))
          
          # mutation: 
          mut_sus_to_res <- as.numeric(rbinom(n=1,size=tmp_sus[j],prob=mut*n_res))
          
          tmp_sus[j] <- offspring_sus-mut_sus_to_res
          tmp_res[j] <- offspring_res+mut_sus_to_res
        }
      }
    }
  }
  
  # make experiment data:
  inoculum <- sort(unique(tmp_res+tmp_sus)) # numbers of parasites at the time of treatment
  n_wells <- apply(as.data.frame(inoculum),1,function(x){sum((tmp_res+tmp_sus)==x)})
  n_pos_wells <- apply(as.data.frame(inoculum),1,function(x){sum(tmp_res[tmp_res+tmp_sus==x]>=1)})
  sim_exp_data <- data.frame(inoculum,n_wells,n_pos_wells)
  
  # save the fraction of positive wells in the experiment:
  tmp_pos_wells <- sum(n_pos_wells)/sum(n_wells)
  
  # estimate the frequency of resistant mutants in this experiment as for the data:
  init.freq <- 2e-8 # initial frequency of resistant mutants guess
  nllh_tmp <- function(log_f){nllh(log_f,sim_exp_data)}
  fit_tmp <- nlm(nllh_tmp,log(init.freq),hessian=TRUE)
  
  # return the estimated frequency of resistant mutants and the fraction of positive wells in the simulated experiment:
  tmp_f <- exp(fit_tmp$estimate)
  
  list(frac_res=tmp_f,frac_pos_wells=tmp_pos_wells)
}

# simulate several experiments and calculate the estimated frequency of resistant mutants in this experiment:
stoch_sim_exps <- function(eq=4.125e-08,inoc,n_inoc,p=pmr,n_gen,n_sim,fitness=fitn_cost,mut=mut_rate,n_res=n_res_mut){
  # input parameters:
  # eq          equilibrium distribution of resistant parasites in the stock
  # inoc        inocula
  # n_inoc      numbers of different inocula (same length as inoc)
  # p           PMR
  # n_gen       number of generations from inoculation to treatment
  # n_sim       number of simulations of experiments
  # fitness     fitness cost for resistance
  # mut         mutation rate
  # n_res       number of resistance mutations
  
  tmp_frac_res <- rep(NA,n_sim)
  tmp_frac_pos_wells <- rep(NA,n_sim)
  for(i in 1:n_sim){
    tmp_exp <- stoch_sim_1_exp(eq,inoc,n_inoc,p,n_gen,fitness,mut,n_res)
    tmp_frac_res[i] <- tmp_exp$frac_res
    tmp_frac_pos_wells[i] <- tmp_exp$frac_pos_wells
  }
  
  list(frac_res=tmp_frac_res,frac_pos_wells=tmp_frac_pos_wells)
}

# sensitivity analysis of the simulations of experiments for different parameter values:
stoch_sim_exps_sens_analysis <- function(exp_inoc,exp_n_inoc,exp_pmr,exp_n_gen,n_sim,param_standard,sens_list){
  # input parameters:
  # exp_inoc        inocula for the experiment
  # exp_n_inoc      numbers of cultures for each inoculum for the experiment (same length as exp_inoc)
  # exp_pmr         PMR for the experiment (same length as exp_inoc)
  # exp_n_gen       number of generations for the experiment (same length as exp_inoc)
  # n_sim           number of simulations of experiments
  # param_standard  standard parameters: a list or dataframe with the parameter names and their "standard values"
  # sens_list       parameter values for the sensitivity analysis: named list with parameter values to try in the sensitivity analysis
  
  # to save sensitivity analysis data:
  sens_analysis_all <- vector(mode='list',length=length(sens_list)) # to save the estimated fraction of resistant parasites in each experiment
  names(sens_analysis_all) <- names(sens_list)
  sens_analysis_pos <- vector(mode='list',length=length(sens_list)) # to save the fraction of positive wells in each experiment
  names(sens_analysis_pos) <- names(sens_list)
  sens_analysis_summary <- vector(mode='list',length=length(sens_list)) # to save the summary of the estimated fraction of resistant parasites in each experiment
  names(sens_analysis_summary) <- names(sens_list)
  
  # sensitivity analysis:
  for(i in 1:length(sens_list)){ # sensitivity analysis for each parameter in the list
    # to save simulation results:
    tmp_sim_result <- matrix(NA,n_sim,length(sens_list[[i]])) # estimated frequency for each simulated experiment
    tmp_sim_pos <- matrix(NA,n_sim,length(sens_list[[i]])) # fraction of positive wells in each simulated experiment
    tmp_sim_result_summary <- data.frame(value=c(),median=c(),CI_lower=c(),CI_upper=c(),mean=c(),geom_mean=c()) # summary of simulation for each parameter value set
    
    for(j in 1:length(sens_list[[i]])){ # for each of set of parameter values
      # current parameter values (standard values except for currently varied parameter):
      tmp_n_res_mut <- ifelse(names(sens_list)[i]=="n_res_mut",sens_analysis_list[[i]][j],param_standard$n_res_mut)
      tmp_mut_rate <- ifelse(names(sens_list)[i]=="mut_rate",sens_analysis_list[[i]][j],param_standard$mut_rate)
      tmp_fitn_cost <- ifelse(names(sens_list)[i]=="fitn_cost",sens_analysis_list[[i]][j],param_standard$fitn_cost)
      tmp_pmr <- ifelse(names(sens_list)[i]=="pmr",sens_analysis_list[[i]][j],param_standard$pmr)
      
      # compute the equilibrium frequency of resistant parasites for the current parameters:
      tmp_eq <- det_mod_equil_par(n_res=tmp_n_res_mut,mut=tmp_mut_rate,fitness=tmp_fitn_cost,p=tmp_pmr)
      
      # make a simulation for the current parameters:
      tmp_exps <- stoch_sim_exps(eq=tmp_eq,inoc=exp_inoc,n_inoc=exp_n_inoc,p=exp_pmr,n_gen=exp_n_gen,n_sim,
                                           fitness = tmp_fitn_cost, mut=tmp_mut_rate, n_res = tmp_n_res_mut)
      tmp_sim_result[,j] <- tmp_exps$frac_res
      tmp_sim_pos[,j] <- tmp_exps$frac_pos_wells
      
      # make summary of the simulation:
      tmp_sim_result_summary <- rbind(tmp_sim_result_summary,
                                      data.frame(value=sens_list[[i]][j],median=median(tmp_sim_result[,j]),
                                                 CI_lower=quantile(tmp_sim_result[,j],probs=(1-alpha_conf)/2),
                                                 CI_upper=quantile(tmp_sim_result[,j],probs=(1+alpha_conf)/2),
                                                 mean=mean(tmp_sim_result[,j]),geom_mean=exp(mean(log(tmp_sim_result[,j]))))) 
    }
    
    rownames(tmp_sim_result_summary) <- NULL
    
    # save simulation results:
    sens_analysis_all[[i]] <- tmp_sim_result
    sens_analysis_pos[[i]] <- tmp_sim_pos
    sens_analysis_summary[[i]] <- tmp_sim_result_summary
  }
  
  # return all simulation results and the summary of the simulation results:
  list(all=sens_analysis_all,summary=sens_analysis_summary,frac_pos=sens_analysis_pos)
}
