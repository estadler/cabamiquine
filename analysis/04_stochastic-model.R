# -------------------------------------------------------------------------
#' Model estimates for the frequency of resistant parasites
#' stochastic, synchronous generations model
#' 
# -------------------------------------------------------------------------


# parameters for the model ------------------------------------------------

# default paramters initialized in setup.R, and:
init_sus_par <- 1 # initial number of susceptible parasites
n_gen <- 100 # number of generations to visualize
n_sim <- 1e5 # number of simulations


# make data for visualization ---------------------------------------------

sim_frac <- stoch_mod_rep_sim(n_gen,n_sim)
names(sim_frac) <- c("generation",glue::glue("frac_{c(1:n_sim)}"))

# replace simulations that gave NAs with other simulation:
while(any(colSums(is.na(sim_frac))>1)){
  sim_frac[,which(colSums(is.na(sim_frac))>1)[1]] <- stoch_mod_frac_res_sim(n_gen)$frac_res
}

# save the simulated data:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(sim_frac,file=glue::glue("output/Simulations-stochastic-model_{Sys.Date()}_{cur_time}.RData"))


# visualization -----------------------------------------------------------

# data for visualization:
vis_sim <- melt(sim_frac,id="generation",measure.vars = names(sim_frac)[2:ncol(sim_frac)],
                value.name = "fraction",variable.name = "simulation")
# make mean and median of stochastic simulations for visualization:
sim_summary <- data.frame(generation=sim_frac$generation,mean=row_geom_mean(sim_frac[,c(2:ncol(sim_frac))]),
                          median=matrixStats::rowMedians(as.matrix(sim_frac[,c(2:ncol(sim_frac))])))
# add the deterministic model to the plot:
det_model <- data.frame(generation=sim_frac$generation,model=det_mod_frac_res_gen(sim_frac$generation))
det_mod_frac_res_equil <- det_mod_equil()

# make the plot:
p_stoch <- ggplot(vis_sim,aes(x=generation,y=fraction,group=simulation)) +
  geom_line(color="grey") +
  # add median:
  geom_line(data=sim_summary,inherit.aes = FALSE,aes(x=generation,y=median),size=1,color="black") + 
  # add deterministic model:
  geom_line(data=det_model,inherit.aes = FALSE,aes(x=generation,y=model),size=1,color="dodgerblue") +
  geom_hline(yintercept=det_mod_frac_res_equil,color="red") +
  # theme, axes etc:
  scale_y_log10() +
  labs(x="Generation", y="Fraction of resistant parasites",title = "Fraction of resistant parasites by generation") +
  theme_bw() + theme(legend.key.size = unit(0.4, 'cm'),
                     panel.grid = element_line(colour="gray95",size = 0.2))

p_stoch

# save the plot as pdf and png:
ggsave(glue::glue("output/Stochastic-model-plot_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width, plot = p_stoch)
ggsave(glue::glue("output/Stochastic-model-plot_{Sys.Date()}_{cur_time}.png"),height=plot_height,width=plot_width, plot = p_stoch)


# make summary of simulated data ------------------------------------------

sim_frac_summary <- data.frame(generation=sim_frac$generation,geom_mean=row_geom_mean(sim_frac[,c(2:ncol(sim_frac))]),
                               mean=rowMeans(sim_frac[,c(2:ncol(sim_frac))]),median=matrixStats::rowMedians(as.matrix(sim_frac[,c(2:ncol(sim_frac))])),
                               CI_lower=apply(sim_frac[,c(2:ncol(sim_frac))],1,function(x){quantile(x,probs=c((1-alpha_conf)/2))}),
                               CI_upper=apply(sim_frac[,c(2:ncol(sim_frac))],1,function(x){quantile(x,probs=c((1+alpha_conf)/2))}))

# save the simulated data:
save(sim_frac_summary,file=glue::glue("output/Simulations-stochastic-model-summary_{Sys.Date()}_{cur_time}.RData"))


# visualize summary data with CI and deterministic model ------------------

# make plot data:
p_data <- cbind(sim_frac_summary,deterministic=det_model$model,equil=rep(det_mod_frac_res_equil,nrow(sim_frac_summary)))
p_data <- melt(p_data[,names(p_data)%in%c("generation","median","deterministic","equil")],id="generation")

p_summary <- ggplot(p_data,aes(x=generation,y=value,color=variable))+
  # confidence region:
  geom_ribbon(data=sim_frac_summary,inherit.aes = FALSE,aes(x=generation,ymin=CI_lower,ymax=CI_upper),fill="black", alpha = 0.1)+
  # stochastic model:
  geom_line(size=1) +
  # legend:
  scale_color_manual(values=c("median"="black","deterministic"="dodgerblue","equil"="red"),
                     labels=c("median of stochastic simulations","deterministic model","equilibrium"),name="") +
  # theme, axes etc:
  scale_y_log10() +
  labs(x="Generation", y="Fraction of resistant parasites",title = "Fraction of resistant parasites by generation") +
  theme_bw() + theme(legend.key.size = unit(0.4, 'cm'),legend.position="bottom",
                     panel.grid = element_line(colour="gray95",size = 0.2))

p_summary

# save the plot:
ggsave(glue::glue("output/Stochastic-model-plot-CI_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width, plot = p_summary)
ggsave(glue::glue("output/Stochastic-model-plot-CI_{Sys.Date()}_{cur_time}.png"),height=plot_height,width=plot_width, plot = p_summary)


# cleanup -----------------------------------------------------------------

rm(init_sus_par,n_gen,n_sim,sim_frac,cur_time,vis_sim,sim_summary,det_model,det_mod_frac_res_equil,p_stoch,
   sim_frac_summary,p_data,p_summary)

