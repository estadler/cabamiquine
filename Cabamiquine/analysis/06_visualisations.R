# -------------------------------------------------------------------------
#' Visualisation of the data estimates and the estimates from the 
#' stochastic model simulations
#' 
# -------------------------------------------------------------------------

# run first the estimation of resistant parasites from the data ("analysis/01_estimates-from-data.R") 
# and the stochastic simulations of the different experiments ("analysis/05_stochastic-sim-experiments.R"),
# then visualise the estimates of the frequency of resistant parasites.


# visualization of data estimates and stochastic simulations -------------

# make a data frame for visualisation:
data_vis <- data.frame(data=rep(c("MIR_Dd2","MIR_3D7","LD_3D7","LD_field","NSG","VIS"),2),
                       group=c(rep("data",6),rep("sim",6)),
                       freq=c(data.estimates$freq_per_par[data.estimates$group=="MIR_Dd2"],
                              data.estimates$freq_per_par[data.estimates$group=="MIR_3D7"],
                              data.estimates$freq_per_par[data.estimates$group=="ld_3D7_x52"],
                              data.estimates$freq_per_par[data.estimates$group=="field"],
                              data.estimates$freq_per_par[data.estimates$group=="NSG"],
                              data.estimates$freq_per_par[data.estimates$group=="VIS"],
                              1/f_MIR_Dd2_summary$median,1/f_MIR_3D7_summary$median,1/f_LD_3D7_summary$median,
                              1/f_LD_field_summary$median,1/f_NSG_summary$median,1/f_VIS_summary$median),
                       CI_lower=c(data.estimates$freq_CI_lower[data.estimates$group=="MIR_Dd2"],
                                  data.estimates$freq_CI_lower[data.estimates$group=="MIR_3D7"],
                                  data.estimates$freq_CI_lower[data.estimates$group=="ld_3D7_x52"],
                                  data.estimates$freq_CI_lower[data.estimates$group=="field"],
                                  data.estimates$freq_CI_lower[data.estimates$group=="NSG"],
                                  data.estimates$freq_CI_lower[data.estimates$group=="VIS"],
                                  1/f_MIR_Dd2_summary$CI_upper,1/f_MIR_3D7_summary$CI_upper,1/f_LD_3D7_summary$CI_upper,
                                  1/f_LD_field_summary$CI_upper,1/f_NSG_summary$CI_upper,1/f_VIS_summary$CI_upper),
                       CI_upper=c(data.estimates$freq_CI_upper[data.estimates$group=="MIR_Dd2"],
                                  data.estimates$freq_CI_upper[data.estimates$group=="MIR_3D7"],
                                  data.estimates$freq_CI_upper[data.estimates$group=="ld_3D7_x52"],
                                  data.estimates$freq_CI_upper[data.estimates$group=="field"],
                                  data.estimates$freq_CI_upper[data.estimates$group=="NSG"],
                                  data.estimates$freq_CI_upper[data.estimates$group=="VIS"],
                                  1/f_MIR_Dd2_summary$CI_lower,1/f_MIR_3D7_summary$CI_lower,1/f_LD_3D7_summary$CI_lower,
                                  1/f_LD_field_summary$CI_lower,1/f_NSG_summary$CI_lower,1/f_VIS_summary$CI_lower))
data_vis$data <- factor(data_vis$data,levels=c("VIS","NSG","LD_field","LD_3D7","MIR_3D7","MIR_Dd2"))
data_vis <- cbind(data_vis,y=as.numeric(data_vis$data))
data_vis$y[data_vis$group=="data"] <- data_vis$y[data_vis$group=="data"]
data_vis$y[data_vis$group=="sim"] <- data_vis$y[data_vis$group=="sim"]-0.1

# visualisation:
my_colours <- c("dodgerblue","dodgerblue","darkorange","darkorange","black","black")
eq <- 1/det_mod_equil() # equilibrium number of parasites to obtain 1 resistant mutant

p_freq <- ggplot(data_vis,aes(x=freq,y=y,shape=group,color=data)) +
  # indicate equilibrium frequency of resistant mutants:
  geom_vline(xintercept=eq,color="red",alpha=0.5) +
  annotate("text",x=eq,y=3.5,angle=90,label="equilibrium frequency\nof resistant parasites\n\n",color="red",alpha=0.5) +
  # visualise data and simulation estimates:
  geom_point() +
  geom_errorbarh(aes(xmin=CI_lower,xmax=CI_upper,linetype=group),height=0) +
  # set colors and shapes, modify legends:
  scale_color_manual(values=my_colours) + 
  scale_shape_manual(values=c(16,15),name=NULL,labels=c("data","stochastic model simulations")) +
  scale_linetype_manual(values=c("solid","dashed"),name=NULL,labels=c("data","stochastic model simulations")) +
  guides(color="none") +
  # theme, axes etc:
  scale_x_log10(breaks=scales::trans_breaks("log10", function(x) 10^x),
                labels=scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "b", outside = TRUE) +   coord_cartesian(clip = "off",xlim = c(4e6,2e10)) +
  scale_y_discrete(breaks=c(1:6),limits = factor(c(1:6)),
                   labels=c(expression(atop("Regrowth"~italic("in vivo")~"VIS","(3D7/Naive individuals)")),
                            expression(atop("Regrowth"~italic("in vivo")~"NSG","(3D7/Naive donor cells)")),
                            expression(atop("Regrowth"~italic("in vitro"),"(Field/Endemic donor cells)")),
                            expression(atop("Regrowth"~italic("in vitro"),"(3D7/Endemic donor cells)")),
                            expression(atop("MIR"~italic("in vitro"),"(3D7/Naive donor cells)")),
                            expression(atop("MIR"~italic("in vitro"),"(Dd2/Naive donor cells)")))) +
  labs(x="Number of parasites to obtain 1 resistant mutant", y="",title = "Estimated frequency of resistant mutants") +
  theme_bw() + theme(legend.key.size = unit(0.4, 'cm'),legend.position="bottom",panel.grid.minor = element_blank(),
                     panel.grid = element_line(colour="gray95",size = 0.2))

p_freq

# save plot:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Frequency-estimates-experiment_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width, plot = p_freq)

# clean up: 
rm(data_vis,my_colours,p_freq,cur_time)


# visualisation of sensitivity analysis results ---------------------------

if(sens_analysis_do){ # only visualize the sensitivity analysis if the sensitivity analysis was done
  # parameters & values for the sensitivity analysis in sens_analysis_list
  for(i in 1:length(sens_analysis_list)){
    # make a data frame for visualization
    freq_sens_analysis <- c(sapply(c(1:length(sens_analysis_list[[i]])),function(x){c(MIR_Dd2_sens_analysis_summary[[i]]$median[x],MIR_3D7_sens_analysis_summary[[i]]$median[x],LD_3D7_sens_analysis_summary[[i]]$median[x],
                                                                                      LD_field_sens_analysis_summary[[i]]$median[x],NSG_sens_analysis_summary[[i]]$median[x],VIS_sens_analysis_summary[[i]]$median[x])}))
    CI_lower_sens_analysis <- c(sapply(c(1:length(sens_analysis_list[[i]])),function(x){c(MIR_Dd2_sens_analysis_summary[[i]]$CI_lower[x],MIR_3D7_sens_analysis_summary[[i]]$CI_lower[x],LD_3D7_sens_analysis_summary[[i]]$CI_lower[x],
                                                                                          LD_field_sens_analysis_summary[[i]]$CI_lower[x],NSG_sens_analysis_summary[[i]]$CI_lower[x],VIS_sens_analysis_summary[[i]]$CI_lower[x])}))
    CI_upper_sens_analysis <- c(sapply(c(1:length(sens_analysis_list[[i]])),function(x){c(MIR_Dd2_sens_analysis_summary[[i]]$CI_upper[x],MIR_3D7_sens_analysis_summary[[i]]$CI_upper[x],LD_3D7_sens_analysis_summary[[i]]$CI_upper[x],
                                                                                          LD_field_sens_analysis_summary[[i]]$CI_upper[x],NSG_sens_analysis_summary[[i]]$CI_upper[x],VIS_sens_analysis_summary[[i]]$CI_upper[x])}))
    data_vis_tmp <- data.frame(data=rep(c("MIR_Dd2","MIR_3D7","LD_3D7","LD_field","NSG","VIS"),1+length(sens_analysis_list[[i]])),
                               group=c(rep("data",6),rep("sim",6*length(sens_analysis_list[[i]]))),
                               value=unlist(c(rep(standard_param[names(standard_param)==names(sens_analysis_list)[i]],6),rep(sens_analysis_list[[i]],rep(6,length(sens_analysis_list[[i]]))))),
                               freq=c(data.estimates$freq_per_par[data.estimates$group=="MIR_Dd2"],data.estimates$freq_per_par[data.estimates$group=="MIR_3D7"],
                                      data.estimates$freq_per_par[data.estimates$group=="ld_3D7_x52"],data.estimates$freq_per_par[data.estimates$group=="field"],
                                      data.estimates$freq_per_par[data.estimates$group=="NSG"],data.estimates$freq_per_par[data.estimates$group=="VIS"],
                                      1/freq_sens_analysis),
                               CI_lower=c(data.estimates$freq_CI_lower[data.estimates$group=="MIR_Dd2"],data.estimates$freq_CI_lower[data.estimates$group=="MIR_3D7"],
                                          data.estimates$freq_CI_lower[data.estimates$group=="ld_3D7_x52"],data.estimates$freq_CI_lower[data.estimates$group=="field"],
                                          data.estimates$freq_CI_lower[data.estimates$group=="NSG"],data.estimates$freq_CI_lower[data.estimates$group=="VIS"],
                                          1/CI_upper_sens_analysis),
                               CI_upper=c(data.estimates$freq_CI_upper[data.estimates$group=="MIR_Dd2"],data.estimates$freq_CI_upper[data.estimates$group=="MIR_3D7"],
                                          data.estimates$freq_CI_upper[data.estimates$group=="ld_3D7_x52"],data.estimates$freq_CI_upper[data.estimates$group=="field"],
                                          data.estimates$freq_CI_upper[data.estimates$group=="NSG"],data.estimates$freq_CI_upper[data.estimates$group=="VIS"],
                                          1/CI_lower_sens_analysis))
    data_vis_tmp$data <- factor(data_vis_tmp$data,levels=c("VIS","NSG","LD_field","LD_3D7","MIR_3D7","MIR_Dd2"))
    data_vis_tmp <- cbind(data_vis_tmp,y=as.numeric(data_vis_tmp$data))
    for(j in 1:length(sens_analysis_list[[i]])){
      data_vis_tmp$y[which(data_vis_tmp$value==sort(unique(data_vis_tmp$value))[j] & data_vis_tmp$group=="sim")] <- data_vis_tmp$y[which(data_vis_tmp$value==sort(unique(data_vis_tmp$value))[j] & data_vis_tmp$group=="sim")]-0.1*j
    }
    if(names(sens_analysis_list)[i]=="fitn_cost"){
      data_vis_tmp$value <- 100*(1-data_vis_tmp$val)
    }
    
    # visualize: 
    my_colours <- c("darkorange","orchid2","palegreen2","springgreen4","dodgerblue")
    my_colours[which(sens_analysis_list[[i]]==as.numeric(standard_param[names(standard_param)==names(sens_analysis_list)[i]]))] <- "black" # data estimate and "standard parameter" values in black
    colour_name <- ifelse(names(sens_analysis_list)[i]=="mut_rate","Mutation rate:",
                          ifelse(names(sens_analysis_list)[i]=="n_res_mut","Number of resistance mutations:",
                                 ifelse(names(sens_analysis_list)[i]=="fitn_cost","Fitness cost [%]:",
                                        ifelse(names(sens_analysis_list)[i]=="pmr","PMR:",names(sens_analysis_list)[i]))))
    
    p_sens <- ggplot(data_vis_tmp,aes(x=freq,y=y,shape=group,color=as.factor(value))) +
      # visualise data and simulation estimates:
      geom_point() +
      geom_errorbarh(aes(xmin=CI_lower,xmax=CI_upper,linetype=group),height=0) +
      # set colors and shapes, modify legends:
      scale_color_manual(values=my_colours,name=colour_name) +
      scale_shape_manual(values=c(16,15),name=NULL,labels=c("data","stochastic model simulations")) +
      scale_linetype_manual(values=c("solid","dashed"),name=NULL,labels=c("data","stochastic model simulations")) +
      # theme, axes etc:
      scale_x_log10(breaks=scales::trans_breaks("log10", function(x) 10^x),
                    labels=scales::trans_format("log10", scales::math_format(10^.x))) +
      annotation_logticks(sides = "b", outside = TRUE) + 
      coord_cartesian(clip = "off", xlim=c(min(c(data_vis_tmp$freq,4e6)),max(c(data_vis_tmp$freq,2e10)))) + 
      scale_y_discrete(breaks=c(1:6),limits = factor(c(1:6)),
                       labels=c(expression(atop("Regrowth"~italic("in vivo")~"VIS","(3D7/Naive individuals)")),
                                expression(atop("Regrowth"~italic("in vivo")~"NSG","(3D7/Naive donor cells)")),
                                expression(atop("Regrowth"~italic("in vitro"),"(Field/Endemic donor cells)")),
                                expression(atop("Regrowth"~italic("in vitro"),"(3D7/Endemic donor cells)")),
                                expression(atop("MIR"~italic("in vitro"),"(3D7/Naive donor cells)")),
                                expression(atop("MIR"~italic("in vitro"),"(Dd2/Naive donor cells)")))) +
      labs(x="Number of parasites to obtain 1 resistant mutant", y="",title = "Estimated frequency of resistant mutants") +
      theme_bw() + theme(legend.key.size = unit(0.4, 'cm'),legend.position="bottom",panel.grid.minor = element_blank(),
                         panel.grid = element_line(colour="gray95",size = 0.2))
    
    p_sens
    
    # save plot:
    cur_time <- format(Sys.time(),"h%H_m%M_s%S")
    ggsave(glue::glue("output/Sensitivity-analysis_{names(sens_analysis_list)[i]}_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width, plot = p_sens)
  }
}

# clean up: 
rm(i,freq_sens_analysis,CI_lower_sens_analysis,CI_upper_sens_analysis,data_vis_tmp,j,my_colours,colour_name,p_sens,cur_time)


# probability of emergence only after treatment ---------------------------

# probability for emergence of resistant mutants only after treatment (number of parasites on the x-axis):
prr <- 1/c(5,10,50,100,500,1000) # parasite reduction ratio
num_par <- 10^seq(5,11,length.out = 1e3) # number of parasites at treatent possible range

# probability of emerging after treatment given that there are no resistant parasites at treatment:
fun_prob_after <- Vectorize(function(par_0){1-dbinom(0,ceiling(par_0*prr/(1-prr)),mut_rate*n_res_mut)})
prob_after <- 100*fun_prob_after(num_par)

# probability of resistant parasites being present at treatment:
frac_before <- data_est_all # fraction of mutants before treatment estimated from all data (2.325211e-09)
fun_prob_before <- Vectorize(function(n){1-dbinom(0,ceiling(n),frac_before)}) # probability of presence before treatment
prob_before <- 100*fun_prob_before(num_par)

# make a data frame for visualisation:
data_vis <- data.frame(par_treatment=rep(num_par,nrow(prob_after)),
                       prob=as.vector(sapply(c(1:nrow(prob_after)),function(x){(100-prob_before)/100*prob_after[x,]})),
                       prr=rep(prr,rep(length(num_par),length(prr))))
data_before <- data.frame(x=num_par,y=prob_before)

my_colours <- c(rev(viridis(dim(prob_after)[1])),"red")

p_emergence <- ggplot(data_vis,aes(x=par_treatment,y=prob,color=factor(prr))) +
  geom_line() + 
  # add probability of presence of resistant mutants at the time of treatment:
  geom_line(data=data_before,inherit.aes = FALSE,aes(x,y),color="red",linetype="dashed") +
  annotate("text",x=1e8,y=80,angle=0,label="Probability of presence\nof resistant mutants at\n time of treatment",color="red") +
  # set colors and shapes, modify legends:
  scale_color_manual(values=my_colours, name="Parasite Reduction Ratio (fold reduction in viable parasites per cycle): ",labels=1/sort(unique(data_vis$prr))) + 
  # theme, axes etc:
  scale_x_log10(breaks=scales::trans_breaks("log10", function(x) 10^x),
                labels=scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "b", outside = TRUE) + 
  coord_cartesian(clip = "off", xlim=c(1e5,1e10), ylim=c(-2,102), expand = FALSE) + 
  labs(x="Parasite number at treatment", y="Probability [%]",title = "Probability of emergence only after treatment") +
  theme_bw() + theme(legend.key.size = unit(0.4, 'cm'),legend.position="bottom",panel.grid.minor = element_blank(),
                     panel.grid = element_line(colour="gray95",size = 0.2))

p_emergence

# save plot:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Probability_emergence_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width, plot = p_emergence)

# clean up:
rm(prr,num_par,fun_prob_after,prob_after,frac_before,fun_prob_before,prob_before,data_vis,data_before,my_colours,p_emergence,cur_time)

