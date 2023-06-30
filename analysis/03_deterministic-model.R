# -------------------------------------------------------------------------
#' Model estimates for the frequency of resistant parasites
#' Deterministic model and equilibrium frequency
#' 
# -------------------------------------------------------------------------


# parameters for the model ------------------------------------------------

# default paramters initialized in setup.R, and:
init_sus_par <- 1 # initial number of susceptible parasites
det_mod_frac_res_equil <- det_mod_equil() # equilibrium of the deterministic model


# visualize the deterministic model ---------------------------------------

# parameters
n_gen <- 100 # number of generations to be visualized
det_data <- data.frame(generation=c(0:n_gen),model=det_mod_frac_res_gen(c(0:n_gen)),equil=rep(det_mod_frac_res_equil,n_gen+1))
p_det_data <- melt(det_data,id="generation")

p_det <- ggplot(p_det_data,aes(x=generation,y=value,color=variable)) +
  geom_line(size=1) + 
  # legend:
  scale_color_manual(values=c("model"="black","equil"="red"),labels=c("deterministic model","equilibrium"),name="") +
  # theme, axes etc:
  scale_y_log10() +
  labs(x="Generation", y="Fraction of resistant parasites",title = "Fraction of resistant parasites by generation") +
  theme_bw() + theme(legend.key.size = unit(0.4, 'cm'),legend.position="bottom",
                     panel.grid = element_line(colour="gray95",size = 0.2))

p_det

# save plot:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Deterministic-model_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width, plot = p_det)


# cleanup -----------------------------------------------------------------

rm(init_sus_par,det_mod_frac_res_equil,n_gen,det_data,p_det_data,p_det,cur_time)
