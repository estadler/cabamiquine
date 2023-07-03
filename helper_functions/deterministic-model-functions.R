
# deterministic model -----------------------------------------------------

# susceptible parasites in a certain generation of the deterministic model:
det_mod_sus_gen <- function(gen){
  init_sus_par*(pmr-n_res_mut*mut_rate)^gen
}

# resistant parasites in a certain generation of the deterministic model:
det_mod_res_gen <- function(gen){
  init_sus_par*(n_res_mut*mut_rate*((pmr-n_res_mut*mut_rate)^gen-(pmr*fitn_cost)^gen))/(pmr-n_res_mut*mut_rate-pmr*fitn_cost)+
    init_res_par*(pmr*fitn_cost)^gen
}

# fraction of resistant parasites by generation in the deterministic model:
det_mod_frac_res_gen <- function(gen){
  det_mod_res_gen(gen)/(det_mod_res_gen(gen)+det_mod_sus_gen(gen))
}

# equilibrium fraction of resistant mutants:
det_mod_equil <- function(){
  (n_res_mut*mut_rate)/(pmr-pmr*fitn_cost)
}

# equilibrium fraction of resistant mutants with parameters:
det_mod_equil_par <- function(n_res,mut,fitness,p){
  (n_res*mut)/(p-p*fitness)
}
