# -------------------------------------------------------------------------
#' Model estimates for the frequency of resistant parasites
#' 
# -------------------------------------------------------------------------


# pairwise comparison of the data estimates -------------------------------

# use a likelihood ratio test for a pairwise comparison the estimates from the data
# null hypothesis of the test: difference between the estimated frequencies of resistant mutants in different data sets is 0
# alternative hypothesis: there is a difference in the estimated frequency of resistant mutants

data.est.comp <- data.frame(data.set=sort(unique(data.estimates$group)))

for(j in c(2:nrow(data.est.comp))){
  # extend data frame:
  data.est.comp <- cbind(data.est.comp,tmp=rep(NA,nrow(data.est.comp)))
  names(data.est.comp)[which(names(data.est.comp)=="tmp")] <- data.est.comp$data.set[j]
  
  for(i in c(1:(j-1))){
    # estimating the frequency of resistant mutants in both data sets together:
    nllh_both <- function(log_f){nllh(log_f,data.all[data.all$group%in%c(data.est.comp$data.set[i],data.est.comp$data.set[j]),])}
    fit_both <- nlm(nllh_both,log(initial_freq_guess),hessian=TRUE)
    
    # estimating the frequency of resistant mutants with a difference between the estimates from the two data:
    nllh_diff <- function(par){
      # par[1]: log-transformed frequency of resistant mutants in the first data set
      # par[2]: exp(par[1])+par[2] is the estimated frequency of resistant mutants in the second data set
      nllh_d(par[1],0,data.all[data.all$group==data.est.comp$data.set[i],])+
        nllh_d(par[1],par[2]/1e10,data.all[data.all$group==data.est.comp$data.set[j],])}
    # note that the difference of the estimates (par[2]) is scaled by 1e10 to improve numerical behavior
    fit_diff <-  nlm(nllh_diff,c(fit_both$estimate,0),hessian=TRUE)
    
    # likelihood ratio test comparison of the estimates:
    p.val <- pchisq(2*(fit_both$minimum-fit_diff$minimum),df=1,lower.tail = FALSE)
    
    # save p-value:
    data.est.comp[i,j] <- p.val
  }
}

# remove last row that only contains NAs:
data.est.comp <- data.est.comp[c(1:(nrow(data.est.comp)-1)),]


# compare field & LD's 52x3D7 data by number of positive cultures ---------

# use the number of cultures and the number of recrudescences only, not the inoculum sizes 

# same underlying frequency of resistant mutants:
nllh_both <- function(log_f){-log(dbinom(data.all$n_pos_wells[data.all$group=="field_sum"],data.all$n_wells[data.all$group=="field_sum"],exp(log_f)))-
    log(dbinom(data.all$n_pos_wells[data.all$group=="ld_3D7_x52"],data.all$n_wells[data.all$group=="ld_3D7_x52"],exp(log_f)))}
fit_both <- nlm(nllh_both,log(sum(data.all$n_pos_wells[data.all$group%in%c("field_sum","ld_3D7_x52")])/sum(data.all$n_wells[data.all$group%in%c("field_sum","ld_3D7_x52")])),hessian=TRUE)

# different underlying frequencies:
nllh_diff <- function(par){-log(dbinom(data.all$n_pos_wells[data.all$group=="field_sum"],data.all$n_wells[data.all$group=="field_sum"],exp(par[1])))-
    log(dbinom(data.all$n_pos_wells[data.all$group=="ld_3D7_x52"],data.all$n_wells[data.all$group=="ld_3D7_x52"],exp(par[1])+par[2]))}
fit_diff <- nlm(nllh_diff,c(log(sum(data.all$n_pos_wells[data.all$group%in%c("field_sum","ld_3D7_x52")])/sum(data.all$n_wells[data.all$group%in%c("field_sum","ld_3D7_x52")])),0),hessian=TRUE)

# likelihood ratio test comparison of the estimates: 
comp_field_ld52x3D7_pval <- pchisq(2*(fit_both$minimum-fit_diff$minimum),df=1,lower.tail = FALSE)


# Save output -------------------------------------------------------------

cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(data.est.comp,glue::glue("output/Data-estimates-comparison_{Sys.Date()}_{cur_time}.xlsx"))
save(comp_field_ld52x3D7_pval,file=glue::glue("output/Data-estimates-comparison_field_LD52x3D7_{Sys.Date()}_{cur_time}.RData"))

# clean up ----------------------------------------------------------------

rm(j,i,nllh_both,fit_both,nllh_diff,fit_diff,p.val,cur_time)
