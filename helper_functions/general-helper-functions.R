
# geometric mean for rows of a data.frame ---------------------------------

row_geom_mean <- function(data){
  apply(data, 1, function(x) ifelse(sum(x==0)==ncol(data),0,exp(mean(log(x)))))
}
