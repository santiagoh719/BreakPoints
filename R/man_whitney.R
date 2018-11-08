man.whi <- function(serie,n_period=10){
  serie <- as.vector(serie)
  n <- length(serie)
  if(n < 2*n_period){
    stop('serie no long enough or n_period to long')
  }
  i_ini <- n_period
  i_fin <- n-n_period
  p_v <- 1
  for(i in i_ini:i_fin){
    aux1 <- serie[1:i]
    aux2 <- serie[(i+1):n]
    p <- wilcox.test(aux1,aux2, paired = F, var.equal = F)$p.value
    if(p < p_v){
      p_v <- p
      i_break <- i+1
    }
  }
  out <- list(breaks = i_break ,p.value =p_v)
  return(out)
}
