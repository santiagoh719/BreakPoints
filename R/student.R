stu <- function(serie,n_period=10){
  if(exists(x = '.Random.seed')){
      old <- .Random.seed
      on.exit( { .Random.seed <<- old } )
  }
  serie <- as.vector(serie)
  n <- length(serie)
  if(n < 2*n_period){
    stop('serie no long enough or n_period too long')
  }
  i_ini <- n_period
  i_fin <- n-n_period
  p_v <- 1
  for(i in i_ini:i_fin){
    aux1 <- serie[1:i]
    aux2 <- serie[(i+1):n]
    p <- t.test(aux1,aux2, paired = F, var.equal = F)$p.value
    if(p <= p_v){
      p_v <- p
      i_break <- i+1
    }
  }
  out <- list(breaks = i_break ,p.value =p_v)
  return(out)
}
