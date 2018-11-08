pettit <- function(serie,n_period){
  serie <- as.numeric(as.vector(serie))
  n <- length(serie)
  if(n < 2*n_period){
    stop('serie no long enough or n_period to long')
  }
  i_ini <- n_period
  i_fin <- n-n_period
  U_v <- 0
  n_row <- n-i_ini
  aux1 <- matrix(serie[1:i_fin],ncol = i_fin,nrow = n_row)
  aux2 <- matrix(serie[(i_ini+1):n],ncol = i_fin,nrow = n_row,byrow = T)
  data <- sign(aux1-aux2)

  for(i in 1:(i_fin-i_ini+1)){
    aa <- i_ini -1 +i
    U <- abs(sum(data[1:aa,i:n_row],na.rm = T))
    if(U > U_v){
      U_v <- U
      i_break <- i + i_ini
    }
  }
  out <- list(breaks = i_break ,p.value = 2 * exp(-6*U_v**2/(n**3 + n**2)))
  return(out)
}
