SNHT  <- function(serie,n_period=10,dstr='norm',simulations = 1000){
  if(exists(x = '.Random.seed')){
      old <- .Random.seed
      on.exit( { .Random.seed <<- old } )
  }
  serie <- as.vector(serie)
  n <- length(serie)
  na_ind <- is.na(serie)
  n_no_na <- n - sum(as.numeric(na_ind))
  if(n < 2*n_period){
    stop('serie no long enough or n_period too long')
  }
  i_ini <- n_period
  i_fin <- n-n_period

  serie_mean <- mean(serie,na.rm = T)
  serie_sd <- sd(serie,na.rm = T)

  if(dstr == 'gamma'){
    delta <- min(serie[!na_ind]) - 1
    par_gamma <- fitdistr(serie[!na_ind]-delta,'gamma')
  }
  serie_com <- serie[!na_ind]
  serie <- (serie - serie_mean)/serie_sd

  t <- rep(0,i_fin)
  for(i in i_ini:i_fin){
    z1 <- mean(serie[1:i],na.rm = T)
    z2 <- mean(serie[(i+1):n],na.rm = T)
    i_no_na <- i - sum(as.numeric(is.na(serie[1:i])))
    t[i] <- i_no_na*z1**2 + (n_no_na-i_no_na)*z2**2
  }
  i_break <- which.max(t)
  t_cri <- max(t)

  a_sim <- vector(mode = 'double',length = simulations)

  #Begin simulations:
  set.seed(14243)
  if(dstr == 'norm'){
    for(j in 1:simulations){
      aux <- rnorm(n_no_na,mean=serie_mean,sd = serie_sd)
      sd_aux <- sd(aux)
      mn_aux <- mean(aux)
      aux <- (aux - mn_aux)/sd_aux
      t <- rep(0,n_no_na)
      for(i in n_period:(n_no_na-n_period-1)){
        z1 <- mean(aux[1:i])
        z2 <- mean(aux[(i+1):n_no_na])
        t[i] <- i*z1**2 + (n_no_na-i)*z2**2
      }
      a_sim[j]<- max(t)
    }
  } else if( dstr == 'gamma'){
    for(j in 1:simulations){
      aux <- rgamma(n=n_no_na,shape=par_gamma$estimate[1],rate = par_gamma$estimate[2])
      aux <- aux + delta
      sd_aux <- sd(aux)
      mn_aux <- mean(aux)
      aux <- (aux - mn_aux)/sd_aux
      t <- rep(0,n_no_na)
      for(i in n_period:(n_no_na-n_period-1)){
        z1 <- mean(aux[1:i])
        z2 <- mean(aux[(i+1):n_no_na])
        t[i] <- i*z1**2 + (n_no_na-i)*z2**2
      }
      a_sim[j]<- max(t)
    }

  } else if (dstr == 'self'){
    for(j in 1:simulations){
      aux <- sample(x = serie_com,replace = T,size = n_no_na)
      sd_aux <- sd(aux)
      if(sd_aux == 0){
        next
      }
      mn_aux <- mean(aux)
      aux <- (aux - mn_aux)/sd_aux
      t <- rep(0,n_no_na)
      for(i in n_period:(n_no_na-n_period-1)){
        z1 <- mean(aux[1:i])
        z2 <- mean(aux[(i+1):n_no_na])
        t[i] <- i*z1**2 + (n_no_na-i)*z2**2
      }
      a_sim[j]<- max(t)
    }

  }else{
    stop('not supported dstr input')
  }
  #Check p.value
  cum_dist_func <- ecdf(a_sim)
  p <- 1-cum_dist_func(t_cri)
  out <- list(breaks = i_break+1 ,p.value = p)
  return(out)
}
