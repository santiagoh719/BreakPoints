
sum_yam <- function(x){
  
  y <- x[(length(x)/2+1):length(x)]
  x <- x[1:(length(x)/2)]
  
  xmin <- mean(x)
  ymin <- mean(y)
  xsd <- sd(x)
  ysd <- sd(y)
  return(abs(xmin-ymin)/(xsd+ysd))
}


yamamoto <- function(serie, alpha = 0.1, n_period = 10){
  
  
  coef <- qt(p = 1-alpha,df = n_period-1) / sqrt(n_period)
  
  quiqui <- c(rep(NA,10),
              rollapply(serie,width=n_period*2,by=1,FUN=sum_yam) / coef,
              rep(NA,9))
  
  if(all(quiqui<=1,na.rm = T)){
    return(list(breaks=NULL, n=NULL))
  }
  
  quiqui1 <- quiqui > 1
  
  index_qui <- 1:length(quiqui)
  index_qui <- index_qui[quiqui1]
  index_qui <- index_qui[!is.na(index_qui)]
  # index_qui <- c(5,10,17:19,40:42,index_qui,62,65)
  if(length(index_qui) == 1){
    return(list(breaks=index_qui, n=1))
  }
  
  aux <- index_qui[2:length(index_qui)] - index_qui[1:(length(index_qui)-1)]
  
  while(any(aux==1)){
    
    id <- which(aux==1)
    
    vect <- id
    
    cont <- 0
    for( iii in id){
      cont <- cont + 1
      if(quiqui[index_qui[iii]] > quiqui[index_qui[iii+1]]){
        vect[cont] <- vect[cont]+1
      }
    }
    
    index_qui <- index_qui[-vect]
    
    if(length(index_qui) == 1){
      return(list(breaks=index_qui, n=1))
    }
    aux <- index_qui[2:length(index_qui)] - index_qui[1:(length(index_qui)-1)]
  }
  
  while(any(aux < n_period)){
    
    id <- which(aux < n_period)
    
    vect <- id
    
    cont <- 0
    for( iii in id){
      cont <- cont + 1
      if(quiqui[index_qui[iii]] > quiqui[index_qui[iii+1]]){
        vect[cont] <- vect[cont]+1
      }
    }
    
    index_qui <- index_qui[-vect]
    
    if(length(index_qui) == 1){
      return(list(breaks=index_qui, n=1))
    }
    aux <- index_qui[2:length(index_qui)] - index_qui[1:(length(index_qui)-1)]
  }
  
  return(list(breaks=index_qui[order(quiqui[index_qui]) < 6], n=length(index_qui)))
}
