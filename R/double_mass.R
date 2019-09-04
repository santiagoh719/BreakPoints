
double_mass <- function(serie,ploting=T,date_axis=NULL,simulations = 10000, alpha = 0.5){
    if(exists(x = '.Random.seed')){
      old <- .Random.seed
    }
    serie <- as.vector(serie)
    n <- length(serie)
    if(!is.null(date_axis) & length(date_axis) != n){date_axis=NULL}
    na_ind <- is.na(serie)
    mn <- mean(serie,na.rm = T)
    desv <- sd(serie,na.rm = T)
    delta <- min(serie[!na_ind])-1
    par_gamma <- fitdistr(serie[!na_ind]-delta,'gamma')
    
    # Make cummulative sum of the serie
    serie[na_ind] <- 0
    serie <- cumsum(serie)
    
    # Begin simulations
      set.seed(14243)
    mat1 <- matrix(NA, nrow = n, ncol = simulations)
    mat2 <- matrix(NA, nrow = n, ncol = simulations)
    for(i in 1:simulations){
      mat1[,i] <- rnorm(n,mn,desv)
      mat1[na_ind,i] <- 0
      mat1[,i] <- cumsum(mat1[,i])
      mat2[,i] <- rgamma(n,shape=par_gamma$estimate[1],rate = par_gamma$estimate[2])
      mat2[,i] <- mat2[,i] + delta
      mat2[na_ind,i] <- 0
      mat2[,i] <- cumsum(mat2[,i])
    }
    # Compute the boundaryes
    q1_1 <- apply(mat1,1,quantile,probs=alpha/2)
    q3_1 <- apply(mat1,1,quantile,probs=1-alpha/2)
    q1_2 <- apply(mat2,1,quantile,probs=alpha/2)
    q3_2 <- apply(mat2,1,quantile,probs=1-alpha/2)
    
    if(exists(x = 'old')){
      old <- .Random.seed
      .Random.seed <- old
    }
    # Make ggplot objet if ploting = T
    if(ploting){
      if(is.null(date_axis)){
        dat <- data.frame(serie=serie,q11=q1_1,q31=q3_1,q12=q1_2,q32=q3_2,x=1:n)
        graf <- ggplot() + geom_line(data=dat, aes(y=serie,x=x),col='red')
        graf <- graf + geom_line(data=dat, aes(y=q11,x=x),linetype="dotted", col='blue')
        graf <- graf + geom_line(data=dat, aes(y=q31,x=x),linetype="dotted", col='blue')
        graf <- graf + geom_line(data=dat, aes(y=q12,x=x),linetype="dotted", col='green')
        graf <- graf + geom_line(data=dat, aes(y=q32,x=x),linetype="dotted", col='green')
      }else{
        dat <- data.frame(serie=serie,q11=q1_1,q31=q3_1,q12=q1_2,q32=q3_2,dates=date_axis)
        graf <- ggplot() + geom_line(data=dat, aes(y=serie,x=dates),col='red')
        graf <- graf + geom_line(data=dat, aes(y=q11,x=dates),linetype="dotted", col='blue')
        graf <- graf + geom_line(data=dat, aes(y=q31,x=dates),linetype="dotted", col='blue')
        graf <- graf + geom_line(data=dat, aes(y=q12,x=dates),linetype="dotted", col='green')
        graf <- graf + geom_line(data=dat, aes(y=q32,x=dates),linetype="dotted", col='green')
      }
    return(graf)
    }else{
      return(list(real=serie,lower_norm=q1_1,upper_norm=q3_1,lower_gamma=q1_2,upper_gamma=q3_2))
    }

}
