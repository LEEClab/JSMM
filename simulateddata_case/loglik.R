loglik <- function(THETA,data){
  
  all.tracks <- data$tracks
  SPS <- all.tracks[,1]
  corridor <- data$corridor
  patch <- data$patch
  di <- data$di
  
  SP  <- data$SP
  ns <- length(SP)
  
  like<-rep(0,ns)
  for (k in 1:ns) {
    theta = THETA[k, ]
    tracks = all.tracks[SPS == SP[k],]
    
    alpha <- exp(theta[1])
    betaC   <- theta[2]
    betaP   <- theta[3]
    li    <- 0
    
    idt <- tracks[,2]
    uid <- unique(idt)
    ta <- tracks[,3]
    
    for(ss in 1:length(uid)){ #for each one of the tracks
      lta <- ta[idt==uid[ss]] #data for the focal track
      
      for(i in 1:(length(lta)-1)){ #for each movement step
        x <- lta[i]
        pr <- exp(-di[x,]/alpha) * exp(betaC * corridor) * exp(betaP * patch)
        pr <- pr / sum(pr)
        pr <- pr + 10^(-10)
        pr <- pr / sum(pr)
        li <- li + log(pr[lta[i + 1]])
      }
    }
    like[k]<-li
  }
  return(like)
}
