jsmm.mcmc <- function(loglikelihood = loglik, data = data, n.iter = n.iter, n.adapt.iter = iter, n.thin = 1, rotate = TRUE){
  
  # The following objects should be changed to use different prior distributions: mu_z, S_z, PSI, nu, rhos, priorRho

  tpm0 <- proc.time() #initial processing time
  
  ## First, we separate some variables from the "data" list
  ns <- data$ns # number of studied species
  np <- data$np # number of species-specific parameters
  SP <- data$SP # sequence of studied species
  includePhylogeny <- data$includePhylogeny
  includeTraits <- data$includeTraits
  CC <- data$CC # Phylogenetic correlation matrix
  TT <- data$TT # Species traits matrix
  
  ## And, then, we set up some variables and prior distributions.
  
  ## To include the intercept into the model, we set TT[,1] <- 1 for all species. In the absence of trait information, only the intercept is included.
  if(!includeTraits){
    TT <- matrix(1, nrow = ns)
    nt <- 1
    colnames(TT) <- "intercept"
  }

  if(includeTraits){
    TT <- TT
    nt <- ncol(TT)
    # traits <- colnames(TT)
  }
 
  ## Prior distributions:
  
  ## Zeta
  mu_z <- matrix(0, nrow = np * nt)
  S_z  <- 100 * diag(np * nt)

  ## Sigma
  PSI <- diag(np)
  nu  <- np

  ## Rho
  ## Pre-compute inverse and determinant of the D-matrix

  if(includePhylogeny){
  
    rhos <- seq(0, 1, 0.01)
    nr   <- length(rhos)
    # priorRho  <- rep(1/nr, nr)
    priorRho <- rep(0.5/(nr - 1), nr)
    priorRho[1] <- 0.5
    lpriorRho <- log(priorRho)
    iWs <- array(dim = c(ns, ns, nr))
    ldetWs <- numeric(nr)
  
    for(i in 1: nr){
      W <- rhos[i] * CC + (1 - rhos[i]) * diag(ns)
      iWs[ , , i] <- solve(W)
      ldetWs[i]   <- log(det(as.matrix(W)))
    }
  }
  
  ## Before running the MCMC we set initial values for the
  ## parameters to be estimated: Theta, Sigma, Zeta, and rho.
  
  if (is.null(data$THETAINIT)){
    Theta <- matrix(0,nrow = ns,ncol=np)
  } else {
    Theta <- data$THETAINIT # this ns x np matrix should be added to the data list in order to use initial values different from 0
  }
  
  Sigma  <- diag(np)  # identity matrix
  iSI <- solve(Sigma) # inverse of matrix Sigma
  Z   <- matrix(0, nrow = nt, ncol = np)
  M   <- TT %*% Z   # mu_k; expected values of THETA
  X   <- TT %x% diag(np)  
  
  if(includePhylogeny){
    RI <- round(nr/2)
    rho <- rhos[RI]
    iW <- iWs[ , ,RI]
    ldetW <- ldetWs[RI]
  }

  if(!includePhylogeny){
    rho <- 1
    iW <- diag(ns)
    # ldetW <- log(det(diag(ns)))
  }
  
  ## Initial likelihoods for THETA
  
  li1 <- loglikelihood(Theta, data)
 
  RES <- as.numeric(t(M)) - as.numeric(t(Theta))
  li2 <- -(1/2) * RES %*% kronecker(iW, iSI) %*% RES
  
  print("inital log-likelihood: ")
  print(c(li1,li2))
  print("sampling starts")
  
  ## Acceptance rates
  
  ac <- array(0, c(ns, np, 2))
  kk <- array(1, c(ns, np)) # sd for the proposal distributions
  acr <- array(NA, c(ns, np))
  
  ns1s2 <- array(0, c(ns))
  s1 <- array(0, c(ns, np))
  s2 <- array(0, c(ns, np, np))
  la <- array(1, c(ns, np))
  vect <- array(0, c(ns, np, np))
  
  for (k in 1:ns){
    vect[k,,] <- diag(np)
  }
  
  ## Posteriors to be stored
  Post_Theta <- array(NA, c(n.iter, ns, np))
  Post_Z <- array(NA, c(n.iter, nt, np))
  Post_Sigma <- array(NA, c(n.iter, np, np))
  Post_rho <- numeric(n.iter)
  Post_LIKE <- array(NA, c(n.iter, ns))
  
  ##  MCMC sampling scheme
  
  for (i in 1:(n.iter + n.adapt.iter)){
    
    print(i)
    
    for (ii in 1:n.thin){
      
      ## UPDATE THETA
      
      for (l in 1:np){
        
        NTHETA <- Theta ###
        
        for(k in 1:ns){
          
          nTHETA <- Theta[k, ]
          mult <- rnorm(1, mean=0, sd = (kk[k,l]*sqrt(la[k,l])))
         
          for (l2 in 1:np){
            
            nTHETA[l2] <- nTHETA[l2] + mult*vect[k,l,l2]
            
          }
          
          NTHETA[k, ] <- nTHETA
        }
        
        nli1 <- loglikelihood(NTHETA, data)
          
        for (k in 1:ns){
          
          N2THETA <- Theta 
          N2THETA[k,] <- NTHETA[k,]
          RES <- as.numeric(t(M)) - as.numeric(t(N2THETA))
          nli2 <- -(1/2) * RES %*% kronecker(iW, iSI) %*% RES
          ac[k,l,1] <- ac[k,l,1] + 1

          if(is.finite(nli1[k]) & is.finite(nli2)){
            
            if(runif(1) < exp(nli1[k] - li1[k] + nli2 - li2)){
              
              Theta[k,] <- NTHETA[k,]
              li1[k] <- nli1[k]
              li2 <- nli2
              ac[k,l,2] <- ac[k,l,2] + 1
              
            }
          }
        }
      }
      
      ## The other parameters are sampled directly from their full conditional distribution 
  
      ## UPDATE Z
  
      XTHE <- as.vector(t(Theta)) # transform to vector by line
      iXSI <- iW %x% iSI
      Vs <- solve(solve(S_z) + t(X) %*% iXSI %*% X)
      Vs <- (Vs + t(Vs))/2
      mus <- Vs %*% (solve(S_z) %*% mu_z + t(X) %*% iXSI %*% XTHE)
      # set.seed(42)
      Z <- matrix(rmvnorm(1, mean = mus, sigma = Vs), ncol = np, byrow = TRUE)
      M <- TT %*% Z
  
      ## UPDATE Sigma
  
      RES <- Theta - M
      A <- t(RES) %*% iW %*% RES
      PSIA <- PSI + A
      PSIA <- (PSIA + t(PSIA))/2
      ##set.seed(42)
      Sigma <- riwish((nu+ns), PSIA) # Inverse Wishart Matrix Distribution
      Sigma <- (Sigma + t(Sigma))/2
      iSI <- solve(Sigma)
  
      ## UPDATE rho
      
      if(includePhylogeny){
        RES <- as.numeric(t(M)) - as.numeric(t(Theta))
        likeRho <- numeric(nr)
        
        for(ii in 1:nr){
          likeRho[ii] <- (-1/2)*(np*ldetWs[ii] + RES%*%(iWs[,,ii] %x% iSI)%*%RES)
        }
        
        postRho <- lpriorRho + likeRho
        pr <- exp(postRho)/sum(exp(postRho))
        RI <- sample(seq(1:nr), size = 1, prob = pr)
        iW <- iWs[,,RI]
        ldetW <- ldetWs[RI]
        rho <- rhos[RI]
        
        RES <- as.numeric(t(M)) - as.numeric(t(Theta))
        li2 <- -(1/2) * RES %*% kronecker(iW, iSI) %*% RES
      }
      
      ## Adaptation
      
      if (i <= n.adapt.iter){
        
        for(k in 1:ns){
          
          q <- 1 + exp(-(i*n.thin)/500)
          w <- 1 - 0.1*exp(-(i*n.thin)/500)
          acr[k,] <- ac[k,,2]/ac[k,,1]
          kk[k,] <- sapply(kk[k,] * q^(acr[k,] - 0.44), trunca)
          s1[k,] <- s1[k,]+w*Theta[k,]
          s2[k,,] <- s2[k,,]+w*(Theta[k,]%*%t(Theta[k,]))
          ns1s2[k] <- ns1s2[k]+w
          
          if (rotate && ((i*n.thin)>50)){
            cov <- (s2[k,,]-(s1[k,]%*%t(s1[k,])/ns1s2[k]))/(ns1s2[k]-1)
            met <- cov + 10^(-5)*diag(np)
            lavect <- eigen(met)
            la[k,] <- abs(lavect$values)
            vect[k,,] <- lavect$vectors
          }
          
        }
        
        ac = ac*w
        
      }
      
    } # END OF THINNING LOOP
    
    ## Store the posteriors
    
    if(i > n.adapt.iter) {
      
      Post_Theta[i-n.adapt.iter,,] <- Theta
      Post_Z[i-n.adapt.iter,,] <- Z
      Post_Sigma[i-n.adapt.iter,,] <- Sigma
      Post_rho[i-n.adapt.iter] <- rho
      Post_LIKE[i-n.adapt.iter,] <- li1
      
    }
    
  }
  
  tpm1 <- proc.time() #final time
  
  ## Organizing outputs
  
  dtpm <- (tpm1[3] - tpm0[3])/60 # time elapsed
  
  rnames_summary <- c("Time elapsed (minutes):", "Total number of iterations:",
                      "Number of adaptation iterations:","Length of thinned posterior",
                      "Number of species", "Number of parameters", "Number of traits")
  results_summary <- matrix(c(as.numeric(dtpm), (n.iter+n.adapt.iter)*n.thin, n.adapt.iter*n.thin, n.iter, ns, np, nt),
                            dimnames = list(rnames_summary,"Model summary"))
  
  posterior <- list(THETA = Post_Theta, ZETA = Post_Z, SIGMA = Post_Sigma, RHO = Post_rho, LIKE = Post_LIKE)
  
  return(list(results_summary = results_summary, posterior = posterior, la = la, vect = vect, ac = ac, kk = kk))
  
}