###################################################################################
##      Joint species movement modelling: how do traits influence movements?     ##
##                     Fitting the JSMM to bird movement data                    ##
###################################################################################

## authors: "Otso Ovaskainen, Danielle Leal Ramos, Eleanor M. Slade, Thomas Merckx, Gleb Tikhonov, Juho Pennanen, Marco Aurélio Pizo, Milton Cezar Ribeiro, and Juan Manuel Morales"
## date: "November, 2018"

## This code can also be viewed in http://htmlpreview.github.io/?https://github.com/LEEClab/JSMM/blob/master/JSMM.html , with more detailed comments.

## Required packages
if (!require("MCMCpack")) install.packages("MCMCpack") # for function riwish
if (!require("mvtnorm")) install.packages("mvtnorm") # for functions dmvnorm and rmvnorm

## Landscape data ##

habitat <- read.table("type_pa.csv", sep = ",")
s <- numeric(length(habitat$V1))
f <- numeric(length(habitat$V1))
s[which(habitat$V1 == 2)] <- 1
f[which(habitat$V1 == 3)] <- 1

xp <- rep(c(1:60), 60)        # x coordinates
yp <- rep(c(60:1), each = 60) # y coordinates
xv <- cbind(xp, yp)
d  <- unname(as.matrix(dist(xv, upper = TRUE, diag = TRUE)) )

## Bird movement data ##

tracks <- read.table("tracks.csv", sep = ",", header = TRUE)

## Phylogenetic correlation matrix (optional)
CC <- read.table("CC.csv", sep = ",", header = FALSE)
CC <- apply(CC, FUN = as.numeric, MARGIN = 1)

## Species traits (optional)
ttraits <- read.table("traits.csv", sep = ",", header = TRUE)

TT <- matrix(0, nrow = nrow(ttraits), ncol = length(unique(ttraits[ , 2])) + 1) # contains the traits (columns) for each species (rows). Cathegorical traits are converted to binomial variables.

colpos = as.numeric(ttraits[ , 2])

for(i in 1:nrow(TT)){
  TT[i, colpos[i]] = 1 
}

TT[ , ncol(TT)] <- log(ttraits[ , 3]) # Log body mass

colnames(TT)<-c("frugivorous","insectivorous","omnivorous","granivorous","log(mass)")

## Note: the sequence of species in TT and CC matrix need to be the same, and the parameters of all species described in TT and CC will be estimated.

## Full data ##

data <- list(
  
  SP = ttraits[ ,1], # vector with species names
  ns = length(ttraits[ ,1]), # number of species
  CC = CC,
  TT = TT,
  tracks = tracks,
  
  semiopen = s,
  forest = f,
  d = d,
  
  includeTraits = TRUE, # change to FALSE to not include species traits in the model
  includePhylogeny = TRUE, # change to FALSE to not include species phylogeny in the model
  
  np = 3, # number of species-specific parameters to be estimated
  parameternames = c("log distance", "semiopen_aff", "forest_aff") # names of species-specific parameters to be estimated
  
)

## Likelihood function ##
## Likelihood function with three arguments: vector theta with species-specific parameters, species movement data (tracks), and environmental covariates (s, f, d).

loglik <- function(THETA,data){
  
  # bird movement data
  all.tracks <- data$tracks
  SPS <- all.tracks[,1]
  SP  <- data$SP
  ns <- length(SP)
  
  #environmental covariates
  semiopen <- data$semiopen
  forest <- data$forest
  d <- data$d
  
  like<-rep(0, ns)
  
  for (k in 1: ns) {
    theta <- THETA[k, ]
    tracks <- all.tracks[SPS == SP[k],]
    
    alpha <- exp(theta[1])
    beta_s   <- theta[2]
    beta_f   <- theta[3]
    li    <- 0
    
    idt <- tracks[ ,2]
    uid <- unique(idt)
    ta <- tracks[ ,3]
    
    for(ss in 1: length(uid)){ #for each one of the tracks
      lta <- ta[idt == uid[ss]] #data for the focal track
      
      for(i in 1: (length(lta)-1)) { #for each movement step
        x <- lta[i]
        pr <- exp(-d[x,]/alpha) * exp(beta_s * semiopen) * exp(beta_f * forest)
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

## JSMM function ##

trunca <- function(x) min(max(x, 10^-5), 10^5) # required function to run the JSMM function below

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

## Fitting the model ##

output <- jsmm.mcmc(loglikelihood = loglik, data = data, n.iter = 100, n.adapt.iter = 10, n.thin = 5, rotate = TRUE)

## Results ##

post <- output$posterior

np <- data$np
nt <- dim(post$ZETA)[2] # number of traits
ns <- data$ns

TT <- data$TT

mass <- TT[,5] # body mass for each bird species
tgroup = 1+(TT[,2]==1)+2*(TT[,3]==1)+3*(TT[,4]==1) # numerical vector with birds diet groups

meZT <- apply(post$ZETA, FUN = mean, MARGIN=2:3) # means for the posterior values of ZETA
M <- TT %*% meZT

tiff(filename="fig_1_birds.tif", res=300, unit="cm", height=17, width=8)
par(mfrow = c(np,1), oma=c(0,0,0,0),mar=c(5,5,3,1))
par(mfrow = c(np,1), oma=c(0,0,0,0),mar=c(2,5,1,1))

for (k in 1:np){  
  me <- apply(post$THETA[,,k], FUN = mean, MARGIN=2)
  lo <- apply(post$THETA[,,k], MARGIN = 2,FUN = quantile, probs = 0.25)
  hi <- apply(post$THETA[,,k], MARGIN = 2,FUN = quantile, probs = 0.75)
  
  plot(1, type="n", xlab="", ylab="",
       xlim=c(2,6),
       ylim=c(min(lo), max(hi)),
       cex=2.5,xaxt="n",yaxt="n",
       bty="n")
  
  box(bty="o", col="gray35")
  axis(1,cex.axis=1.6, col = "gray35",col.axis="gray35", las = 1,tck=-0.01,mgp = c(2, 0.5, 0),
       lwd = 1,at = seq(-2,6,1))
  axis(2, cex.axis = 1.6, col = "gray35", col.axis="gray35", las = 1, tck=-0.01, mgp = c(3, 0.4, 0), lwd = 1)
  #  title(xlab =if(case==2){"log mass"}else{"log wing span"}, line = 2.3, cex.lab=2)
  title(xlab ="", line = 2.3, cex.lab=2)
  # title(ylab = data$parameternames[k], line = 2.5,cex.lab=2)
  title(ylab = "", line = 2.5,cex.lab=2)
  
  xvar = mass
  
  for(i in 1: length(unique(tgroup))) {
    
    pdit = unique(tgroup)[i]
    
    T1 <- matrix(0, ncol = nt, nrow = 2)
    
    T1[1, nt] <- min(xvar)
    T1[2, nt] <- max(xvar)
    T1[1, pdit] <- 1
    T1[2, pdit] <- 1
    
    pre <- T1 %*% meZT
    
    lines(pre[,k] ~ c(min(xvar), max(xvar)), col = as.numeric(pdit),lwd = 1)
    points(me[tgroup == pdit] ~ xvar[tgroup == pdit], pch = 16, col = as.numeric(pdit))
    segments(xvar[tgroup == pdit], lo[tgroup == pdit],xvar[tgroup == pdit], hi[tgroup == pdit], col = as.numeric(pdit))
    
    # colors birds: (1,frugivorous,black), (2,insectivorous,red), (3,omnivorous, green), (4,granivorous, blue)
    
  }
}

dev.off()

post <- output$posterior

ndt = 4
traitnames = colnames(TT)[1:4]
TC <- matrix(0,ncol=nt, nrow=ndt) #traits 1-4
TC[1, 1] <- 1
TC[2, 2] <- 1
TC[3, 3] <- 1
TC[4, 4] <- 1

npost = dim(post$ZETA)[1]
pdif = array(0,c(np,ndt,ndt))
for (i in 1:npost){
  tmp = TC %*% post$ZETA[i,,]
  for (j in 1:np){
    for (k1 in 1:ndt){
      for (k2 in 1:ndt){
        pdif[j,k1,k2] = pdif[j,k1,k2]+(tmp[k1,j]>tmp[k2,j])/npost
      }
    }
  }
}

for (j in 1:np){
  print(data$parameternames[j])
  for (k1 in 1:(ndt-1)){
    for (k2 in (k1+1):ndt){
      if (min(1-pdif[j,k1,k2],pdif[j,k1,k2])<0.15){
        print(c(traitnames[k1],">",traitnames[k2]))
        print(pdif[j,k1,k2])
      }
    }
  }
}

## Probabilities of effects of traits and phylogeny > 0

print("effect of xvar (mass of birds)")
for (i in 1:np){
  print(data$parameternames[i])
  print(mean(post$ZETA[,nt,i]>0))
}

print("effect of diet")
for (i in 1:np){
  print(data$parameternames[i])
  print(mean(post$ZETA[,4,i]>0))
}

## Phylogeny
quantile(post$RHO,c(0.025,0.5,0.975))
mean(post$RHO>0)