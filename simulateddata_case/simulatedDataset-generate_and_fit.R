###################################################################################
##      Joint species movement modelling: how do traits influence movements?     ##
## Testing the validity of the posterior sampling scheme with simulated datasets ##
###################################################################################

## authors: "Otso Ovaskainen, Danielle Leal Ramos, Eleanor M. Slade, Thomas Merckx, Gleb Tikhonov, Juho Pennanen, Marco Aurélio Pizo, Milton Cezar Ribeiro, and Juan Manuel Morales"
## date: "November, 2018"

## 1. Before running the code below, please set the working directory to a folder containing all the files downloaded from the "simulateddata_case" folder.

## 2. Required packages & functions

library(mvtnorm) # for function dmvnorm and rmvnorm
library(MCMCpack) # for function riwish
library(Matrix)
library(doParallel) # for parallel computing

source("trunca.r")
source("jsmm_mcmc.r")
source("loglik.r")

## 3. Creating datasets with know parameter values

for (repl in 1:100){ # to create replicates
  n <- 100 # number of species
  m <- 110 # number of steps; we will delete first 10 positions later
  
  for (nl in c(FALSE, TRUE)){ # nl is "non-linear"; it's set to TRUE when the movement distance is assumed to behave non-linearly with body size

    lbody <- runif(n, min = 2, max = 6) # sampling random values of log body sizes from uniform distribution within the range [2,6] 
    lbody2 <- lbody^2
    frugivore <- rbinom(n, size=1, prob=0.5) # each species with equal probability to be frugivorous (1) or insectivorous (0)
   
    if (nl){
      zeta <- matrix(c(-7, 1, 4, -0.5, 0, 0, 0.5, 0, 0, 2, 0, 0), nrow = 4, ncol = 3) # true parameter if movement distance is assumed to behave non-linearly with body size
    } else {
      zeta <- matrix(c(-2, 1, 0.5, 0, 0, 0, 0.5, 0, 0, 2, 0, 0), nrow = 4, ncol = 3) # true parameter if movement distance is assumed to behave linearly with body size
    }

    T <- cbind(1, frugivore, lbody, lbody2) # n×3 trait matrix;columns correspond to intercept (all values set to one), an indicator variable describing if the species is a frugivore (with value one) or not (with value zero), log-body size, and squared log-body size
    mu <- T %*% zeta
    var <- c("log distance", "affinity to semi-open areas", "forest affininity")
    PC <- diag(n)

    for (i in seq(2, n, 2)) {
      PC[i-1, i] <- 0.9
    }

    for (i in seq(2, n, 2)) {
      PC[i, i-1] <- 0.9
    } 
    
    Sigma <- matrix(c(1, 0, 0, 0, 0.5, 0, 0, 0, 0.2), nrow = 3, ncol = 3)
    rho <- 0.7 # strength of the phylogenetic signal
    V <- kronecker(rho * PC + (1-rho) * diag(n), Sigma)
    theta <- matrix(rmvnorm(1, mean = as.vector(t(mu)), sigma = V), ncol = 3, byrow = TRUE)
    
    ## Environmental attributes
    habitat <- read.table("type_pa.csv", sep = ",")
    s <- numeric(length(habitat$V1)) # semi-open
    f <- numeric(length(habitat$V1)) # forest
    s[which(habitat$V1 == 2)] <- 1
    f[which(habitat$V1 == 3)] <- 1

    ## Matrix of distances
    xp <- rep(c(1:60), 60)        # x coordinates
    yp <- rep(c(60:1), each = 60) # y coordinates
    xv <- cbind(xp, yp)
    d  <- unname(as.matrix(dist(xv, upper = TRUE, diag = TRUE)) ) # matrix of distances between pixels
    
    ## Creating maps
    ## 1. Creating horizontal X and vertical Y coordinates; each pixel = 10m
    Xspl <- xp*10
    Yspl <- rep(c(60:1), each=60)
    Yspl <- Yspl*10
    PID <- c(1:3600)
    
    ## 2. Creating a data frame with Xspl and Yspl tracks coordinates and values: pixel id, and habitat
    
    DTSET <- matrix(nrow = m, ncol = n) # matrix with tracks for each specie (column)
    DTSET[1, ] <- sample(PID, n, replace = TRUE)
    
    for (i in 1 : n) {
      ltas <- DTSET[,i]
      alpha <- exp(theta[i,1])
      beta_s <- theta[i,2]
      beta_f <- theta[i,3]
      
      for (ii in 2 : m) {
        x <- ltas[ii - 1]
        pm <- exp(-d[x, ]/alpha) * exp(beta_s * s) * exp(beta_f * f)
        pm <- pm/sum(pm)
        ltas[ii] <- sample(PID, size = 1, prob = pm)
      }

      DTSET[,i] <- ltas
      
    }
    
    CCALL <- PC
    TTALL <- T
    
    SPALL <- 1 : n
    parameternames <- c("log distance", "semiopen_aff", "forest_aff")
    np <- length(parameternames)
        
    for (mm in  c(10,100)) {

      for (nn in c(10,100)) {
        
        spsample <- sample(c(1:100), 100)[1:nn]
        SP <- SPALL[spsample]
        ns <- nn
        TT <- TTALL[spsample, ]
        CC <- CCALL[spsample, spsample]
        
        tracks <- DTSET
        tracks <- tracks[((nrow(tracks) - mm) + 1): nrow(tracks), spsample]
        TA <- as.vector(unlist(tracks))
        SPS <- rep(SP, each = nrow(tracks))
        TID <- rep(c(1:nn), each = nrow(tracks))
        
        tracks <- cbind(SPS, TID, TA)
        
        trueValues <- list(zeta = zeta, rho=rho, Sigma = Sigma, theta = theta[spsample,])

        data <- list(SP = SP, ns = ns,  np = np, semiopen = s, forest = f, d = d, 
                     CC = CC, TT = TT, tracks = tracks, includeTraits = TRUE, includePhylogeny = TRUE,
                     parameternames = parameternames, trueValues = trueValues)
        save(data, file = paste("bird_data", if(nl){"_NL"}else{"_L"},
                                "_", toString(repl),
                                "n_", toString(nn),
                                "m_", toString(mm),
                                ".Rdata",sep=""))
      }
    }
  }
}

## 4. Fitting the JSMM to the dataset. This step is similar to the fitting of real bird data: http://htmlpreview.github.io/?https://github.com/LEEClab/JSMM/blob/master/JSMM.html .

cl <- makeCluster(15)
registerDoParallel(cl)
   
tpm <- foreach(repl = 1:15, .export = c('loglik'), .packages = c('Matrix','mvtnorm','MCMCpack')) %dopar%
   {
      for (mm in c(10,100)) {
         for (nn in c(10,100)) { 

            ## REPLICATE MODELS (WITH DATA GENERATED WITH LINEAR AND NON-LINEAR MODELS) FITTED TO LINEAR MODEL

            load(file = paste("bird_data_L_", toString(repl), "n_", toString(nn), "m_", toString(mm), ".Rdata", sep = ""))
            data$TT <- data$TT[ , c(1, 2, 3)]
            
            out <- jsmm.mcmc(loglikelihood = loglik, data = data, n.iter = 100, n.adapt.iter = 10, n.thin = 5, rotate = TRUE)
            outfile <- paste("bird_posterior_1000_200_1_", toString(repl), "n_", toString(nn), "m_", toString(mm), ".Rdata", sep="")
            save(out, file = outfile)
            
            ## REPLICATE MODELS (WITH DATA GENERATED WITH LINEAR AND NON-LINEAR MODELS) FITTED TO NON-LINEAR MODEL

            for (nl in FALSE){#c(FALSE,TRUE)){
               load(file = paste("bird_data_", if(nl){"NL_"}else{"L_"}, toString(repl), "n_", toString(nn), "m_", toString(mm), ".Rdata", sep = ""))
               out <- jsmm.mcmc(loglikelihood = loglik, data = data, n.iter = 100, n.adapt.iter = 10, n.thin = 5, rotate = TRUE)
               outfile <- paste("bird_posterior_1000_200_1_",if(nl){"NL_"} else{"L_"}, toString(repl), "n_", toString(nn), "m_", toString(mm), ".Rdata", sep="")
               save(out, file = outfile)
            }
         }
      }
   }
stopCluster(cl)
