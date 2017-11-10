## model is: y_i ~ N(mu, Sigma)
##           2x1    2x1   2x2
##
##     with priors:
##         mu ~ N(m, M^-1)
##         Sigma ~ InvWishart(s, S)
##
## Inputs:
##         Y:  n x 2 matrix of partially observed data. Missing values are
##             given by an NA.
##
##         m:  prior mean for mu (2-vector)
##
##         M:  prior precision for mu (2 x 2 matrix)
##   
##         s: prior degrees of freedom for Sigma (scalar)
##
##         S: prior scale matrix for Sigma (2 x 2 matrix)
##
##         burnin: number of burn-in iterations for Gibbs sampler
##
##         ngibbs: number of iterations after burnin
##
##
## NOTE: there is no error checking in the code below. The user is responsible
##       for ensuring that valid arguments are passed to the function below.
##       Use at your own risk.
##
## Kevin Quinn
## 10/29/2017
##

library(MCMCpack) ## needed for InvWishart() (also loads coda library)
library(mvtnorm)

MCMCMissingBivariateNormal <- function(Y, m, M, s, S, burnin=500, ngibbs=2000){

    ## subset Y by missingness pattern
    no.miss.indic <- !is.na(Y[,1]) & !is.na(Y[,2])
    n.no.miss <- sum(no.miss.indic)
    y1.miss.indic <- is.na(Y[,1]) & !is.na(Y[,2])
    n.y1.miss <- sum(y1.miss.indic)
    y2.miss.indic <- !is.na(Y[,1]) & is.na(Y[,2])
    n.y2.miss <- sum(y2.miss.indic)

    Y.no.miss <- Y[no.miss.indic,]
    Y.y1.miss <- Y[y1.miss.indic,]
    Y.y2.miss <- Y[y2.miss.indic,]
    ## no information when all data missing so exclude

    ## pre-compute number of obervations
    n <- n.no.miss + n.y1.miss + n.y2.miss

    
    tot.iter <- burnin + ngibbs
    mu.store <- matrix(NA, ngibbs, 2)
    Sigma.store <- matrix(NA, ngibbs, 3) ## each row is vech(Sigma)


    ## starting value for Sigma
    Sigma <- matrix(c(1,0,0,1), 2, 2, byrow=TRUE)
    Sigma.inv <- chol2inv(chol(Sigma))

    for (iter in 1:tot.iter){

        ## sample Y^(mis) | mu, Sigma, Y^(obs)
        ## first y1 | mu, Sigma, y2
        mean.fcd <- Sigma[1,2]/Sigma[2,2]*Y.y1.miss[,2]
        var.fcd <- Sigma[1,1] - (Sigma[1,2]^2)/Sigma[2,2]
        #cat("y1: var.fcd", var.fcd, "\n")
        Y.y1.miss[,1] <- rnorm(n.y1.miss, m=mean.fcd, s=sqrt(var.fcd))
            
        ## now y2 | mu, Sigma, y1
        mean.fcd <- Sigma[1,2]/Sigma[1,1]*Y.y2.miss[,1]
        var.fcd <- Sigma[2,2] - (Sigma[1,2]^2)/Sigma[1,1]
        #cat("y2: var.fcd", var.fcd, "\n")
        Y.y2.miss[,2] <- rnorm(n.y2.miss, m=mean.fcd, s=sqrt(var.fcd))

        Y <- rbind(Y.no.miss, Y.y1.miss, Y.y2.miss)
        y.mean <- apply(Y, 2, mean)
        
        ## sample mu | Sigma, Y
        cov.fcd <- chol2inv(chol(M + n*Sigma.inv))
        mean.fcd <- cov.fcd %*% (M %*% m + n*Sigma.inv %*% y.mean )
        mu <- rmvnorm(1, mean=mean.fcd, sigma=cov.fcd)
        
        
        ## sample Sigma | mu, Y
        Y.minus.mu <- Y
        Y.minus.mu[,1] <- Y.minus.mu[,1] - mu[1]
        Y.minus.mu[,2] <- Y.minus.mu[,2] - mu[2]
        SSE <- crossprod(Y.minus.mu)
        Sigma <- riwish(s+n, S+SSE)
        Sigma.inv <- chol2inv(chol(Sigma))
       

        if (iter %% 1000 == 0){
            cat("\niteration", iter, "of", tot.iter, "\n")
            print(mu)
            print(Sigma)
        }
        

        ## store draws if after burn-in period
        if (iter > burnin){
            mu.store[iter-burnin,] <- mu
            Sigma.store[iter-burnin,] <- vech(Sigma)
        }
        
    }


    ## put draws into coda object and return
    output <- cbind(mu.store, Sigma.store)
    output <- mcmc(data=output, start=burnin+1, end=tot.iter)
    colnames(output) <- c("mu1", "mu2", "Sigma11", "Sigma12", "Sigma22")
    return(output)
}


