## model is: y_i ~ N(mu, Sigma)
##           2x1    2x1   2x2
##
##     with priors:
##         mu ~ N(m, M^-1)
##         Sigma ~ InvWishart(s, S)
##
## Inputs:
##         Y:  n x 2 matrix of observed data
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

MCMCBivariateNormal <- function(Y, m, M, s, S, burnin=500, ngibbs=2000){

    tot.iter <- burnin + ngibbs
    mu.store <- matrix(NA, ngibbs, 2)
    Sigma.store <- matrix(NA, ngibbs, 3) ## each row is vech(Sigma)

    ## pre-compute means of Y
    y.mean <- apply(Y, 2, mean)
    ## pre-compute number of obervations
    n <- nrow(Y)

    ## starting value for Sigma
    Sigma <- matrix(c(1,0,0,1), 2, 2, byrow=TRUE)
    Sigma.inv <- chol2inv(chol(Sigma))
    
    for (iter in 1:tot.iter){

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










