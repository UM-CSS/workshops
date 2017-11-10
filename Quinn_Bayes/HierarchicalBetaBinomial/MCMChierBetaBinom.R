## MCMC for Beta-Binomial Hierarchical model with exponential prior on
## alpha and beta
##
## The model is:
##
## [y_i | theta_i] ~ Binomial(m_i, theta_i)     1=1,...,n
## [\theta_i | alpha, beta] ~ Beta(alpha, beta)
## alpha ~ Exponential(a)
## beta ~ Exponential(b)
##
##
## Inputs:
##   m: n-vector of sample sizes for each of the n first-level units
##
##   y: n-vector of outcome counts for each of the n first-level units
##
##   a: rate parameter for the Exponential distribution prior for alpha
##
##   b: rate parameter for the Exponential distribution prior for beta
##
##   method: MCMC method. "marginal" samples (alpha, beta) from the
##             posterior distributions marginalized over theta and then
##             samples theta conditional on (alpha, beta). "conditional"
##             iteratively samples from [alpha, beta | theta] and
##             [theta | alpha, beta].
##
##   burnin: the number of initial "burn-in" iterations to be discarded
##
##   thin: the thinning interval. Every thin-th draw will be stored.
##
##   mcmc: the number of MCMC iterations after the burn-in period
##
##   tune: a positive tuning factor to adjust the Metropolis-Hastings
##         sampling. tune should be set so that the acceptance rate is
##         between 0.20 and 0.50.
##
##   verbose: how often should information be printed to the screen.
##            Defaults to 0 (never). Setting verbose=t would send output
##            to the screen every t-th iteration.
##
##   startvals: 2-vector of starting values for alpha and beta. Defaults to
##              NULL which sets the starting values equal to the prior mean
##              of alpha and beta.
##
##
##
## Output: a coda mcmc object containing the draws from the
##         posterior distribution of (alpha, beta, theta)
##
##
## Requires the mvtnorm and coda packages
##
## Kevin Quinn
## 2/10/2012
##

MCMChierBetaBinom <- function(m, y, a=.01, b=.01,
                              method=c("marginal", "conditional"),
                              burnin=500, thin=1, mcmc=10000, tune=1.5,
                              verbose=0, startvals=NULL){

  library(mvtnorm)
  library(coda)


  ## some utility functions
  dbetabinom <- function(x, n, a, b, log=FALSE){
    logfun <- lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b)
    if (log==TRUE){
      return(logfun)
    }
    else{
      return(exp(logfun))
    }
  }


  ## log p(alpha, beta | y)
  log.marg.posterior <- function(alphabeta, m, y, a, b){
    alpha <- alphabeta[1]
    beta <- alphabeta[2]
    piece1 <- dexp(alpha, a, log=TRUE) + dexp(beta, b, log=TRUE)
    piece2 <- sum(dbetabinom(y, m, alpha, beta, log=TRUE))
    logfunval <- piece1 + piece2
    return(logfunval)
  }
  
  ## log p(alpha, beta | y, theta) = log p(alpha, beta | theta)
  ## equality above follows from (alpha, beta) \indep y | theta
  log.cond.posterior <- function(alphabeta, a, b, theta){
    alpha <- alphabeta[1]
    beta <- alphabeta[2]
    piece1 <- dexp(alpha, a, log=TRUE) + dexp(beta, b, log=TRUE)
    piece2 <- sum(dbeta(theta, alpha, beta, log=TRUE))
    logfunval <- piece1 + piece2
    return(logfunval)
  }



  
  ## input checking
  if (length(m) != length(y)){
    stop("m is not the same length as y\n")
  }
  if (min(m) <= 0){
    stop("min(m) not positive\n")
  }
  if (min(y) < 0){
    stop("min(y) is negative\n")
  }
  if (a <= 0){
    stop("a is non-positive\n")
  }
  if (b <= 0){
    stop("b is non-positive\n")
  }
  if (is.null(startvals)){
    startvals <- c(1/a, 1/b)
  }
  if (length(startvals) != 2){
    stop("startvals not of length 2\n")
  }
  method <- match.arg(method, choices=c("marginal", "conditional"))

  
  ## defining some constants
  n <- length(y)
  tot.iter <- burnin + mcmc


  
  ## storage for the posterior draws
  alpha.samp <- rep(NA, mcmc/thin)
  beta.samp <- rep(NA, mcmc/thin)
  theta.samp <- matrix(NA, mcmc/thin, n)


  ## calculate the (marginal) posterior mode and Hessian for use in
  ## constructing proposal distributions for (alpha, beta)
  out <- optim(startvals, fn=log.marg.posterior, method="L-BFGS-B",
               lower=rep(1e-6, 2),
               control=list(fnscale=-1), hessian=TRUE,
               m=m, y=y, a=a, b=b)

  
  alphabeta.mean <- out$par
  alphabeta.var <- solve(-1*out$hessian) * tune^2




  ## The Marginal Sampling Method
  ##   sample from p(alpha, beta | y)  (Metropolis step)
  ##   then sample from
  ##   p(theta | y, alpha, beta)  (direct sampling from the conditional)
  ##
  accepts <- 0
  if (method == "marginal"){
    alphabeta.cur <- startvals
    counter <- 1
    for (iter in 1:tot.iter){

      
      ## sample from [alpha, beta | y]
      ## (random walk Metropolis step)
      alphabeta.cand <- rmvnorm(1, alphabeta.cur, alphabeta.var, method="chol")
      if (min(alphabeta.cand) > 0){
        ratio <- log.marg.posterior(alphabeta.cand, m, y, a, b) - log.marg.posterior(alphabeta.cur, m, y, a, b)
        
        if (runif(1) < exp(ratio) ){
          alphabeta.cur <- alphabeta.cand
          accepts <- accepts + 1
        }
      }
      
      
      if (iter > burnin & iter %% thin == 0){
        alpha.samp[counter] <- alphabeta.cur[1]
        beta.samp[counter] <- alphabeta.cur[2]
        ## sample from [theta | alpha, beta, y]
        for (i in 1:n){
          theta.samp[counter,i] <- rbeta(1, alphabeta.cur[1]+y[i],
                                         alphabeta.cur[2]+m[i]-y[i])
        }
        counter <- counter + 1
      }
      
      if ( verbose > 0 & iter %% verbose == 0){
        cat("iteration", iter, "of", tot.iter, "\n")
        cat("acceptance rate:", accepts/iter, "\n\n")
      }
            
    }
  } ## end marginal sampling method

  

  ## The Conditional Sampling Method
  ##    (Metropolis within Gibbs)
  ##   sample from p(alpha, beta | y, theta) (Metropolis step)
  ##   then sample from
  ##   p(theta | y, alpha, beta) (direct sampling from the conditional)
  ##
  accepts <- 0
  if (method == "conditional"){
    theta <- (y+1)/(m+2)  ## starting values for theta
    alphabeta.cur <- startvals
    counter <- 1
    for (iter in 1:tot.iter){


      ## sample from [alpha, beta | theta, y]
      ## (random walk Metropolis step)
      alphabeta.cand <- rmvnorm(1, alphabeta.cur, alphabeta.var, method="chol")
      if (min(alphabeta.cand) > 0){
        ratio <- log.cond.posterior(alphabeta.cand, a, b, theta) - log.cond.posterior(alphabeta.cur, a, b, theta)

        
        if (runif(1) < exp(ratio) ){
          alphabeta.cur <- alphabeta.cand
          accepts <- accepts + 1
        }
      }

      ## sample from [theta | alpha, beta, y]
      for (i in 1:n){
        theta[i] <- rbeta(1, alphabeta.cur[1]+y[i],
                          alphabeta.cur[2]+m[i]-y[i])
      }
      
      
      if (iter > burnin & iter %% thin == 0){
        alpha.samp[counter] <- alphabeta.cur[1]
        beta.samp[counter] <- alphabeta.cur[2]
        for (i in 1:n){
          theta.samp[counter,] <- theta
        }
        counter <- counter + 1
      }
      
      if ( verbose > 0 & iter %% verbose == 0){
        cat("iteration", iter, "of", tot.iter, "\n")
        cat("acceptance rate:", accepts/iter, "\n\n")
      }            
    }
  } ## end conditional sampling method



  outmat <- cbind(alpha.samp, beta.samp, theta.samp)
  colnames(outmat) <- c("alpha", "beta", paste("theta", 1:n, sep=""))
  outmat <- mcmc(outmat, start=burnin+1, end=tot.iter, thin=thin)
  
  return(outmat)
  
}


