## fits a DLM via univariate Gibbs sampling
##
## The model is:
##
## y_t = theta_t + epsilon_t,            epsilon_t ~ N(0, sigma2)
## theta_t = theta_{t-1} + omega_t,      omega_t ~ N(0, W_t)  
##
## Kevin Quinn
## 10/25/2017
##


library(MCMCpack) ## needed for rinvgamma()
set.seed(82382)

## load data
load("TrumpDailyApproval.Rda")
y <- TrumpDailyApproval$approval
y.obs <- y
T <- length(y)

## prior for mu0
theta0.mean <- 50
theta0.var  <- 100

## prior for sigma2
alpha <- 5
beta <- 5

## evolution variance for theta
Wt <- 0.1

## initialize parameter values
sigma2 <- 1
theta <- rep(mean(y, na.rm=TRUE), length(y))



## MCMC parameters
burnin <- 1000
gibbs  <- 5000
tot.iter <- burnin + gibbs

## storage
theta.storemat <- matrix(NA, gibbs, T)
sig.storemat <- matrix(NA, gibbs, 1)
count <- 0

### The Gibbs Sampler
for (iter in 1:tot.iter){

    ## impute missing y values
    y <- rnorm(T, m=theta, s=sqrt(sigma2))
    y[!is.na(y.obs)] <- y.obs[!is.na(y.obs)]
    
    
    ## sample theta
    for (i in 1:T){
        if (i==1){
            ## sample theta0
            theta.var <- 1/(1/theta0.var + 1/Wt)
            theta.mean <- theta.var*(theta0.mean/theta0.var + theta[1]/Wt)
            theta0 <- rnorm(1,theta.mean, sqrt(theta.var))
            
            ## sample theta[1]
            theta.var  <- 1/(2/Wt + 1/sigma2) 
            theta.mean <- theta.var*(theta0/Wt + y[i]/sigma2 + theta[i+1]/Wt)
            theta[i] <- rnorm(1, theta.mean, sqrt(theta.var))
        }
    else if (i==T){
        theta.var <- 1/(1/Wt + 1/sigma2)
        theta.mean <- theta.var*(y[i]/sigma2 + theta[i-1]/Wt)
        theta[i] <- rnorm(1,theta.mean, sqrt(theta.var))
    }
    else{
        theta.var  <- 1/(2/Wt + 1/sigma2) 
        theta.mean <- theta.var*(theta[i-1]/Wt + y[i]/sigma2 + theta[i+1]/Wt)
        theta[i] <- rnorm(1, theta.mean, sqrt(theta.var))
    }
    }
    
    ## sample sigma2
    e <- y - theta
    sse <- crossprod(e)
    sigma2 <- rinvgamma(1, (alpha+T)/2, (beta+sse)/2)
    
    ## store draws
    if (iter > burnin){
        count <- count + 1
        theta.storemat[count,] <- theta
        sig.storemat[count,] <- sigma2
    }
    
    ## print intermediate results
    if (iter%%500==0){
        cat("iterations = ", iter, "\n")
    }
}


## plot data and posterior quantities
thetahat <- apply(theta.storemat, 2, mean)
thetavar <- apply(theta.storemat, 2, var)
thetastd <- sqrt(thetavar)
plot(TrumpDailyApproval$day, y.obs,
     ylim=c(30,50),
     xlab="Time", ylab="Trump Approval",
     main="Quinnipiac University Trump Approval Polls: Voters")
xx <- c(TrumpDailyApproval$day, rev(TrumpDailyApproval$day))
yy <- c(apply(theta.storemat, 2, quantile, .025),
        rev(apply(theta.storemat, 2, quantile, .975)))
polygon(xx, yy, col="wheat", border="white")
lines(TrumpDailyApproval$day, thetahat, col="red", lwd=2)
points(TrumpDailyApproval$day, y.obs)




## turn output into a coda mcmc object
library(coda)
post.1 <- mcmc(cbind(sig.storemat, theta.storemat))
theta.names <- paste("theta", 1:T, sep="")
varnames(post.1) <- c("sigma2", theta.names)





