## fits a DLM via Gibbs sampling (efficient approach of Carter and Kohn 1994)
## notation as in West and Harrison 1997
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


## known parameter values
Ft <- 1 
Gt <- 1
Vt <- 1   ## observation equation error variance (sigma2)
Wt <- 0.1  ## evolution equation error variance 

## initial values
m0 <- 50
C0 <- 100
sigma2 <- 1
theta <- rep(mean(y, na.rm=TRUE), length(y))
a <- matrix(0,T,1)
R <- a
f <- a
Q <- a
e <- a
A <- a
m <- a
C <- a

## prior for sigma2
alpha <- 5
beta <- 5

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
    ## run the Kalman filter forward through time
    a[1] <- m0
    R[1] <- C0 + Wt
    f[1] <- a[1]
    Q[1] <- R[1] + Vt
    e[1] <- y[1] - f[1]
    A[1] <- R[1]/Q[1]
    m[1] <- a[1] + A[1]*e[1]
    C[1] <- R[1] - A[1] * Q[1] * A[1]
    for (i in 2:T){
        a[i] <- m[i-1]
        R[i] <- C[i-1] + Wt
        f[i] <- a[i]
        Q[i] <- R[i] + Vt
        e[i] <- y[i] - f[i]
        A[i] <- R[i]/Q[i]
        m[i] <- a[i] + A[i]*e[i]
        C[i] <- R[i] - A[i] * Q[i] * A[i]
    }
    
    ## sample theta[T]
    theta[T] <- rnorm(1, m[T], sqrt(C[T]))
    
    ## sample theta backwards through time
    for (i in (T-1):1){
        B <- C[i]/R[i+1]
        h <- m[i] + B * (theta[i+1] - a[i+1])
        H <- C[i] - B * R[i+1] * B
        theta[i] <- rnorm(1, h, sqrt(H))
    }
    
    ## sample sigma2
    e <- y - theta
    sse <- crossprod(e)
    sigma2 <- rinvgamma(1, (alpha+T)/2, (beta+sse)/2)
    Vt <- sigma2
    
    
    ## store draws
    if (iter > burnin){
        count               <- count + 1;
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
post.eff <- mcmc(cbind(sig.storemat, theta.storemat))
theta.names <- paste("theta", 1:T, sep="")
varnames(post.eff) <- c("sigma2", theta.names)

## acf plots (assumes dlm.R has already been run and output is post.1)
par(mfrow=c(2,3))
acf(post.1[,1], main="sigma2")
acf(post.1[,7], main="theta6")
acf(post.1[,101], main="theta100")
acf(post.eff[,1], main="sigma2")
acf(post.eff[,7], main="theta6")
acf(post.eff[,101], main="theta100")
