## Boston Exit Poll Data
##
## Data from Cobb, Greiner, and Quinn (2012)
## "Can Voter ID Laws Be Administered in a Race-Neutral Manner?
## Evidence from the City of Boston in 2008" QJPS
##
##
## Data from Cathedral High School Polling Location
##
set.seed(776829)
library(coda)

## Overall fraction asked for ID
n <- 314   ## number of voters sampled
y <- 86    ## number of sampled voters who report being asked for ID

## parameters for beta prior
a <- 1
b <- 1

## posterior is Beta(y+a, n-y+b)
post.mean.exact <- (y+a)/(n+a+b)
post.var.exact <- ((y+a)*(n-y+b))/( (n+a+b)^2 * (n+a+b+1) )
post.sd.exact <- sqrt(post.var.exact)
post.025.exact <- qbeta(0.025, y+a, n-y+b)
post.975.exact <- qbeta(0.975, y+a, n-y+b)

burnin <- 500 ## number of burn-in iterations
M <- 10000 ## number of iterations after burnin-in
tot.iter <- burnin + M ## number of total iterations
post.sim <- rep(NA, M) ## storage vector for the post-burn-in MCMC draws


## proposals are generated from Beta(prop.a, prop.b)
prop.a <- 1.0   
prop.b <- 1.0

theta <- 0.5 ## starting value

n.accepts <- 0 ## number of proposals that are accepted

## the Metropolis-Hastings sampling
for (iter in 1:tot.iter){

    ## generate candidate
    theta.can <- rbeta(1, shape1=prop.a, shape2=prop.b)

    ## calculate the acceptace probability
    log.piece1 <- dbeta(theta.can, y+a, n-y+b, log=TRUE)
    log.piece2 <- dbeta(theta, y+a, n-y+b, log=TRUE)
    log.piece3 <- dbeta(theta, prop.a, prop.b, log=TRUE)
    log.piece4 <- dbeta(theta.can, prop.a, prop.b, log=TRUE)
    
    ratio <- exp(log.piece1 - log.piece2 + log.piece3 - log.piece4)
    u <- runif(1)
    if (u < ratio){
        theta <- theta.can
        n.accepts <- n.accepts + 1
    }

    if (iter %% 1000 == 0){
        cat("iteration", iter, "of", tot.iter, "\n")
        cat("theta =", theta, "\n")
        cat("acceptance rate =", round(n.accepts/tot.iter, 3), "\n\n")
    }
    
    ## store results
    if (iter > burnin){
        post.sim[iter-burnin] <- theta
    }

    post.sim <- mcmc(post.sim, start=burnin+1, end=tot.iter)
    names(post.sim) <- "theta"
}
cat("acceptance rate =", round(n.accepts/tot.iter, 3), "\n")

## summarize posterior samples
post.mean.sim <- mean(post.sim)
post.var.sim <- var(post.sim)
post.sd.sim <- sd(post.sim)
post.025.sim <- quantile(post.sim, prob=0.025)
post.975.sim <- quantile(post.sim, prob=0.975)

## plotting the posterior density
## (exact and histogram of draws from the posterior)
hist(post.sim, nclass=50, prob=TRUE, xlim=c(0,1), col="orange", xlab="theta")
thetavals <- seq(from=0, to=1, by=.001)
lines(thetavals, dbeta(thetavals, y+a, n-y+b), lwd=3, col="blue")

cat("###################################################################\n")
cat("Summary of the Posterior Distribution (Exact)\n")
cat("posterior mean of theta:        ", round(post.mean.exact, 3), "\n", sep="")
cat("posterior variance of theta:    ", round(post.var.exact, 5), "\n", sep="")
cat("posterior SD of theta:          ", round(post.sd.exact, 3), "\n", sep="")
cat("central 95% credible interval: (", round(post.025.exact, 3), ", ",
    round(post.975.exact, 3), ")\n", sep="")
cat("###################################################################\n\n\n")

cat("###################################################################\n")
cat("Summary of the Posterior Distribution (Simulation)\n")
cat("posterior mean of theta:        ", round(post.mean.sim, 3), "\n", sep="")
cat("posterior variance of theta:    ", round(post.var.sim, 5), "\n", sep="")
cat("posterior SD of theta:          ", round(post.sd.sim, 3), "\n", sep="")
cat("central 95% credible interval: (", round(post.025.sim, 3), ", ",
    round(post.975.sim, 3), ")\n", sep="")
cat("###################################################################\n")
