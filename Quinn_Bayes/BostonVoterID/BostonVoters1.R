## Boston Exit Poll Data
##
## Data from Cobb, Greiner, and Quinn (2012)
## "Can Voter ID Laws Be Administered in a Race-Neutral Manner?
## Evidence from the City of Boston in 2008" QJPS
##
##
## Data from Cathedral High School Polling Location
##


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


M <- 50000 ## number of Monte Carlo simulations
post.sim <- rbeta(M, y+a, n-y+b)
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
