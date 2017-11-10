## Boston Exit Poll Data
##
## Data from Cobb, Greiner, and Quinn (2012)
## "Can Voter ID Laws Be Administered in a Race-Neutral Manner?
## Evidence from the City of Boston in 2008" QJPS
##
##
## Data from Cathedral High School Polling Location
##



## black fraction asked for ID
n.b <- 56    ## number of black voters sampled
y.b <- 20    ## number of sampled black voters who report being asked for ID

## white fraction asked for ID
n.w <- 199   ## number of white voters sampled
y.w <- 42    ## number of sampled white voters who report being asked for ID

## parameters for beta prior for theta.b
a.b <- 1
b.b <- 1

## parameters for beta prior for theta.w
a.w <- 1
b.w <- 1


M <- 100000 ## number of Monte Carlo simulations

## posterior for theta.b is Beta(y.b+a.b, n.b-y.b+b.b)
post.sim.b <- rbeta(M, y.b+a.b, n.b-y.b+b.b)
post.mean.sim.b <- mean(post.sim.b)
post.var.sim.b <- var(post.sim.b)
post.sd.sim.b <- sd(post.sim.b)
post.025.sim.b <- quantile(post.sim.b, prob=0.025)
post.975.sim.b <- quantile(post.sim.b, prob=0.975)

## posterior for theta.w is Beta(y.w+a.w, n.w-y.w+b.w)
post.sim.w <- rbeta(M, y.w+a.w, n.w-y.w+b.w)
post.mean.sim.w <- mean(post.sim.w)
post.var.sim.w <- var(post.sim.w)
post.sd.sim.w <- sd(post.sim.w)
post.025.sim.w <- quantile(post.sim.w, prob=0.025)
post.975.sim.w <- quantile(post.sim.w, prob=0.975)



cat("\n\n###################################################################\n")
cat("Summary of the Posterior Distribution (Blacks)\n")
cat("posterior mean of theta:        ", round(post.mean.sim.b, 3), "\n", sep="")
cat("posterior variance of theta:    ", round(post.var.sim.b, 5), "\n", sep="")
cat("posterior SD of theta:          ", round(post.sd.sim.b, 3), "\n", sep="")
cat("central 95% credible interval: (", round(post.025.sim.b, 3), ", ",
    round(post.975.sim.b, 3), ")\n", sep="")
cat("###################################################################\n\n")

cat("###################################################################\n")
cat("Summary of the Posterior Distribution (Whites)\n")
cat("posterior mean of theta:        ", round(post.mean.sim.w, 3), "\n", sep="")
cat("posterior variance of theta:    ", round(post.var.sim.w, 5), "\n", sep="")
cat("posterior SD of theta:          ", round(post.sd.sim.w, 3), "\n", sep="")
cat("central 95% credible interval: (", round(post.025.sim.w, 3), ", ",
    round(post.975.sim.w, 3), ")\n", sep="")
cat("###################################################################\n\n")


cat("Pr(theta.black > theta.white) = ",
    round(mean(post.sim.b > post.sim.w), 3), "\n\n")
