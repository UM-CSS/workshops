## multinomial mixture analysis of the Federalist authorship problem
##
## Fits 2-component multinomial mixture model to stop word frequencies of
## of Federalist papers with the constraint that the known Hamilton and
## Madison essays are in separate mixture components with certainty.
## Assumes that only Hamilton or Madison could have been the author of a
## disputed essay.
##
## Kevin Quinn
## 10/24/2017
##
set.seed(928472)

## load MCMCpack (needed for Dirichlet distribution)
library(MCMCpack)


## load Hamilton data
ham <- read.csv("FP-stopwords-Hamilton.csv", stringsAsFactors=FALSE)
##      we're assuming independence so we can sum the word counts
ham.vec <- apply(ham[,2:127], 2, sum)
## number of known Hamilton essays
ham.n <- nrow(ham)

## load Madison data
mad <- read.csv("FP-stopwords-Madison.csv", stringsAsFactors=FALSE)
##      we're assuming independence so we can sum the word counts
mad.vec <- apply(mad[,2:127], 2, sum)
## number of known Madison essays
mad.n <- nrow(mad)


## load disputed authorship data
dis <- read.csv("FP-stopwords-disputed.csv", stringsAsFactors=FALSE)
## extract the word counts from the disputed essays
dis.counts <- dis[,2:127]
## number of disputed essays
dis.n <- nrow(dis)

## prior parameter for the Dirichlet distributions (assumed same for both)
a <- rep(1/126, 126)

## prior paramaters for beta distribution governing pi (marginal probability
## of a randomly chosen essay being authored by Hamilton)
b <- 0.5
c <- 0.5

## number of burn-in iterations
burnin <- 100
## number of Monte Carlo samples
m <- 200
## total iterations
tot.iter <- burnin + m


## starting values
theta.ham <- ham.vec/sum(ham.vec)
theta.mad <- mad.vec/sum(mad.vec)
pi <- 0.5
post.prob.ham <- rep(NA, dis.n)
z <- rep(NA, dis.n)



## storage for the posterior draws
theta.ham.store <- matrix(NA, m, length(theta.ham))
theta.mad.store <- matrix(NA, m, length(theta.mad))
pi.store <- rep(NA, m)
post.prob.ham.store <- matrix(NA, m, length(post.prob.ham))
z.store <- matrix(NA, m, length(z))



## the Gibbs sampling
for (iter in 1:tot.iter){

    ## sample the latent cluster memberships, z, for the disputed docs
    ## z_i = 1 if essay i is Hamilton, 0 if essay i is Madison
    for (i in 1:dis.n){
        dis.i <- dis.counts[i,]
        prob.ham.i <- (pi * dmultinom(dis.i, prob=theta.ham))/
            (pi * dmultinom(dis.i, prob=theta.ham) +
             (1-pi) * dmultinom(dis.i, prob=theta.mad))
        z[i] <- rbinom(1, 1, prob.ham.i)
        post.prob.ham[i] <- prob.ham.i
    }


    ## sample theta.ham
    holder <- dis.counts[z==1,]
    data <- apply(holder, 2, sum) + ham.vec
    theta.ham <- rdirichlet(1, a + data)
    

    ## sample theta.mad
    holder <- dis.counts[z==0,]
    data <- apply(holder, 2, sum) + mad.vec
    theta.mad <- rdirichlet(1, a + data)


    ## sample pi
    pi <- rbeta(1, sum(z) + ham.n + b,  sum(z==0) + mad.n + c)
    

    if (iter %% 20 == 0){
        cat("iteration", iter, "of", tot.iter, "\n")
    }

    if (iter > burnin){
        theta.ham.store[iter-burnin,] <- theta.ham
        theta.mad.store[iter-burnin,] <- theta.mad
        pi.store[iter-burnin] <- pi
        z.store[iter-burnin,] <- z
        post.prob.ham.store[iter-burnin,] <- post.prob.ham
    }

    
} ## end Gibbs sampling loop






## calculate posterior probability of disputed essay's authorship
cat("\n", "Essay              E[Pr(Hamilton|y)]         95% CI \n")
cat("------------------------------------------------------------\n")
for (i in 1:nrow(dis)){
    essay <- dis$X[i]
    cat(essay, "   Pr(Hamilton) = ",
        signif(mean(post.prob.ham.store[,i]), 3),
        "   (",
        signif(quantile(post.prob.ham.store[,i], prob=0.025), 3),
        ", ",
        signif(quantile(post.prob.ham.store[,i], prob=0.975), 3),
        ") ", "\n", sep="")
}
cat("------------------------------------------------------------\n")

