## simple analysis of the Federalist authorship problem
##
## Fits multinomial models to stop word frequencies of known Hamilton and
## Madison essays and then calculates the posterior probability that each
## of the disputed essays was a Hamilton essay. Assumes that only Hamilton
## or Madison could have been the author of a disputed essay.
##
## Kevin Quinn
## 10/24/2017
##
set.seed(583104)

## load MCMCpack (needed for Dirichlet distribution)
library(MCMCpack)


## load Hamilton data
ham <- read.csv("FP-stopwords-Hamilton.csv", stringsAsFactors=FALSE)

## load Madison data
mad <- read.csv("FP-stopwords-Madison.csv", stringsAsFactors=FALSE)

## load disputed authorship data
dis <- read.csv("FP-stopwords-disputed.csv", stringsAsFactors=FALSE)



## prior parameter for the Dirichlet distributions (assumed same for both)
a <- rep(1/126, 126)

## prior probabilities of authorship
prior.ham <- 0.5
prior.mad <- 0.5

## number of Monte Carlo samples
m <- 500


## fit multinomial model to Hamilton data and save m draws from posterior
##      we're assuming independence so we can sum the word counts
ham.vec <- apply(ham[,2:127], 2, sum)
ham.post.draws <- rdirichlet(m, a+ham.vec)



## fit multinomial model to Madison data and save m draws from posterior
##      we're assuming independence so we can sum the word counts
mad.vec <- apply(mad[,2:127], 2, sum)
mad.post.draws <- rdirichlet(m, a+mad.vec)





## calculate posterior probability of disputed essay's authorship
cat("\n", "Essay              E[Pr(Hamilton|y)]         95% CI \n")
cat("------------------------------------------------------------\n")
for (i in 1:nrow(dis)){
    dis.i <- dis[i,]
    essay <- dis.i$X
    data.i <- dis.i[2:127]
    prob.ham.draws.i <- rep(NA, m)
    for (j in 1:m){
        theta.ham <- ham.post.draws[j,]
        theta.mad <- mad.post.draws[j,]
        prob.ham.draws.i[j] <- (prior.ham * dmultinom(data.i, prob=theta.ham))/
            (prior.ham * dmultinom(data.i, prob=theta.ham) +
             prior.mad * dmultinom(data.i, prob=theta.mad))
    }
    cat(essay, "   Pr(Hamilton) = ",
        signif(mean(prob.ham.draws.i), 3),
        "   (",
        signif(quantile(prob.ham.draws.i, prob=0.025), 3),
        ", ",
        signif(quantile(prob.ham.draws.i, prob=0.975), 3),
        ") ", "\n", sep="")
}
cat("------------------------------------------------------------\n")

