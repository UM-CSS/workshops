load("ClusteredBinomialData.Rda")

source("MCMChierBetaBinom.R")

out <- MCMChierBetaBinom(m=data$n.respondents, y=data$n.approve,
                         a=0.01, b=0.01, method="marginal",
                         burnin=5000, thin=1, mcmc=10000,
                         tune=1.5, verbose=500)


