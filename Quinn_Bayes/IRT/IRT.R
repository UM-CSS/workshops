library(MCMCpack)


data(SupremeCourt)

posterior1 <- MCMCirt1d(t(SupremeCourt),
                        theta.constraints=list(Scalia="+", Ginsburg="-"),
                        B0.alpha=.2, B0.beta=.2,
                        burnin=500, mcmc=100000, thin=20, verbose=500,
                        store.item=TRUE)

geweke.diag(posterior1)
plot(posterior1)
summary(posterior1)
