source("BivariateNormalModel.R")
source("MissingBivariateNormalModel.R")

load("fulldata.Rda")

out.full <- MCMCBivariateNormal(Y=fulldata, m=c(0,0),
                                M=matrix(c(.01,0,0,.01), 2, 2, byrow=TRUE),
                                s=5, S=matrix(c(5,0,0,5), 2, 2, byrow=TRUE),
                                burnin=5000, ngibbs=10000)




load("MCARdata.Rda")

out.mcar <- MCMCMissingBivariateNormal(Y=MCARdata, m=c(0,0),
                                M=matrix(c(.01,0,0,.01), 2, 2, byrow=TRUE),
                                s=5, S=matrix(c(5,0,0,5), 2, 2, byrow=TRUE),
                                burnin=5000, ngibbs=10000)


load("MARdata.Rda")

out.mar <- MCMCMissingBivariateNormal(Y=MARdata, m=c(0,0),
                                M=matrix(c(.01,0,0,.01), 2, 2, byrow=TRUE),
                                s=5, S=matrix(c(5,0,0,5), 2, 2, byrow=TRUE),
                                burnin=5000, ngibbs=10000)




load("NIdata.Rda")

out.ni <- MCMCMissingBivariateNormal(Y=NIdata, m=c(0,0),
                                M=matrix(c(.01,0,0,.01), 2, 2, byrow=TRUE),
                                s=5, S=matrix(c(5,0,0,5), 2, 2, byrow=TRUE),
                                burnin=5000, ngibbs=10000)






