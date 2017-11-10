## Castaneda v. Partida
## (based on examples 5.10.1 and 5.10.2 in DeGroot and Schervish
##     _Probability and Statistics_
##

n <- 220   ## number of grand jurors chosen
y <- 100   ## number of Mexican-American jurors chosen

theta.true <- 0.791  ## true Mexican-American fraction of population


## prior should be such that E[theta] = 0.791
## (no discrimination in expectation)
##
## E[theta] = a / (a + b)
## so we should only look at priors for which 
## a = 3.785 * b
##
## ( verify (3.785 * b) / (3.785 * b + b) = 0.791 for all b )
##

b <- seq(from=1, to=750, by=1)
a <- 3.785 * b
ab <- a+b

prior.025 <- rep(NA, length(b))
prior.975 <- rep(NA, length(b))
post.025 <- rep(NA, length(b))
post.975 <- rep(NA, length(b))

for (i in 1:length(b)){
  prior.025[i] <- qbeta(0.025, a[i], b[i])
  prior.975[i] <- qbeta(0.975, a[i], b[i])
  post.025[i] <- qbeta(0.025, a[i]+y, b[i]+n-y)
  post.975[i] <- qbeta(0.975, a[i]+y, b[i]+n-y)

}


plot(ab, prior.025, type="n", lty=2, col="orange", ylim=c(0,1), ylab="theta",
     xlab="a + b")

polygon(c(ab, rev(ab)), c(prior.025, rev(prior.975)), col=rgb(1,.647,0,.65),
        border=NA)
polygon(c(ab, rev(ab)), c(post.025, rev(post.975)), col=rgb(0,0,1,.65),
        border=NA)

abline(h=theta.true, col="white", lwd=1.25)
legend(x=1500, y=0.4, legend=c("Prior 95% CI", "Posterior 95% CI"),
       fill=c(rgb(1,.647,0,.65), rgb(0,0,1,.65))) 
