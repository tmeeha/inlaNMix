# set parameters
nrep.max <- 5 # J surveys
n <- 75 # I sites
b0 <- 0.5 # lambda intercept
b1 <- 2.0 # effect of x1 on lambda
a0 <- 1.0 # p intercept
a2 <- 0.5 # effect of x2 on p
size <- 3.0 # size of theta
overdispersion <- 1 / size # for negative binomial

# make empty vectors and matrix
x1 <- c(); x2 <- c(); lambdas <- c(); Ns <- c(); y <- matrix(NA, n, nrep.max)

# fill vectors and matrix
for(i in 1:n) {
  x1.i <- runif(1) - 0.5
  lambda <- exp(b0 + b1 * x1.i)
  N <- rnbinom(1, mu = lambda, size = size)
  x2.i <- runif(1) - 0.5
  eta <- a0 + a2 * x2.i
  p <- exp(eta) / (exp(eta) + 1)
  nr <- sample(1:nrep.max, 1)
  y[i, 1:nr] <- rbinom(nr, size = N, prob = p)
  x1 <- c(x1, x1.i); x2 <- c(x2, x2.i)
  lambdas <- c(lambdas, lambda); Ns <- c(Ns, N)
}

# bundle counts, lambda intercept, lambda cov
library(INLA)
Y <- inla.mdata(y, 1, x1)
idx1 <- 1:n

# run inla and summarize output
result <- inla(Y ~ 1 + x2 + f(idx1, model='iid'),
         data = list(Y=Y, x2=x2, idx1=idx1),
         family = "nmixnb",
         control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.01,
                              prec.intercept = 0.01),
         control.family = list(hyper = list(theta1 = list(param = c(0, 0.01)),
                                            theta2 = list(param = c(0, 0.01)),
                                            theta3 = list(prior = "flat",
                                                          param = numeric()))),
         control.compute=list(config = TRUE, waic=T))
summary(result)

# get and evaluate fitted values
source("inla.nmix.lambda.fitted.R")
lam.fits <- inla.nmix.lambda.fitted(result, 500)$fitted.summary
plot(lam.fits$median.lambda, lambdas)
round(sum(lam.fits$median.lambda), 0); sum(Ns)
