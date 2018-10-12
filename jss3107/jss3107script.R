


# code for Meehan, Michel, and Rue #############################################
# demonstration of R-INLA for analyzing N-mixture models
# 08/24/2017
# ##############################################################################



# directories, libraries, and helpers ##########################################
# uncomment the next two installation steps as needed
# install.packages(c("runjags", "unmarked", "ggplot2"))
# install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")
library(runjags)
library(INLA) # may require testing version
library(unmarked)
library(ggplot2)
set.seed(12345)

# multiplot function
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} # end

# plotting theme
theme_acbs <- function (base_size = 11, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey20"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = rel(0.9), angle = 0),
          axis.text.y = element_text(size = rel(0.9), angle = 0),
          strip.background = element_rect(fill = "grey85", colour = "grey20"),
          legend.key = element_rect(fill = "white", colour = NA),
          complete = TRUE)
} # end

# data generating function
sim.nmix <- function(n.sites = 72, # number of study sites
                     n.surveys = 3,                   # short term replicates
                     n.years = 9,                     # odd number years
                     b0 = 2.0,                        # intercept log lambda, b0
                     b1 = 2.0,                        # slope log lambda, x1
                     b2 = -3.0,                       # slope log lambda, x2
                     b3 = 1.0,                        # slope log lambda, x3
                     a0 = 1.0,                        # intercept logit p
                     a1 = -2.0,                       # slope logit p, x1
                     a4 = 1.0,                        # slope logit p, x4
                     th = 3.0                         # size of overdisperison
){

  # make empty N and Y arrays
  if(n.years %% 2 == 0) {n.years <- n.years + 1}
  N.tr <- array(dim = c(n.sites, n.years))
  Y <- array(dim = c(n.sites, n.surveys, n.years))
  Y.m <- array(dim = c(n.sites, n.surveys, n.years))

  # create abundance covariate values
  x1 <- array(as.numeric(scale(runif(n = n.sites, -0.5, 0.5), scale = F)),
              dim = c(n.sites, n.years))
  x2 <- array(as.numeric(scale(runif(n = n.sites, -0.5, 0.5), scale = F)),
              dim = c(n.sites, n.years))
  yrs <- 1:n.years; yrs <- (yrs - mean(yrs)) / (max(yrs - mean(yrs))) / 2
  x3 <- array(rep(yrs, each = n.sites), dim = c(n.sites, n.years))

  # fill true N array
  lam.tr <- exp(b0 + b1 * x1 + b2 * x2 + b3 * x3)
  for(i in 1:n.sites){
    for(k in 1:n.years){
      N.tr[i, k] <- rnbinom(n = 1, mu = lam.tr[i, k], size = th)
    }}

  # create detection covariate values
  x1.p <- array(x1[,1], dim = c(n.sites, n.surveys, n.years))
  x4 <- array(as.numeric(scale(runif(n = n.sites * n.surveys * n.years,
    -0.5, 0.5), scale = F)), dim = c(n.sites, n.surveys, n.years))

  # average x4 per site-year for ex 1
  x4.m <- apply(x4, c(1, 3), mean, na.rm = F)
  out1 <- c()
  for(k in 1:n.years){
    chunk1 <- x4.m[ , k]
    chunk2 <- rep(chunk1, n.surveys)
    out1 <- c(out1, chunk2)
  }
  x4.m.arr <- array(out1, dim = c(n.sites, n.surveys, n.years))

  # fill Y.m count array using site-year x4.m for ex 1
  p.tr1 <- plogis(a0 + a1 * x1.p + a4 * x4.m.arr)
  for (i in 1:n.sites){
    for (k in 1:n.years){
      for (j in 1:n.surveys){
        Y.m[i, j, k] <- rbinom(1, size = N.tr[i, k], prob = p.tr1[i, j, k])
      }}}

  # fill Y count array using site-surv-year x4 for ex 2
  p.tr2 <- plogis(a0 + a1 * x1.p + a4 * x4)
  for (i in 1:n.sites){
    for (k in 1:n.years){
      for (j in 1:n.surveys){
        Y[i, j, k] <- rbinom(1, size = N.tr[i, k], prob = p.tr2[i, j, k])
      }}}

  # format Y.m for dataframe
  Y.m.df <- Y.m[ , , 1]
  for(i in 2:n.years){
    y.chunk <- Y.m[ , , i]
    Y.m.df <- rbind(Y.m.df, y.chunk)
  }

  # format covariates for dataframe
  x1.df <- rep(x1[ , 1], n.years)
  x2.df <- rep(x2[ , 1], n.years)
  x3.df <- rep(x3[1, ], each = n.sites)
  x1.p.df <- rep(x1.p[ , 1, 1], n.years)
  x4.df <- c(x4.m)

  # make dataframes for output
  inla.df <- unmk.df <- data.frame(y1 = Y.m.df[ , 1], y2 = Y.m.df[ , 2],
                                   y3 = Y.m.df[ , 3], x1 = x1.df, x2 = x2.df, x3 = x3.df,
                                   x1.p = x1.p.df, x4.m = x4.df)

  # return data
  return(list(inla.df = inla.df, unmk.df = unmk.df, n.sites = n.sites,
              n.surveys = n.surveys, n.years = n.years, x1 = x1[ , 1],
              x2 = x2[ , 1], x3 = x3[1, ], x4 = x4, x4.m = x4.m, x4.m.arr = x4.m.arr,
              Y = Y, Y.m = Y.m, lam.tr = lam.tr, N.tr = N.tr, x1.p = x1.p[ , 1, 1]
  ))
} # end

# inla.nmix.lambda.fitted function
inla.nmix.lambda.fitted <- function(result, sample.size=1000,
                                    return.posteriors=FALSE,
                                    scale="exp"){
  # Checks and warnings
  library(INLA)
  fam <- result$.args$family
  if(length(grep(pattern = "nmix", x = fam)) == 0) {
    stop("This function is only for models with 'nmix' or 'nmixnb' likelihoods")
  }
  if(missing(result)) stop("Please specify a model result")
  s.check <- as.numeric(scale == "exp") + as.numeric(scale == "log")
  if(s.check == 0) stop("Scale must be set to 'exp' or 'log'")
  if(sample.size < 500) warning("Please increase the sample size")

  # Get counts and lambda covariates from 'inla.mdata' object
  mdata.obj <- result$.args$data[[1]]
  counts <- as.data.frame(mdata.obj[grep(pattern="Y", names(mdata.obj))])
  lambda.covs <- as.data.frame(mdata.obj[grep(pattern="X", names(mdata.obj))])
  lambda.covs <- as.matrix(lambda.covs)
  n.counts <- ncol(counts)
  n.lambda.covs <- ncol(lambda.covs)
  n.data <- nrow(counts)

  # Get samples from hyperpars
  hyperpar.samples <- inla.hyperpar.sample(sample.size, result)
  s.names <- rownames(hyperpar.samples)

  # Discard hyperpar.samples columns that are not nmix
  hyperpar.samples <- hyperpar.samples[ , grep("beta",
                                               colnames(hyperpar.samples))]
  n.samp.covs <- ncol(hyperpar.samples)
  if(n.lambda.covs != n.samp.covs) {
    stop("This function can not handle multiple 'nmix' components.")
  }

  # Combine lambda covariates and hyperparameter posteriors
  fitted.posteriors <- matrix(-1.0000, nrow=n.data, ncol=sample.size)
  # For each site or site-by-year combination
  for(i in 1:n.data){
    obs <- lambda.covs[i,]
    # For each sample from the hyperparameter posterior
    for(j in 1:sample.size){
      post <- hyperpar.samples[j, ]
      fitted <- sum(obs * post)
      fitted.posteriors[i,j] <- fitted
    }
  }

  # Clean up the resulting matrix
  index <- 1:n.data
  fitted.posteriors <- as.data.frame(fitted.posteriors)
  row.names(fitted.posteriors) <- NULL
  names(fitted.posteriors) <- s.names

  # Create posterior summaries of fitted values
  if(scale=="exp"){
    fitted.posteriors <- exp(fitted.posteriors)
  }
  fitted.meds <- round(apply(fitted.posteriors, 1, median), 4)
  fitted.means <- round(apply(fitted.posteriors, 1, mean), 4)
  fitted.sds <- round(apply(fitted.posteriors, 1, sd), 4)
  fitted.q025 <- round(apply(fitted.posteriors, 1, quantile, probs=0.025), 4)
  fitted.q500 <- round(apply(fitted.posteriors, 1, quantile, probs=0.500), 4)
  fitted.q975 <- round(apply(fitted.posteriors, 1, quantile, probs=0.975), 4)
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  fitted.modes <- round(apply(fitted.posteriors, 1, Mode), 4)
  fitted.summary <- data.frame(mean.lambda=fitted.means, sd.lambda=fitted.sds,
                               quant025.lambda=fitted.q025,
                               median.lambda=fitted.q500,
                               quant975.lambda=fitted.q975,
                               mode.lambda=fitted.modes)
  fitted.summary <- cbind(index, fitted.summary)
  fitted.posteriors <- cbind(index, fitted.posteriors)

  # Create returned object
  if(return.posteriors==TRUE){
    out <- list(fitted.summary=fitted.summary,
                fitted.posteriors=fitted.posteriors)
  } else {
    out <- list(fitted.summary=fitted.summary)
  }
  return(out)

} # end
# ##############################################################################



# create dataset for examples 1 and 2 ##########################################
set.seed(12345)
sim.data <- sim.nmix()
# ##############################################################################



# inla analysis, example 1 #####################################################
# view data
str(sim.data$inla.df, digits.d=2)
round(head(sim.data$inla.df),3)

# make inla.mdata object
inla.data <- sim.data$inla.df
y.mat <- as.matrix(inla.data[,c("y1", "y2", "y3")])
counts.and.count.covs <- inla.mdata(y.mat, 1,
                                    inla.data$x1, inla.data$x2,
                                    inla.data$x3)
str(counts.and.count.covs)

# run inla model
out.inla.1 <- inla(counts.and.count.covs ~ 1 + x1.p + x4.m,
  data = list(counts.and.count.covs = counts.and.count.covs,
             x1.p = inla.data$x1.p, x4.m = inla.data$x4.m),
  family = "nmixnb",
  control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.01,
                      prec.intercept = 0.01),
  control.family = list(hyper = list(theta1 = list(param = c(0, 0.01)),
                                    theta2 = list(param = c(0, 0.01)),
                                    theta3 = list(param = c(0, 0.01)),
                                    theta4 = list(param = c(0, 0.01)),
                                    theta5 = list(prior = "flat",
                                                  param = numeric()))),
  verbose = TRUE,
  control.compute=list(config = TRUE))

# view summary
summary(out.inla.1, digits = 3)

# get fitted values
out.inla.1.lambda.fits <- inla.nmix.lambda.fitted(result = out.inla.1,
  sample.size = 5000, return.posteriors = FALSE)$fitted.summary
head(out.inla.1.lambda.fits)

# compare fitted and true lambdas
summary(out.inla.1.lambda.fits$median.lambda)
summary(c(sim.data$lam.tr))
cor(out.inla.1.lambda.fits$median.lambda, c(sim.data$lam.tr))
sum(out.inla.1.lambda.fits$median.lambda)
sum(c(sim.data$lam.tr))
sum(c(sim.data$N.tr))
# ##############################################################################



# jags analysis, example 1 #####################################################
jags.model.string <- "
model {
# priors
a0 ~ dnorm(0, 0.01)       # detection intercept
a1 ~ dnorm(0, 0.01)       # detection cov 1 effect
a4 ~ dnorm(0, 0.01)       # detection cov 3 effect
b0 ~ dnorm(0, 0.01)       # lambda intercept
b1 ~ dnorm(0, 0.01)       # lambda cov 1 effect
b2 ~ dnorm(0, 0.01)       # lambda cov 2 effect
b3 ~ dnorm(0, 0.01)       # year effect
th ~ dunif(0, 5)          # overdispersion size
# abundance component
for (k in 1:n.years){
for (i in 1:n.sites){
N[i, k] ~ dnegbin(prob[i, k], th) # negative binomial
prob[i, k] <- th / (th + lambda[i, k]) # overdispersion
log(lambda[i, k]) <- b0 + (b1 * x1[i]) +
(b2 * x2[i]) + (b3 * x3[k])
# detection component
for (j in 1:n.surveys){
Y.m[i, j, k] ~ dbin(p[i,j,k], N[i,k])
p[i, j, k] <- exp(lp[i,j,k]) / (1 + exp(lp[i,j,k]))
lp[i, j, k] <- a0 + (a1 * x1.p[i]) + (a4 * x4.m[i, k])
} # close j loop
} # close i loop
} # close k loop
} # close model loop
"

# parameters to monitor
params <- c("a0", "a1", "a4", "b0", "b1","b2",
            "b3", "th")

# jags data
jags.data <- list(Y.m = sim.data$Y.m, x1 = sim.data$x1,
                  x2 = sim.data$x2, x3 = sim.data$x3,
                  x1.p = sim.data$x1.p, x4.m = sim.data$x4.m,
                  n.sites = sim.data$n.sites, n.surveys = sim.data$n.surveys,
                  n.years = sim.data$n.years)

# initial values
N.init <- sim.data$Y.m # initial count values
N.init[is.na(N.init)] <- 1 # clean up NA's
N.init <- apply(N.init, c(1, 3), max) + 1 # zero values cause trouble
inits <- function() list(N = N.init, b0 = rnorm(1, 0, 0.01),
                         a0 = rnorm(1, 0, 0.01), a1 = rnorm(1, 0, 0.01),
                         b2 = rnorm(1,0,0.01), b1 = rnorm(1, 0, 0.01),
                         a4 = rnorm(1, 0, 0.01), b3 = rnorm(1, 0, 0.01),
                         th = runif(1, 0.5, 2.5))

# set run parameters
nc <- 3; na <- 3000; nb <- 3000; ni <- 6000; nt <- 10

# run jags
ptm.1 <- proc.time()
out.jags.1 <- run.jags(model = jags.model.string, data = jags.data,
                       monitor = params, n.chains = nc, inits = inits,
                       burnin = nb, adapt = na, sample = ni, thin = nt,
                       modules = "glm on", method = "parallel")

# view summary
round(summary(out.jags.1), 3)[ , c(1:5, 9, 11)]

# computing time
jags.time.1 <- proc.time() - ptm.1
round(jags.time.1, 2)[3]
# ##############################################################################



# unmarked analysis, example 1##################################################
# format count data
y.unmk <- as.matrix(sim.data$unmk.df[,c("y1", "y2", "y3")])

# format site covariates
site.covs.unmk <- data.frame(x1 = sim.data$unmk.df[ , "x1"],
  x2 = sim.data$unmk.df[ , "x2"], x3 = sim.data$unmk.df[ , "x3"])

# format observation covariates
obs.covs.unmk <- list(x1.p = matrix(sim.data$unmk.df[ , "x1.p"],
  sim.data$n.sites * sim.data$n.years, sim.data$n.surveys),
  x4.m = matrix(sim.data$unmk.df[ , "x4.m"],
  sim.data$n.sites * sim.data$n.years,
  sim.data$n.surveys))

# make unmarked pcount frame
unmk.data <- unmarkedFramePCount(y = y.unmk, siteCovs = site.covs.unmk,
                                 obsCovs = obs.covs.unmk)

# run unmk model
ptm <- proc.time()
out.unmk.1 <- pcount(~ 1 + x1.p + x4.m
                     ~ 1 + x1 + x2 + x3,
                     data = unmk.data, mixture = "NB")

# view summary
summary(out.unmk.1)

# computing time
round(unmk.time.1 <- proc.time() - ptm, 2)[3]
# ##############################################################################



# explore results, example 1 ###################################################
# get jags out
jags.out.tab <- round(as.data.frame(summary(out.jags.1)[,c(2, 1, 3)]), 3)

# get unmk out
e1 <- as.matrix(out.unmk.1@estimates@estimates$det@estimates)
ci1 <- as.matrix(confint(out.unmk.1, type="det"))
e2 <- as.matrix(out.unmk.1@estimates@estimates$state@estimates)
ci2 <- as.matrix(confint(out.unmk.1, type="state"))
e3 <- as.matrix(out.unmk.1@estimates@estimates$alpha@estimates)
ci3 <- as.matrix(confint(out.unmk.1, type="alpha"))
e4 <- rbind(e1,e2,e3)
ci4 <- rbind(ci1,ci2,ci3)
unmk.out.tab <- cbind(e4,ci4)

# get inla out
in1 <- out.inla.1$summary.fixed[,c(4, 3, 5)]
in2 <- out.inla.1$summary.hyperpar[,c(4, 3, 5)][-5,]
in3 <- 1/out.inla.1$summary.hyperpar[,c(4, 3, 5)][5,]
in3 <- in3[,c(1,3,2)]
inla.out.tab <- rbind(in1,in2,in3)

# combine out
all.out.tab <- round(cbind(jags.out.tab, unmk.out.tab, inla.out.tab),2)
colnames(all.out.tab) <- c("JAGS Median","JAGS Lower CrL", "JAGS Upper CrL",
                           "UNMK Estimate","UNMK Lower CL", "UNMK Upper CL",
                           "INLA Median","INLA Lower CrL", "INLA Upper CrL")
rownames(all.out.tab) <- c("P Int", "P Cov 1", "P Cov 3", "Lambda Int",
                           "Lambda Cov 1", "Lambda Cov 2", "Lambda Year",
                           "Overdispersion")

# get output from jags
jags.mcmc <- combine.mcmc(out.jags.1, thin=3, return.samples=6000)
jags.df <- as.data.frame(jags.mcmc)

# make density for b0
d1 <- as.data.frame(
  out.inla.1$marginals.hyperpar$`beta[1] for NMix observations`)
d1$x <- (d1$x)
d2 <- as.data.frame(density(jags.df$b0)[c(1,2)])
d2$x <- (d2$x)
col1 <- 4
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- (c(d3,d4,d5))
d7 <- data.frame(y=c(0,0,0), x=d6)
p1 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab(expression(beta[0])) +
  ylab(" ") + scale_y_continuous(breaks=c(3,6,9)) + theme_acbs() +
  geom_vline(xintercept=(2), col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) +
  theme(axis.title.x = element_text(size=12)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))

# make density for b1
d1 <- as.data.frame(
  out.inla.1$marginals.hyperpar$`beta[2] for NMix observations`)
d2 <- as.data.frame(density(jags.df$b1)[c(1,2)])
col1 <- 5
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p2 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab(expression(beta[1]))  +
  ylab("Posterior density") + theme_acbs() +
  geom_vline(xintercept=2, col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) +
  theme(axis.title.x = element_text(size=12)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))

# make density for b2
d1 <- as.data.frame(
  out.inla.1$marginals.hyperpar$`beta[3] for NMix observations`)
d2 <- as.data.frame(density(jags.df$b2)[c(1,2)])
col1 <- 6
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p3 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab(expression(beta[2]))  +
  ylab(" ") + theme_acbs() + theme(axis.title.x = element_text(size=12)) +
  geom_vline(xintercept=-3, col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))

# make density for b3
d1 <- as.data.frame(
  out.inla.1$marginals.hyperpar$`beta[4] for NMix observations`)
d2 <- as.data.frame(density(jags.df$b3)[c(1,2)])
col1 <- 7
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p4 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab(expression(beta[3]))  +
  ylab(NULL) + theme_acbs() + theme(axis.title.x = element_text(size=12)) +
  geom_vline(xintercept=1, col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))

# make density for a0
d1 <- as.data.frame(out.inla.1$marginals.fixed$`(Intercept)`)
d1$x <- (d1$x)
d2 <- as.data.frame(density(jags.df$a0)[c(1,2)])
d2$x <- (d2$x)
col1 <- 1
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- (c(d3,d4,d5))
d7 <- data.frame(y=c(0,0,0), x=d6)
p5 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab(expression(alpha[0]))  +
  ylab(NULL) + theme_acbs() + geom_vline(xintercept=(1), col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) +
  theme(axis.title.x = element_text(size=12)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))

# make density for a1
d1 <- as.data.frame(out.inla.1$marginals.fixed$x1)
d2 <- as.data.frame(density(jags.df$a1)[c(1,2)])
col1 <- 2
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p6 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab(expression(alpha[1]))  + ylab(NULL) + theme_acbs() +
  geom_vline(xintercept=-2, col="gray20", lty=1)+
  geom_line(data=d7, aes(y=y, x=x)) +
  theme(axis.title.x = element_text(size=12)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))

# make density for a4
d1 <- as.data.frame(out.inla.1$marginals.fixed$x4)
d2 <- as.data.frame(density(jags.df$a4)[c(1,2)])
col1 <- 3
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- c(d3,d4,d5)
d7 <- data.frame(y=c(0,0,0), x=d6)
p7 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab(expression(alpha[4]))  + ylab(NULL) + theme_acbs() +
  geom_vline(xintercept=1, col="gray20", lty=1)+
  geom_line(data=d7, aes(y=y, x=x)) +
  theme(axis.title.x = element_text(size=12)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))

# make density for overdisp
d1 <- as.data.frame(
  out.inla.1$marginals.hyperpar$`overdispersion for NMix observations`)
d1$x <- 1 / d1$x; d1$y <- d1$y / 8
d2 <- as.data.frame(density(jags.df$th)[c(1,2)])
col1 <- 8
d3 <- all.out.tab[col1, 5]
d4 <- all.out.tab[col1, 4]
d5 <- all.out.tab[col1, 6]
d6 <- exp(c(d3,d4,d5))
d7 <- data.frame(y=c(0,0,0), x=d6)
p8 <- ggplot(data=d2, aes(x=x, y=y)) + geom_path(col="gray60", lty=1) +
  geom_path(data=d1, aes(x=x, y=y), col="black", lty=2 ) +
  xlab(expression(theta)) + ylab(NULL) + theme_acbs() +
  geom_vline(xintercept=3, col="gray20", lty=1) +
  geom_line(data=d7, aes(y=y, x=x)) +
  theme(axis.title.x = element_text(size=12)) +
  geom_point(data=d7, aes(y=y, x=x), col="black", shape=c(32,19,32))

# plot all
png("fig1.png", width = 6, height = 5, units = 'in', res = 600)
multiplot(p1,p2,p3,p4,p5,p6,p7,p8,cols=3)
dev.off()
# ##############################################################################



# simulation code, example 2 ###################################################
sim.fun <- function(){

  # same inputs as before for most parameters
  n.sites = 72; n.surveys = 3; n.years = 9
  b0 = 2.0; b1 = 2.0; b2 = -3.0; b3 = 1.0
  a0 = 1.0; a1 = -2.0
  th = 3.0

  # now vary the effect size for detection covariate 3
  a4 = runif(1, -3.0, 3.0)

  # keep track of input values for later analysis
  real.vals <- c(a0, a1, a4, b0, b1, b2, b3, th)

  # set jags run parameters
  nc <- 3; na<- 500; nb <- 100; ni <- 1000; nt <- 1

  # parameters to monitor
  params <- c("a0", "a1", "a4", "b0", "b1","b2",
              "b3", "th")

  # make data
  sim.data <- sim.nmix(n.sites = n.sites, n.surveys = n.surveys,
                       n.years = n.years,
                       b0 = b0, b1 = b1, b2 = b2,
                       b3 = b3,
                       a0 = a0, a1 = a1, a4 = a4,
                       th = th)

  # bundle jags with site-survey-year x4
  jags.data.big <- list(Y = sim.data$Y, x1 = sim.data$x1,
    x2 = sim.data$x2, x3 = sim.data$x3,
    x1.p = sim.data$x1.p, x4 = sim.data$x4,
    n.sites = sim.data$n.sites, n.surveys = sim.data$n.surveys,
    n.years = sim.data$n.years)

  # initial values for first jags run
  N.init <- sim.data$Y # initial count values
  N.init[is.na(N.init)] <- 1 # clean up NA's
  N.init <- apply(N.init, c(1, 3), max) + 1 # zero values cause trouble
  inits <- function() list(N = N.init, b0 = rnorm(1, 0, 0.01),
                           a0 = rnorm(1, 0, 0.01), a1 = rnorm(1, 0, 0.01),
                           b2 = rnorm(1,0,0.01), b1 = rnorm(1, 0, 0.01),
                           a4 = rnorm(1, 0, 0.01), b3 = rnorm(1, 0, 0.01),
                           th = runif(1, 0.5, 2.5))

  # new jags model with expanded x4
  jags.model.string.expand <- "
  model {
  # priors
  a0 ~ dnorm(0, 0.01)       # detection intercept
  a1 ~ dnorm(0, 0.01)     # detection cov 1 effect
  a4 ~ dnorm(0, 0.01)     # detection cov 3 effect
  b0 ~ dnorm(0, 0.01)     # lambda intercept
  b1 ~ dnorm(0, 0.01)   # lambda cov 1 effect
  b2 ~ dnorm(0, 0.01)   # lambda cov 2 effect
  b3 ~ dnorm(0, 0.01)        # year effect
  th ~ dunif(0, 5)    # overdispersion size
  # abundance component
  for (k in 1:n.years){
  for (i in 1:n.sites){
  N[i, k] ~ dnegbin(prob[i, k], th) # negative binomial
  prob[i, k] <- th / (th + lambda[i, k]) # overdispersion
  log(lambda[i, k]) <- b0 + (b1 * x1[i]) +
  (b2 * x2[i]) + (b3 * x3[k])
  # detection component, note new i, j, k, subscript for x4
  for (j in 1:n.surveys){
  Y[i, j, k] ~ dbin(p[i,j,k], N[i,k])
  p[i, j, k] <- exp(lp[i,j,k]) / (1 + exp(lp[i,j,k]))
  lp[i, j, k] <- a0 + (a1 * x1.p[i]) + (a4 * x4[i, j, k])
  } # close j loop
  } # close i loop
  } # close k loop
  } # close model loop
  "
  # call jags for site-survey-year x4
  out.jags.big <- run.jags(model=jags.model.string.expand, data=jags.data.big,
                           monitor=params, n.chains=nc, inits=inits,
                           burnin=nb, adapt=na, sample=ni,
                           thin=nt, modules="glm on", method="parallel")

  # collect results
  bo <- as.data.frame(round(summary(out.jags.big),2))[,c(1:5,9,11)]
  bo$simtype <- "big"; bo$a4 <- a4; bo$par <- row.names(bo)
  bo$true.vals <- real.vals

  # bundle jags data with site-year x4
  jags.data.small <- list(Y = sim.data$Y, x1 = sim.data$x1,
    x2 = sim.data$x2, x3 = sim.data$x3,
    x1.p = sim.data$x1.p, x4 = sim.data$x4.m.arr,
    n.sites = sim.data$n.sites, n.surveys = sim.data$n.surveys,
    n.years = sim.data$n.years)

  # call jags for site-year x4
  out.jags.small <- run.jags(model=jags.model.string.expand,
                             data=jags.data.small,
                             monitor=params, n.chains=nc, inits=inits,
                             burnin=nb, adapt=na, sample=ni,
                             thin=nt, modules="glm on", method="parallel")
  so <- as.data.frame(round(summary(out.jags.small),2))[,c(1:5,9,11)]
  so$simtype <- "small"; so$a4 <- a4; so$par <- row.names(so)
  so$true.vals <- real.vals

  # combine results
  all.out <- rbind(bo,so)
  all.out$diffs <- all.out$Median - all.out$true.vals
  all.out$range <- all.out$Lower95 - all.out$Upper95
  row.names(all.out) <- NULL
  return(all.out)
}

# first round
ptm <- proc.time()
sim.out <- sim.fun()

# subsequent rounds
nsims <- 49
for(w in 2:nsims){
  iter.out <- sim.fun()
  sim.out <- rbind(sim.out, iter.out)
  cat(paste("\n\n\n this is iteration", w, "\n\n\n"))
}
proc.time() - ptm
round(simu.time.1 <- proc.time() - ptm, 2)[3]
# ##############################################################################



# explore results, example 2 ###################################################
sim.out$parf <- factor(sim.out$par, labels=c("alpha[0]",
  "alpha[1]", "alpha[4]", "beta[0]", "beta[1]", "beta[2]",
  "beta[3]", "theta"))
# plot
png("fig2.png", width = 6, height = 5, units = 'in', res = 600)
ggplot(data=sim.out, aes(x=a4, y=diffs, colour=simtype, linetype=simtype)) +
  geom_point(pch=1, size=1.3) +
  facet_wrap(~factor(parf), scales="fixed",labeller=label_parsed) +
  xlab(expression(paste("Value of ", alpha[4], " used to simulate data"))) +
  ylab("Difference between posterior median and true parameter value") +
  geom_smooth(span=3, se=F, size=0.6) +
  theme_acbs() + scale_color_manual(values=c("black","gray50")) +
  theme(strip.background = element_rect(colour = "black", fill = "white"),
    legend.position="none", axis.title.x = element_text(size=12)) +
  scale_linetype_manual(values=c(1,1))
dev.off()
# ##############################################################################



# inla analysis, example 3 #####################################################
# get real mallard data from unmarked package
data(mallard)

# format unmarked data
mallard.umf <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                   obsCovs = mallard.obs)
mallard.umf[1:6,]

# format mallard data
length <- mallard.site[ , "length"]
elev <- mallard.site[ , "elev"]
forest <- mallard.site[ , "forest"]
mean.ivel <- rowMeans(mallard.obs$ivel, na.rm = T) # average per site
mean.ivel[is.na(mean.ivel)] <- mean(mean.ivel, na.rm = T)
mean.date <- rowMeans(mallard.obs$date, na.rm = T) # average per site
mean.date.sq <- mean.date^2

# dataframe representation
mallard.inla.df <- data.frame(y1 = mallard.y[ , 1], y2 = mallard.y[ , 2],
  y3 = mallard.y[ , 3], length, elev, forest, mean.ivel, mean.date,
  mean.date.sq)
round(head(mallard.inla.df), 3)
str(mallard.inla.df)

# make inla.mdata object
counts.and.count.covs <- inla.mdata(mallard.y, 1, length, elev, forest)

# run inla model
out.inla.2 <- inla(counts.and.count.covs ~ 1 + mean.ivel +
mean.date + mean.date.sq,
  data = list(counts.and.count.covs = counts.and.count.covs,
    mean.ivel = mean.ivel, mean.date = mean.date,
    mean.date.sq = mean.date.sq),
  family = "nmixnb",
  control.fixed = list(mean = 0, mean.intercept = 0, prec = 0.01,
    prec.intercept = 0.01),
  control.family = list(hyper = list(theta1 = list(param = c(0, 0.01)),
    theta2 = list(param = c(0, 0.01)), theta3 = list(param = c(0, 0.01)),
    theta4 = list(param = c(0, 0.01)), theta5 = list(prior = "flat",
    param = numeric()))))
summary(out.inla.2, digits = 3)
# ##############################################################################



# unmarked analysis, example 3 #################################################
# run unmarked model
out.unmk.2 <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest,
                     mixture = "NB", mallard.umf)
summary(out.unmk.2)
# ##############################################################################



# explore results, example 3 ###################################################
# get unmk out
e1 <- as.matrix(out.unmk.2@estimates@estimates$det@estimates)
ci1 <- as.matrix(confint(out.unmk.2, type="det"))
e2 <- as.matrix(out.unmk.2@estimates@estimates$state@estimates)
ci2 <- as.matrix(confint(out.unmk.2, type="state"))
e3 <- exp(as.matrix(out.unmk.2@estimates@estimates$alpha@estimates))
ci3 <- exp(as.matrix(confint(out.unmk.2, type="alpha")))
e4 <- rbind(e1,e2,e3)
ci4 <- rbind(ci1,ci2,ci3)
unmk.out.tab.2 <- cbind(e4,ci4); colnames(unmk.out.tab.2) <- c("X1", "X2", "X3")

# get inla out
in1 <- out.inla.2$summary.fixed[,c(4, 3, 5)]
in2 <- out.inla.2$summary.hyperpar[,c(4, 3, 5)][-5,]
in3 <- 1/out.inla.2$summary.hyperpar[,c(4, 3, 5)][5,]
in3 <- in3[,c(1,3,2)]
inla.out.tab.2 <- rbind(in1,in2,in3);
colnames(inla.out.tab.2) <- c("X1", "X2", "X3")

# combine out
all.out.tab.2 <- round(rbind(unmk.out.tab.2, inla.out.tab.2),2)
colnames(all.out.tab.2) <- c("Estimate","Lower", "Upper")
all.out.tab.2$Parameter <- factor(rep(c("logit(p) Intercept",
                                        "Survey intensity",
                                        "Survey date", "Date squared",
                                        "log(lambda) Intercept",
                                        "Transect length", "Transect elevation",
                                        "Forest cover",
                                        "Overdispersion"), 2))
all.out.tab.2$Technique <- factor(rep(c("unmarked", "R-INLA"), each=9))
rownames(all.out.tab.2) <- NULL
pd <- position_dodge(0.2)

# plot
png("fig3.png", width = 6, height = 3.5, units = 'in', res = 600)
ggplot(data=all.out.tab.2, aes(y=Estimate, x=Parameter, col=Technique)) +
  geom_point(position=pd) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position=pd) +
  theme_acbs() + scale_color_manual(values=c("black","gray50")) +
  coord_flip() + geom_hline(yintercept=0, col="gray30", lty=2) +
  ylab("Estimated parameter value (95% CI or CrI)") +
  xlab("Model variable")
dev.off()
################################################################################


