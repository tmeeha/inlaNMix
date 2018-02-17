


## inla.nmix.lambda.fitted() function ##########################################
##
## The inla.nmix.lambda.fitted() function is a helper function for taking an
## 'nmix' or 'nmixnb' model, and computing expected abundance values for each
## site or site-by-year combination in the data, using the linear predictor for ## lambda (hereafter, 'fitted values').
##
## The uncertainty associated with fitted values derives from repeated sampling
## of INLA posteriors for the parameters of the linear predictor, and repeated
## solving of the linear predictor equation for each site or site-by-year
## combination. By default, fitted values are exponentiated to the count
## scale, and summaries of approximate fitted posteriors are returned by the
## function. Full estimated posteriors of fitted values are also available as
## an optional output.
##
## Tim Meehan, 16 February 2018. ###############################################



# Begin function ###############################################################
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

  # Discard overdispersion marginal posterior if 'nmixnb'
  if(fam == "nmixnb"){
    hyperpar.samples <- hyperpar.samples[,-(ncol(hyperpar.samples))]
  }
  n.samp.covs <- ncol(hyperpar.samples)
  if(n.lambda.covs != n.samp.covs) {
    stop("The number of hyperparameters and covariates does not match")
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

  } # End function ###############################################################


