Rd <- file.path("~/GitHub/inlaNMix/meehan_et_al_inla_nmix_jss/inla.nmix.lambda.fitted.Rd")

outfile <- file.path("~/GitHub/inlaNMix/meehan_et_al_inla_nmix_jss/inla.nmix.lambda.fitted.html")

tools::Rd2HTML(Rd, out = outfile, package = "INLA")

