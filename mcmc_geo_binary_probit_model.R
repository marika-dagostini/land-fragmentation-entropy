library(PrevMap)
library(geoR)
# number of observations
n = nrow(data)
# number of regression coefficients
n.mun = length(unique(data$MUNICIPALITY))
# priors specification
control.prior = control.prior(
  beta.mean = rep(0, n.mun),
  beta.covar = diag(1, n.mun, n.mun),
  uniform.phi = c(0, 10),
  uniform.sigma2 = c(0, 100)
)
# tuning parameter used in the MCMC algorithm
control.mcmc = control.mcmc.Bayes(
  n.sim = 10,
  burnin = 0,
  thin = 1,
  h.theta1 = 0.05,
  h.theta2 = 0.05,
  L.S.lim = c(1, 50),
  epsilon.S.lim = c(0.01, 0.02),
  start.beta = rep(0, n.mun),
  start.sigma2 = 1,
  start.phi = 0.15,
  start.S = rep(0, n),
  start.nugget = NULL,
  binary = TRUE
)
# two - levels geostatistical binary probit model
fit.Bayes = binary.probit.Bayes(
  formula = frag ~ 1 + factor(COMUNE),
  coords = ~ x.pxl + y.pxl,
  data = data,
  ID.coords = ID.coords,
  control.prior = control.prior,
  control.mcmc = control.mcmc,
  kappa = 5 / 2
)