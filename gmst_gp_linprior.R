# gmst_gp.R

# Upload global mean surface temperature (gmst)
hadcrut4file = 'https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/time_series/HadCRUT.4.6.0.0.annual_ns_avg.txt'
hadcrut4_all = read.table(hadcrut4file)
hadcrut4 = hadcrut4_all[ , 1:2]
colnames(hadcrut4) = c('Year', 'Anomaly')

# ---------------------------------------------------------
# Model gmst as a Gaussian process
# ---------------------------------------------------------


# take a sample of HadCRUT4 data, so that fitting isn't so slow for the minute.
# call the sample y
 ix = seq(from = 1, to = nrow(hadcrut4), by = 10)
 y = hadcrut4[ix, 2]

# Use full data set
#y = hadcrut4[,2]
T = length(y)
t = 1:T


library(R2jags)
library(MASS) # Useful for mvrnorm function

model_code = '
model
{
  # Likelihood
  y ~ dmnorm(Mu, Sigma.inv)
  Sigma.inv <- inverse(Sigma)

  # Set up mean and covariance matrix
  for(i in 1:T) {
  #Mu[i] <- alpha
  Mu[i] <- alpha + beta * t[i]
  Sigma[i,i] <- pow(sigma, 2) + pow(tau, 2)

  for(j in (i+1):T) {
  Sigma[i,j] <- pow(tau, 2) * exp( - rho * pow(t[i] - t[j], 2) )
  Sigma[j,i] <- Sigma[i,j]
  }
  }
  
  # JAGS uses precision (= 1/variance) not standard deviation, so
  # dnorm (0,0.01) is equiv to dnorm(0, 10^2) in normal R

  alpha ~ dnorm(0, 0.01) # default dnorm(0, 0.01)
  sigma ~ dunif(0, 10) # default dunif(0,10)
  tau ~ dunif(0, 10) # default dunif(0, 10)
  rho ~ dunif(0.01, 5) # default dunif(0.1, 5)
  beta ~ dnorm(0, 0.1)

}
  '

model_data = list(T = T, y = y, t = t)

# Choose the parameters to watch
model_parameters =  c("alpha", "sigma", "tau", "rho", "beta")

# Run the model - This can be slow with lots of data
model_run = jags(data = model_data,
                 parameters.to.save = model_parameters,
                 model.file=textConnection(model_code),
                 n.chains=4, # Number of different starting positions
                 n.iter=1000, # Number of iterations
                 n.burnin=200, # Number of iterations to remove at start
                 n.thin=2) # Amount of thinning

alpha = model_run$BUGSoutput$sims.list$alpha
tau = model_run$BUGSoutput$sims.list$tau
sigma = model_run$BUGSoutput$sims.list$sigma
rho = model_run$BUGSoutput$sims.list$rho
beta = model_run$BUGSoutput$sims.list$beta

pdf(file = 'parameter_estimates_linprior.pdf', width = 7, height = 5)
par(mfrow = c(2,3))
hist(alpha, breaks=30)
hist(tau, breaks=30)
hist(sigma, breaks=30)
hist(rho, breaks=30)
hist(beta, breaks=30)
dev.off()

# Plot the GP

T_new = 100 # number of interpolation points for looking at the way the GP mean
            # function fits. 100 seems a good number, which captures the variability

t_new = seq(0,T,length=T_new)
Mu = rep(mean(alpha), T) + (t* mean(beta))
Mu_new = rep(mean(alpha), T_new) + (t_new * mean(beta))
Sigma_new = mean(tau)^2 * exp( -mean(rho) * outer(t, t_new, '-')^2 )
Sigma_star = mean(sigma)^2 * diag(T_new) + mean(tau)^2 * exp( - mean(rho) * outer(t_new,t_new,'-')^2 )
Sigma = mean(sigma)^2 * diag(T) + mean(tau)^2 * exp( - mean(rho) * outer(t,t,'-')^2 )

# Use fancy equation to get predictions
pred_mean = Mu_new + t(Sigma_new)%*%solve(Sigma, y - Mu)
pred_var = Sigma_star - t(Sigma_new)%*%solve(Sigma, Sigma_new)

pdf(file = 'interpolated_gp_linprior.pdf', width = 7, height = 5)
plot(t,y)
lines(t_new, pred_mean, col='red')

pred_low = pred_mean - 1.95 * sqrt(diag(pred_var))
pred_high = pred_mean + 1.95 * sqrt(diag(pred_var))
lines(t_new, pred_low, col = 'red', lty = 2)
lines(t_new, pred_high, col = 'red', lty = 2)
dev.off()

# Extrapolate
T_ext = 30
t_ext = seq(T+1,T*1.5,length = T_ext) # this needs to be at same res as data
Mu = rep(mean(alpha), T) + (t* mean(beta))
Mu_ext = rep(mean(alpha), T_ext) + (t_ext * mean(beta))
Sigma_ext = mean(tau)^2 * exp( -mean(rho) * outer(t, t_ext, '-')^2 )
Sigma_ext_star = mean(sigma)^2 * diag(T_ext) + mean(tau)^2 * exp( - mean(rho) * outer(t_ext,t_ext,'-')^2 )
Sigma = mean(sigma)^2 * diag(T) + mean(tau)^2 * exp( - mean(rho) * outer(t,t,'-')^2 )

# Use fancy equation to get predictions
ext_mean = Mu_ext + t(Sigma_ext)%*%solve(Sigma, y - Mu)
ext_var = Sigma_ext_star - t(Sigma_ext)%*%solve(Sigma, Sigma_ext)
ext_low = ext_mean - 1.95 * sqrt(diag(ext_var))
ext_high = ext_mean + 1.95 * sqrt(diag(ext_var))

pdf(file = 'extrapolate_mean_linprior.pdf', width = 7, height = 5)
plot(t,y, xlim = c(0,T*1.5), ylim = range(y, ext_high, ext_low))

# plot the interpolated best estimate and uncertainty
lines(t_new, pred_mean, col = 'red', lty = 1)
lines(t_new, pred_low, col = 'red', lty = 2)
lines(t_new, pred_high, col = 'red', lty = 2)

# plot the extrapolated best estimate and uncertainty
lines(t_ext, ext_mean, col='blue')
lines(t_ext, ext_low, col = 'blue', lty = 2)
lines(t_ext, ext_high, col = 'blue', lty = 2)
dev.off()


# recreate the plot from the previous code chunk
pdf(file = 'sample_uncertainty_linprior.pdf', width = 7, height = 5)
plot(t,y, xlim = c(0,T*1.5), ylim = range(y, ext_high, ext_low))

# Take samples from the markov chain to show the possible solutions.
for(i in 1500:1600){
  Mu = rep(mean(alpha[i]), T) + (t* mean(beta[i]))
  Mu_new = rep(mean(alpha[i]), T_new) + (t_new * mean(beta[i]))
  Sigma_new = tau[i]^2 * exp( -rho[i] * outer(t, t_new, '-')^2 )
  Sigma_star = sigma[i]^2 * diag(T_new) + tau[i]^2 * exp( - rho[i] * outer(t_new,t_new,'-')^2 )
  Sigma = sigma[i]^2 * diag(T) + tau[i]^2 * exp( - rho[i] * outer(t,t,'-')^2 )
  
  # Use fancy equation to get predictions
  pred_samp = Mu_new + t(Sigma_new)%*%solve(Sigma, y - Mu)
  lines(t_new, pred_samp, col='grey')
}

for(i in 1:200){
  T_ext = 30
  t_ext = seq(T+1,T*1.5,length=T_ext)
  t_ext = seq(T+1,T*1.5,length = T_ext) # this needs to be at same res as data
  Mu = rep(mean(alpha[i]), T) + (t* mean(beta[i]))
  Mu_ext = rep(mean(alpha[i]), T_ext) + (t_ext * mean(beta[i]))
  Sigma_ext = tau[i]^2 * exp( -rho[i] * outer(t, t_ext, '-')^2 )
  Sigma_ext_star = sigma[i]^2 * diag(T_ext) + tau[i]^2 * exp( - rho[i] * outer(t_ext,t_ext,'-')^2 )
  Sigma = sigma[i]^2 * diag(T) + tau[i]^2 * exp( - rho[i] * outer(t,t,'-')^2 )
  
  # Use fancy equation to get predictions
  ext_mean = Mu_ext + t(Sigma_ext)%*%solve(Sigma, y - Mu)
  lines(t_ext, ext_mean, col='blue')
}

# plot the interpolated best estimate and uncertainty
lines(t_new, pred_mean, col = 'red', lty = 1)
lines(t_new, pred_low, col = 'red', lty = 2)
lines(t_new, pred_high, col = 'red', lty = 2)

# plot the extrapolated best estimate and uncertainty
lines(t_ext, ext_mean, col='blue')
lines(t_ext, ext_low, col = 'blue', lty = 2)
lines(t_ext, ext_high, col = 'blue', lty = 2)
lines(t_new, pred_samp, col='grey')
dev.off()



  
  