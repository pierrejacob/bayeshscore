##################################################################################################
# This implements example 3.2. in O'Hagan (1995) with analytical formulas
# using conjugacy of the chosen priors.
##################################################################################################
set.seed(19)
# set hyperparameters
muprior = 0
sigma2prior = 100
nu0 = 1
s02 = 1
#--------------------------------------------------------------------------------------------
# Generate data
nobservations = 100
observations = list()
DGP = list(c(0,5), # Case 1: DGP = N(0,5), Model 2 is correct
           c(1,1), # Case 2: DGP = N(1,1), Model 1 is correct
           c(2,3), # Case 3: DGP = N(2,3), both model 1 and 2 are misspecified
           c(0,1)) # Case 4: DGP = N(0,1), both model 1 and 2 are correct
for (i in 1:4) {
  observations[[i]] = matrix(rnorm(nobservations, DGP[[i]][1], sqrt(DGP[[i]][2])), nrow = 1)
}

#--------------------------------------------------------------------------------------------
# Sanity check (exact computation)
#--------------------------------------------------------------------------------------------
plot_exact = list()
results_exact = data.frame()
for (i in 1:4){
  mu_star = DGP[[i]][1]
  sigma2_star = DGP[[i]][2]
  expected_H_model1_exact = rep(0,4)
  expected_H_model2_exact = rep(0,4)
  sigma2_post = rep(NA, nobservations)
  mu_post = rep(NA, nobservations)
  nu_post = rep(NA, nobservations)
  s2_post = rep(NA, nobservations)
  hexact1 = rep(NA, nobservations)
  hexact2 = rep(NA, nobservations)
  for (t in 1:nobservations) {
    sigma2_post[t] = 1/(t + 1/sigma2prior)
    mu_post[t] = (sum(observations[[i]][,1:t]) + (1/sigma2prior)*muprior)*sigma2_post[t]
    nu_post[t] = nu0 + t
    s2_post[t] = (nu0*s02 + sum(observations[[i]][,1:t]^2))/nu_post[t]
    hexact1[t] = (observations[[i]][,t]^2 - 2*observations[[i]][,t]*mu_post[t] + sigma2_post[t] + mu_post[t]^2 - 2) + sigma2_post[t]
    hexact2[t] = observations[[i]][,t]^2*(2/(nu_post[t]*s2_post[t]^2) + 1/s2_post[t]^2) - 2/s2_post[t] + (2/(nu_post[t]*s2_post[t]))*(observations[[i]][,t])^2
  }
  hexact1 = cumsum(hexact1)
  hexact2 = cumsum(hexact2)
  # expected_H_model1_exact[i] = (sigma2_star-2)
  # expected_H_model2_exact[i] = (1/sigma2_star^2)*(sigma2_star+mu_star^2-2*sigma2_star)
  expected_H_model1_exact[i] = (sigma2_star+mu_star^2-2*mu_star*mu_post[nobservations]+mu_post[nobservations]^2-2)
  thetas2_postmean = nu_post[nobservations]*s2_post[nobservations]/(nu_post[nobservations]-2)
  expected_H_model2_exact[i] = (1/thetas2_postmean^2)*(sigma2_star+mu_star^2-2*thetas2_postmean)
  results_exact = rbind(results_exact, data.frame(time = 1:nobservations,
                                                  case = factor(i),
                                                  hscore1 = hexact1,
                                                  hscore2 = hexact2,
                                                  hfactor = hexact2 - hexact1,
                                                  Ehscore1 = (1:nobservations)*expected_H_model1_exact[i],
                                                  Ehscore2 = (1:nobservations)*expected_H_model2_exact[i],
                                                  Ehfactor = (1:nobservations)*(expected_H_model2_exact[i]-expected_H_model1_exact[i])))
  local({
    i = i;
    results_exact = results_exact;
    hexact1 = hexact1;
    hexact2 = hexact2;
    expected_H_model1_exact = expected_H_model1_exact;
    expected_H_model2_exact = expected_H_model2_exact;
    results_all = results_all;
    index = 1:nobservations
    # plot exact H-scores
    plot_exact[[i]] <<-ggplot() +
      geom_line(aes(index,expected_H_model1_exact[i]),col="red",linetype="dashed") +
      geom_line(aes(index,expected_H_model2_exact[i]),col="blue",linetype="dashed") +
      geom_line(aes(index,hexact1/index),col="red") +
      geom_line(aes(index,hexact2/index),col="blue") +
      # geom_point(aes(index,hexact1/index),col="red") +
      # geom_point(aes(index,hexact2/index),col="blue") +
      xlab("time index") + ylab("H-score / time")
  })
}
do.call(grid.arrange,c(plot_exact, ncol = 2))

# plot exact H-factor
g1 = ggplot(subset(results_exact,case==1)) +
  geom_line(aes(time,hfactor,col=case),size=1) +
  geom_line(aes(time,Ehfactor,col=case),linetype="dashed") +
  xlab("time index") + ylab("H-factor")
g2 = ggplot(subset(results_exact,case==2)) +
  geom_line(aes(time,hfactor,col=case),size=1) +
  geom_line(aes(time,Ehfactor,col=case),linetype="dashed") +
  xlab("time index") + ylab("H-factor")
g3 = ggplot(subset(results_exact,case==3)) +
  geom_line(aes(time,hfactor,col=case),size=1) +
  geom_line(aes(time,Ehfactor,col=case),linetype="dashed") +
  xlab("time index") + ylab("H-factor")
g4 = ggplot(subset(results_exact,case==4)) +
  geom_line(aes(time,hfactor,col=case),size=1) +
  geom_line(aes(time,Ehfactor,col=case),linetype="dashed") +
  xlab("time index") + ylab("H-factor")
grid.arrange(g1,g2,g3,g4,ncol=2)
