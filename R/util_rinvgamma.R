#'@rdname rinvgamma
#'@title rinvgamma
#'@description rinvgamma
#'@export

#------ Define sampler of inverse gamma distribution -----------#
rinvgamma = function(N,a,b){
  return (1/rgamma(N,a,rate = b))
}
# # Sanity check: density of inverse Gamma distribution
# a = 10
# b = 1
# g11 <- ggplot(data.frame(sim=rinvgamma(10000,a,b)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dinvgamma(y,a,b,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)
