#'@rdname dinvgamma
#'@title dinvgamma
#'@description dinvgamma
#'@export

#------ Define density of inverse gamma distribution -----------#
dinvgamma = function(y,a,b,log){
  if (log==TRUE) {
    return (a*log(b)-lgamma(a)-(a+1)*log(y)-b/y)
  }
  else {
    return (exp(a*log(b)-lgamma(a)-(a+1)*log(y)-b/y))
  }
}
# # Sanity check: density of inverse Gamma distribution
# a = 10
# b = 1
# g11 <- ggplot(data.frame(sim=rinvgamma(10000,a,b)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dinvgamma(y,a,b,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)
