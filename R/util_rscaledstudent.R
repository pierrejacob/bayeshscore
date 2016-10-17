#'@rdname rtscaled
#'@title rtscaled
#'@description rtscaled
#'@export

#------ Define sampler of scaled Student t distribution -----------#
rtscaled = function(N,df,s2){
  return (sqrt(s2)*rt(N,df))
}
# Sanity check: density of scaled t distribution
# df = 5
# s2 = 2
# g11 <- ggplot(data.frame(sim=rtscaled(10000,df,s2)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dtscaled(y,df,s2,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)
