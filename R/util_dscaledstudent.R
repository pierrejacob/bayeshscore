#'@rdname dtscaled
#'@title dtscaled
#'@description dtscaled
#'@export
#------ Define density of scaled Student t distribution -----------#
dtscaled = function(y,df,s2,log){
  if (log==TRUE) {
    return (lgamma((df+1)/2)-lgamma(df/2)-(1/2)*log(df*pi*s2)-((df+1)/2)*log(1+(1/df)*(y)^2/s2))
  }
  else {
    return (exp(lgamma((df+1)/2)-lgamma(df/2)-(1/2)*log(df*pi*s2)-((df+1)/2)*log(1+(1/df)*(y)^2/s2)))
  }
}
# Sanity check: density of scaled t distribution
# df = 5
# s2 = 2
# g11 <- ggplot(data.frame(sim=rtscaled(10000,df,s2)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dtscaled(y,df,s2,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)
