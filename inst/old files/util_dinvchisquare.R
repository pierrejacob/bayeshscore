#'@rdname dinvchisq
#'@title dinvchisq
#'@description dinvchisq
#'@export

#------ Define density of scaled inverse chi-square distribution -----------#
dinvchisq = function(y,df,s2=1,log){
  return (dinvgamma(y,a = df/2,b = df*s2/2,log))
}
# # Sanity check: density of scaled inverse chi-square distribution
# df = 20
# s2 = 2
# g11 <- ggplot(data.frame(sim=rinvchisq(10000,df,s2)), aes(x = sim)) + geom_density() +
#   stat_function(fun = function(y)dinvchisq(y,df,s2,FALSE), colour = "red",size=1.5,linetype=2)
# grid.arrange(g11,ncol = 1, nrow = 1)
