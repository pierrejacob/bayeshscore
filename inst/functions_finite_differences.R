#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# This file contains some generic functions to compute derivatives
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


# This function computes the partial first derivative of the log of a function g 
# with respect to the k-th variable using forward difference
d1log = function(g,y,k = 1,h = 10^(-8)) {
  d = length(y)
  e = vector("numeric",length = d)
  e[k] = h
  gright = g(y+e)
  gcenter = g(y)
  if (gright==0) {
    if (gcenter==0) {
      return (0)
    }
    else {
      return (((gright-gcenter)/h)/gcenter)
    }
  }
  else {
    return ((log(gright)-log(gcenter))/h)
  }
}

# This function computes the partial second derivative of the log of a function g
# with respect to the k-th variable
d2log = function(g,y,k = 1, h = 10^(-4)) {
  d = length(y)
  e = vector("numeric",length = d)
  e[k] = h
  gright = g(y+e)
  gcenter = g(y)
  gleft = g(y-e)
  if (gcenter!=0) {
    return (-(((gright-gcenter)/h)/gcenter)^2+(1/gcenter)*(gright-2*gcenter+gleft)/(h^2))
  }
  else {
    if ((gright==0)&(gleft==0)) {
      return (0)
    }
    else {
      return (0)
    }
  }
}