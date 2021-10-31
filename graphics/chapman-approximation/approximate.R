library(plotly)
R = 10
# this plot zooms in on behavior that describes sunsets, which is very sensitive and approximates Ch∥*erf(x)
# domain = seq(R-10, R+3, 0.1) 
# this plot zooms in on behavior that describes noontime scattering, which approximates exp(-x)
# domain = seq(R-10, R+3, 0.1) 
# this plot views the entire appreciable domain of the problem
xdomain = seq(-(R+4), (R+4), 0.3)
ydomain = seq(-(R+4), (R+4), length=length(xdomain))
x = rep(xdomain, length(ydomain))
y = rep(ydomain, each=length(xdomain))
a = 0.6
# `Ch` is the modified "Chapman" function, the function by which we nudge a naive integration-by-parts, 
# if z==r, then Ch=Ch∥ where "Ch∥" is the Chapman function that is described by Christian Schuler (2018)
Ch = function(x,r) { (1+1/(2*r)) * sqrt(pi*r/2) + a*x } 
# `Ixgt0` is a simplified version of integral that is only applicable when x>0
Ixgt0 = function(x,r,R) { pmin(exp(R-r),1)/(x/r + 1/Ch(x,r)) }
# `I` is 
I = function(x,y,R) { sign(x)*(Ixgt0(sqrt(pmax(R*R-y*y,0)),y,R) - Ixgt0(abs(x),sqrt(pmax(x*x+y*y,0)),R)) }
# long-hand form, with different handling for values outside sensible range of input
I = function(x,y,R) { sign(x)*((exp(R-y))/(sqrt(pmax(R*R-y*y,0))/y + 1/Ch(sqrt(pmax(R*R-y*y,0)),y)) - (exp(R-sqrt(pmax(x*x+y*y,0))))/(abs(x)/sqrt(pmax(x*x+y*y,0)) + 1/Ch(abs(x),sqrt(pmax(x*x+y*y,0))))) }
# this plot makes it easier to see the outline of the earth and exponential dropoff:
plot_ly(x=x,y=y,z=Ixgt0(abs(x),sqrt(x*x+y*y),R),type='scatter3d') 
# plot_ly(x=x,y=y,z=I(x,y,R),type='scatter3d')