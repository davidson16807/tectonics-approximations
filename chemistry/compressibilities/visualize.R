# import libraries
library(plotly)

"
'visualize.R' documents our attempt to approximate charts for the 'compressibility' factor, `Z`
Z is the factor by which the ideal gas law fails to preduct the equation of state for real gases. 

Our model has the aspirational goal to predict behavior in exotic conditions, so our estimate for Z must be valid for a wide range of p and T,
and although we currently do not plan to run it repeatedly for every grid cell in a raster, we do not want to preclude the possibility, 
so performance is a secondary objective. Ideally, we would also like the approximation to be easy to implement, so we want to avoid interpolation tables.

Tackacs (1989) gives an excellent review of past approximations for Z but most either require a sacrafice to one of our goals.
Fast approximations either use interpolation tables or simple functions that only perform well over small regions.
In any case, we also tend to be megalomaniacal and would just like to see if we can do it ourselves. Therefore, we attempt our own approximation.

We use two data sets to represent known values for Z on a computer.
The first was created by digitizing the well known charts from Nelson and Obert (1955). The script we used to digitize them can be found under 'digitize.R'.
The second was created using OCR on tables from a more extensive, more recent, but less well known paper by Varsanyi (1986).
The generation of the Varsanyi dataset from the pdf is documented within 'varsanyi1986.R'.
Varsanyi's dataset appears ideal in that it provides the actual tables that underly his charts, 
there are over 5000 data points, and the values are reported to 3 decimal places,
but it is also valuable in that he provides estimates for different compounds that vary by their critical compressibility, 'Zc'.
The estimates we use from Varsanyi were taken by observing a sample of air, which has a Zc of 278.
We chose estimates for this compound rather than others since air is very representative of what we will most often try to model.
"
nelson = read.csv('nelson1955.csv')
varsanyi = read.csv('varsanyi1986.csv')
pTZ = rbind(varsanyi[varsanyi$Zc==0.278 & is.finite(varsanyi$Z), c('p','T','Z')],  nelson[nelson$p>25 | nelson$T>10 | (0.6<nelson$T&nelson$T<0.8),c('p','T','Z')])

# fundamental input
p = pTZ$p
T = pTZ$T
Z = pTZ$Z

# shorthands
V = p/T  # this is representative but not equal to volume, since volume is proportional to p/T
T1 = 1/T # shorthand

# color indicators for graphs
hot=Z*0+3
warm=Z*0+2
cool=Z*0+1
cold=Z*0+0

# basic graphs
plot_ly(x=p, y=T,  z=Z) #naive first plot

"
After doing some unstructured investigation, I discovered a planar relation between p/T, 1/T, and Z.
This makes it much easier to find an adequate approximation for most of the graph by linear regression.
"
plot_ly(x=V, y=T1, z=Z) 
test = p>3; summary(lm(Z[test]~V[test]+T1[test])) # linear regression on the planar relationship
L = 1.12 + 0.101*V - 0.928*T1 # "L" is the result of this linear regression
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L), color = c(cold, hot)) # plot the regression

"
From the same graph of V and T1 we also see what looks like a sigmoid relationship between T1 and the minimum value of Z for a given V.
We fit a sigmoidal function to it manually, but if needed we may be able to run an optimization algorithm on it.
"
S = 1/(1+exp(8.5*(T1-0.84))) # "S" is our handmade approximation to the sigmoidal relationship
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,S), color = c(cold, hot)) # plot the sigmoid

"
When L<S, we know that S no long approximates Z and we need to interpolate to a better approximation.
Fortunately, We know that Z→1 as V→0, so the other function we must find is merely the constant 1, and Z is the interpolation between L and 1, 
i.e. there is an interpolant, either `I1` or `I2` where `Z = I1*L + (1-I1)*1` or `Z = (1-I2)*L + I2*1`
"
plot_ly(x=c(V,V,V), y=c(T1,T1,T1), z=c(Z,L,S), color = c(cold, cool, hot)) # plot where S and L intersect

"
You will notice from the above graph that the purple line starts its dip back up to 1 around where the yellow line intersects the teal.
There is no magic here: we have defined the yellow line (S) to be the minimum, and the teal line (L) simply traces one of the two sides of the curve.
We can presume in these circumstances that there is always an intersection between L and S for any value T1 (though it may be nonphysical, such as where p<0).
We can find the value for V at this intersection algebraically, and then derive the corresponding Z.
By these means, we find we are always guaranteed to know the precise coordinates of the intersection for our given model of L and S.
We then have a good sense for just how fast the line ought to return to 1 at V=0, and we may even be able to derive a function for it.

L = la + b*Vi - c*T1
L = S
la + b*Vi - c*T1 = S
Vi = (S + c*T1 - la) / b

where `Vi` is the V at our intersection. 

So we now know we must start interpolation somewhere around where V=Vi and end pretty much right at where V=0.
We tried experimenting with exponentials for our interpolant, but found a simple linearstep() function gave better results.
We even tried smoothstep() but that also gave worse results, causing an unwanted decrease in dZ/dV at p=0
This is fine anyway, since linearstep() should be faster than an exponential.
"
Vi = (1/(1+exp(8.5*(T1-0.84)))  + 1.0351556*T1 - 1.1589765) / 0.1035014 # NOTE: the values here are updated from running optimization on L
Zi = 1.1589765 + 0.1035014*Vi - 1.0351556*T1
plot_ly(x=c(V,Vi), y=c(T1,T1), z=c(Z,Zi), color = c(cold, hot)) # plot where S and L intersect

# I2 = exp(-1*V/pmax(0, Vi))
# I2 = pmin(1, pmax(0, 1-1*V/pmax(0, Vi)))
I2 =  pmin(1, pmax(0,  1-1*V/pmax(0, Vi) )) 
I2 * I2 * (3.0 - 2.0 * I2); # hermite interpolation/"smoothstep"
Z2=I2+(1-I2)*L
plot_ly(x=c(V,V,Vi), y=c(T1,T1,T1), z=c(Z,Z2,Zi), color = c(I2, warm, hot)) # plot where S and L intersect

"
I then  ran optimization on our approximation for Z. It made some improvements to L, which I reincorporated in the equations above, 
and I did discover in the process that I could improve accuracy if I allowed the upper bound of the sigmoid to reach past 1.
In this sense, the sigmoid S no longer represents the minimum value taken on by Z, but rather the value at which we start interpolation.
"

Zmpe = function(params) {
	return(mean(((Z-Zvalue(params))/Z)**2))
}
Zme = function(params) {
	return(mean(((Z-Zvalue(params)))**2))
}
Zmxe = function(params) {
	return(max(abs(Z-Zvalue(params))))
}
Zmxpe = function(params) {
	return(max(abs(Z-Zvalue(params))/Z))
}
Zvalue = function(params) {
	a = params[1]; b = params[2]; c = params[3];  # L parameters
	# a = 1.1589765; b = 0.1035014; c = -1.0351556;  # L parameters
	f = params[4]; g = params[5]; h = params[6]; # S parameters
	# f = 8.5; g = 0.84; # S parameters
	L = a + b*V + c*T1
	S = h/(1+exp(f*(T1-g))) 
	Vi = (S -c*T1 - a) / b 
	I2 =  pmin(1, pmax(0,  1-1*V/pmax(0, Vi) )) 
	I2 = I2 * I2 * (3.0 - 2.0 * I2); # hermite interpolation, A.K.A. "smoothstep"
	Z2=I2+(1-I2)*L
	return(Z2)
}
Zstart = c(1.159, 0.104, -1.035, 8.5, 0.84, 1.17)
result = optim(  Zstart,  Zmpe, control=list(maxit=1e9))
print(result)
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Zvalue(result$par)), color = c(cold, hot))


plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Zvalue(Zstart)), color = c(cold, hot))


"
In trying to improve accuracy, I noticed the plane L was a major source of error, especially towards upper and lower T.
I initially thought this was because there was some warping at the edges which could be accounted for by exponents or a higher order polynomial.
However, when plotting the planar part of the graph on its own I realized it matched the same sigmoidal trend with respect to T1 that we see in S.
This is best viewed along a diagonal where Z∝-V
"
test = p>1 & T1>1; plot_ly(x=V[test], y=T1[test], z=Z[test])
"
So L could be modeled as a planar relation that's been offset by a sigmoid term. 
This sigmoid term is not to be confused with the sigmoid S, since S defined the value of Z at which we start interpolating along V,
whereas this sigmoid indicates the amount we offset along Z. 
We first try isolating the planar relation by performing linear regression on the lower portion of the sigmoid relation, around where T<1.
We could conceptually use the upper portion as well, but we lack good data for it. 
"
test = p>1 & T1>1; summary(lm(Z[test]~V[test]+T1[test]))
L0 = 0.354 + 0.149*V**0.9 -0.319*T1
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L0), color = c(cold, hot))
"
The only thing left to parameterize is the upper bound for the term (analogous to `h` above).
We find this value by taking the residual between the known values for Z and the planar relationship we just discovered.
"
test = p>3; plot_ly(x=V[test], y=T1[test], z=(Z-L0)[test]) 
"
Since we're lazy and optimization comes later, we spitball the h analog to be 0.6, the f to be 7.7, and the g to be 0.5
"
LT1 = 0.354 + -0.319*T1 + 0.6/(1+exp(7.7*(T1-0.5))) # "LT1" describes the sum of components responsible for calculating L based on T1
L = 0.149*V**0.9 + LT1
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L), color = c(cold, hot))


LT1= 0.354 + -0.319*T1 + 0.6/(1+exp(7.7*(T1-0.5)))
L  = 0.149*V**0.9 + LT1
S = 1.3/(1+exp(8.5*(T1-0.84))) # "S" is our handmade approximation to the sigmoidal relationship
Vi = ((S - LT1) / 0.149)#**(1/0.9)
I2 = pmin(1, pmax(0,  1-1*V/pmax(0, Vi) )) 
I2 = I2 * I2 * (3.0 - 2.0 * I2); # hermite interpolation, A.K.A. "smoothstep"
Z2 = I2+(1-I2)*L

plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Z2), color = c(cold, hot)) # plot where S and L intersect


Zmpe = function(params) {
	return(mean(((Z-Zvalue(params))/Z)**2))
}
Zmxe = function(params) {
	return(max(abs(Z-Zvalue(params))))
}
Zmxpe = function(params) {
	return(max(abs(Z-Zvalue(params))/Z))
}



Zvalue = function(params) {
	a = params[1]; b = params[2]; c = params[3]; d = params[4];  # L parameters
	# d = 0.9 # we allow hardcoding this since optimization tends to mess it up in nonsensical ways
	f = params[5]; g = params[6]; h = params[7]; # S parameters
	# h = 1.17 # we allow hardcoding this since optimization tends to mess it up in nonsensical ways
	i = params[8]; j = params[9]; k = params[10]; # L sigmoid term
	# h = 0.6 # we allow hardcoding this since optimization tends to mess it up in nonsensical ways
	LT1= a + c*T1 + k/(1+exp(i*(T1-j)))
	L  = b*V**d + LT1
	S  = h/(1+exp(f*(T1-g))) # "S" is our handmade approximation to the sigmoidal relationship
	Vi = ((S - LT1) / b)#**(1/d) 
	# NOTE: we comment out the exponent above for three reasons: 
	# 1.) it should not matter at the low-Z regions where this formula starts to affect output (in order words it never does anything)
	# 2.) having it around reduces performance
	# 3.) having it around threatens us with NaNs if (S-LT1)/b happens to be negative
	I2 = pmin(1, pmax(0,  1-1*V/pmax(0, Vi) )) 
	# I2 = I2 * I2 * (3.0 - 2.0 * I2); # hermite interpolation, A.K.A. "smoothstep"
	Z2 = I2+(1-I2)*L
	return(Z2)
}


"
What follows are our starting parameters and some of the best parameter sets we found
"

Zstart = c( 0.354, 0.104, -0.319, 0.9, 8.5, 0.84, 1.17, 7.7, 0.5, 0.6 )
result  = optim(  Zstart,  Zmpe, control=list(maxit=1e9))
print(result)
paste('c(', paste(round(result$par, digits=3), collapse=', '), ')')
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Zvalue(result$par)), color = c(cold, hot))

"
We arrive at the 'canonical' parameter set that we'll put in our application for gases with Zc=0.278. It's good enough.
The following are its guarantees:

* Applicable for use from 0.6<T<15 and 0<p<40.  This exceeds the range of virtually all approximations reviewed by Takacs (1989).
* Produces reasonable behavior outside the observed data range (e.g. no Taylor series explosions)
* 0.35% mean absolute percentage error, surpassing all but 3 of the approximations reviewed by Takacs (1989).
  The remaining 3 approximations that surpassed this one either made use of large memory footprints (Gray-Sims) 
  or narrower applicable ranges (Carlile-Gillett, Dranchuk-A. Kassem).
* max error not exceeding 0.26, appearing highest at regions where 2.0<T<2.5 and p>35.
* Only 2 exponentials and 1 power used, which should place runtime on the order of the fast approximations (Papay, Burnett)
* Only 10 parameters, with straight forward interpretations, making it easy to implement and manipulate without requiring reoptimization.
"
Zbest = c( 0.153, 0.145, -0.143, 0.9, 8.38, 0.833, 1.11, 8.988, 0.672, 0.773 )
Zmpe(Zbest)
Zme(Zbest)
Zbestworst = Zmxe(Zbest)
Zbestworst
T1[abs(Zvalue(Zbest)-Z) > Zbestworst - 0.01]
V[abs(Zvalue(Zbest)-Z) > Zbestworst - 0.01]
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Zvalue(Zbest)), color = c(cold, hot))



"
We now load the equivalent data set for gases with other values of Zc. Let's start with the lower bound, Zc=0.244.
"
pTZ = rbind(varsanyi[varsanyi$Zc==0.244 & is.finite(varsanyi$Z), c('p','T','Z')])

# fundamental input
p = pTZ$p
T = pTZ$T
Z = pTZ$Z
V = p/T 
T1 = 1/T 
hot=Z*0+3
warm=Z*0+2
cool=Z*0+1
cold=Z*0+0

test = p>1 & T1>1.2; plot_ly(x=(V**0.92)[test], y=T1[test], z=Z[test])
test = p>1 & T1>1.2; summary(lm(Z[test]~(V**0.92)[test]+T1[test]))
L0   = 0.139 -0.103*T1 + 0.103*V**0.92;
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L0), color = c(cold, hot))
test = p>3; plot_ly(x=V[test], y=T1[test], z=(Z-L0)[test]) 
LS   = 0.95/(1+exp(5*(T1-0.6))) # "LS" describes the contribution of T1 to L that can be described using a sigmoid. It oddly seems higher for Zc=0.244
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z-L0, LS), color = c(cold, hot))
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L0+LS), color = c(cold, hot))
LT1  = 0.139 -0.103*T1 + LS # "LT1" describes the sum of components responsible for calculating L based on T1
L = 0.103*V**0.92 + LT1
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L), color = c(cold, hot))
S  = 1/(1+exp(8*(T1-0.8)))
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,S), color = c(cold, hot))
Vi = ((S - LT1) / 0.103)#**(1/0.9)
I2 = pmin(1, pmax(0,  1-1*V/pmax(0, Vi) )) 
Z2 = I2+(1-I2)*L
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Z2), color = c(cold, warm)) # plot where S and L intersect

# C++       Zhi0   ZhiV    ZhiT1  ZhiVk  ZmidH Zmid0 Zmidmax ZhiSH ZhiS0 ZhiSmax
# R         L      L(V)    L(T1)  L(V,k) S(H)  S(0)  S(max)  LS(H) LS(0) LS(max)
Zstart = c( 0.139, 0.103, -0.103, 0.92,  8,    0.8,  1,      5,    0.6,  0.95 )
result  = optim(  Zstart,  Zme, control=list(maxit=1e9))
print(result)
paste('c(', paste(round(result$par, digits=3), collapse=', '), ')')
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Zvalue(result$par)), color = c(cold, hot))

# C++       Zhi0   ZhiV    ZhiT1  ZhiVk  ZmidH Zmid0 Zmidmax ZhiSH ZhiS0 ZhiSmax
# R         L      L(V)    L(T1)  L(V,k) S(H)  S(0)  S(max)  LS(H) LS(0) LS(max)
# Zbest = c( 0.139, 0.093, -0.082, 0.964, 4.67, 0.534, 2.283, 6.367, 0.52, 1.199 )
Zbest =   c( 0.139, 0.103, -0.103, 0.92,  8,    0.8,   1,     5,     0.6,  0.95 )
Zmpe(Zbest)
Zme(Zbest)
Zbestworst = Zme(Zbest)
Zbestworst
T1[abs(Zvalue(Zbest)-Z) > Zbestworst - 0.01]
V[abs(Zvalue(Zbest)-Z) > Zbestworst - 0.01]
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Zvalue(Zbest)), color = c(cold, hot))








"
We now load the equivalent data set for gases with other values of Zc. Let's start with the lower bound, Zc=0.316
"

pTZ = rbind(varsanyi[varsanyi$Zc==0.316 & is.finite(varsanyi$Z), c('p','T','Z')])

# fundamental input
p = pTZ$p
T = pTZ$T
Z = pTZ$Z
V = p/T 
T1 = 1/T 
hot=Z*0+3
warm=Z*0+2
cool=Z*0+1
cold=Z*0+0

plot_ly(x=V, y=T1, z=Z) 
test = p>1 & T1<0.126; plot_ly(x=(V**1.0)[test], y=T1[test], z=Z[test])
test = p>1 & T1<0.126; summary(lm(Z[test]~(V**1.0)[test]+T1[test]))
L0   = 1.016 -0.138*T1 + 0.103*V**1.0;
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L0), color = c(cold, hot))
test = p>3; plot_ly(x=V[test], y=T1[test], z=(Z-L0)[test]) 
LS   = -0.8/(1+exp(-5*(T1-0.65))) # "LS" describes the contribution of T1 to L that can be described using a sigmoid. 
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L0+LS), color = c(cold, hot))
LT1  = 1.016 -0.138*T1 + LS # "LT1" describes the sum of components responsible for calculating L based on T1
L = 0.103*V**1.0 + LT1
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,L), color = c(cold, hot))
S  = 1/(1+exp(8*(T1-0.85)))
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,S), color = c(cold, hot))
Vi = ((S - LT1) / 0.103)#**(1/0.9)
I2 = pmin(1, pmax(0,  1-1*V/pmax(0, Vi) )) 
Z2 = I2+(1-I2)*L
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Z2), color = c(cold, warm)) # plot where S and L intersect

# C++       Zhi0   ZhiV    ZhiT1  ZhiVk  ZmidH Zmid0 Zmidmax ZhiSH ZhiS0 ZhiSmax
# R         L      L(V)    L(T1)  L(V,k) S(H)  S(0)  S(max)  LS(H) LS(0) LS(max)
Zstart = c( 1.016, 0.103, -0.138, 1.0,   8,    0.85, 1,     -5,    0.65, -0.8 )
result  = optim(  Zstart,  Zme, control=list(maxit=1e9))
print(result)
paste('c(', paste(round(result$par, digits=3), collapse=', '), ')')
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Zvalue(result$par)), color = c(cold, hot))

# C++       Zhi0   ZhiV    ZhiT1  ZhiVk  ZmidH Zmid0 Zmidmax ZhiSH ZhiS0 ZhiSmax
# R         L      L(V)    L(T1)  L(V,k) S(H)  S(0)  S(max)  LS(H) LS(0) LS(max)
Zbest = c( 1.229, 0.07, -0.967, 1.129, 5.387, 0.772, 1.248, -4.89, 1.087, -0.209 )
Zmpe(Zbest)
Zme(Zbest)
Zbestworst = Zme(Zbest)
Zbestworst
T1[abs(Zvalue(Zbest)-Z) > Zbestworst - 0.01]
V[abs(Zvalue(Zbest)-Z) > Zbestworst - 0.01]
plot_ly(x=c(V,V), y=c(T1,T1), z=c(Z,Zvalue(Zbest)), color = c(cold, hot))

