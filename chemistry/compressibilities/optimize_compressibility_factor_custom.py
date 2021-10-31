from math import *
import csv
import random
import numpy as np

from optimize import genetic_algorithm

with open('pTZ.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(csvreader, None) # skip header
    observations = [( np.array([float(p),float(T)]), float(Z))
                      for (i,p,Z,T) in csvreader ]
    Aout = np.array([Z     for (p,T), Z in observations if p >= 5 or (p>=1.2 and T<1.05) ])
    Ain  = np.array([(p,T) for (p,T), Z in observations if p >= 5 or (p>=1.2 and T<1.05) ])
    Iout = np.array([Z     for (p,T), Z in observations])
    Iin  = np.array([(p,T) for (p,T), Z in observations])

def max_absolute_error(estimated, observed):
    return np.max(np.abs(observed-estimated))

def max_percent_absolute_error(estimated, observed):
    return np.max(np.abs(observed-estimated/observed))

def mean_percent_absolute_error(estimated, observed):
    return np.mean(np.abs(observed-estimated/observed))

def mean_absolute_error(estimated, observed):
    return np.mean(np.abs(observed-estimated))

def A(Aparams, Ain):
    p = Ain[:,0]
    T = Ain[:,1]
    a = (Aparams[0])
    b = (Aparams[1])
    g = (Aparams[2])
    h = (Aparams[3])
    k = (Aparams[4])
    # NOTE: m is the slope of a crude linear relation between p and Z 
    # between the min and max observations for Z at a given T.
    # In other words, `m=40-2.6/(Z1-Z0); lm(m~T)` 
    # where Z0 is the minimum Z observed for a p and Z1 is the Z observed at p=40.
    # we approximate the minimum as occuring at p=2.65 for the samples we consider
    m = h*T+k
    t = a*T+b
    return a + b*p/T + g/T
    # return np.exp(-b**2 / T**g) + a * p**k/T**h
    # return np.exp(-b*b*np.power(T, -2)) + a*np.power(T, -1) * np.power(p,1)

def Acost1(Aparams):
    return max_absolute_error(A(Aparams, Ain), Aout)

def Acost2(Aparams):
    return mean_absolute_error(A(Aparams, Ain), Aout)

def Acode(Aparams):
    a = (Aparams[0])
    b = (Aparams[1])
    g = (Aparams[2])
    h = (Aparams[3])
    k = (Aparams[4])
    return f'{a:.3f} {b:+.3f}*p/T + {g:+.3f}/T'
# exp(-{b:.3f}**2/T**{g:.3f}) + {a:.3f} * p**{k:.3f}/T**{h:.3f}
# exp(-{b:.3f}**2/T**2) + {a:.3f}/T**2 * p**1

def Atext(Aparams):
    arraytext = ','.join(f'{Aparams[i]:.3f}' for i in range(len(Aparams)))
    return( f'''#
# Aguess = np.array([{arraytext}])
# max error:  {Acost1(Aparams)} 
# {Acode(Aparams)}
# mean error: {Acost2(Aparams)} ''')

Aguess = np.array([1.104, 0.101, -0.924, 13.249, -4.144])
# Aguess = np.array([0.103,1.245,2.083,1.030,0.994]) # best found for the other model
Asolutions = [Aguess + np.array([random.gauss(0,0.1) for j in range(len(Aguess))]) for i in range(50000)]
Asolutions = genetic_algorithm([Acost1], Atext, Asolutions, survival_rate=0.8, mutant_deviation=0.1)

# Aguess = np.array([1.579,-0.588,2.619,13.293,-3.782])
# max error:  0.13133828264154812 
# exp(-(1.579*T-0.588)**-2.619) + (p-2.65)/(13.293*T-3.782)
# mean error: 0.060288342807478766 


def I(Iparams):
    p = Iin[:,0]
    T = Iin[:,1]
    c = (Iparams[0])
    d = (Iparams[1])
    f = (Iparams[2])
    m = (Iparams[3])
    n = (Iparams[4])
    o = (Iparams[5])
    # return np.exp(-c*p/T -d*T -f*T )
    return np.exp(- p**m/(T**n + c) - d*p/T )
    # return np.exp(-c*c*T**-f * p**d - m*m*T**-n * p**o)
    # return np.exp(-c*T**-2 * p**2 - m*T**-1 * p**1)

def Icost1(Iparams):
    Iestimate = I(Iparams)
    return max_absolute_error(Iestimate + (1-Iestimate)*A(Abest,Iin), Iout)

def Icost2(Iparams):
    Iestimate = I(Iparams)
    return mean_absolute_error(Iestimate + (1-Iestimate)*A(Abest,Iin), Iout)

def Icode(Iparams):
    c = (Iparams[0])
    d = (Iparams[1])
    f = (Iparams[2])
    m = (Iparams[3])
    n = (Iparams[4])
    o = (Iparams[5])
    # Icodetext = f'exp(-{c:.3f}*p/T -{d:.3f}*T -{f:.3f}/T )'
    Icodetext = f'exp(-p**{m:.3f}/(T**{n:.3f} + {c:.3f}) - {d:.3f}*p/T )'
    # Icodetext = f'exp(-{c**2:.3f}/T**{f:.3f}*p**{d:.3f}-{m**2:.3f}/T**{n:.3f}*p**{o:.3f})'
    # Icodetext = f'exp(-{c:.3f}/T**2*p**2-{m:.3f}/T*p)'
    Acodetext = Acode(Abest)
    return f'{Icodetext} + (1-{Icodetext}) * ({Acodetext})'

def Itext(Iparams):
    arraytext = ','.join(f'{Iparams[i]:.3f}' for i in range(len(Iparams)))
    return( f'''#
# Iguess = np.array([{arraytext}])
# {Icode(Iparams)}
# max error:  {Icost1(Iparams)} 
# mean error: {Icost2(Iparams)} ''')

Abest = np.array([1.104, 0.101, -0.924, 13.249, -4.144])
Iguess = np.array([0.5, 0.45, 0,  4, 16, 0])
Isolutions = [Iguess] + [Iguess + np.array([random.gauss(0,1) for j in range(len(Iguess))]) for i in range(50000)]
Isolutions = genetic_algorithm([Icost1], Itext, Isolutions, survival_rate=0.8, mutant_deviation=1)



# Aguess = np.array([0.103,1.245,2.083,1.030,0.994])
# exp(-1.245**2/T**2.083) + 0.103/T**1.030 * p**0.994
# cost: 0.05995847040165758 
#
# Iguess = np.array([-0.007,0.751,0.357,-1.277,0.964,1.776])
# [0.103,1.245,-0.007,0.751,0.357,2.083]
# exp(-0.000/T**0.357*p**0.751-1.630/T**0.964*p**1.776) + (1-exp(-0.000/T**0.357*p**0.751-1.630/T**0.964*p**1.776)) * ( exp(-1.550/T**2.083) + 0.103/T**1.030 * p**0.994 )
# cost: 0.22018513974648002 
#
# We initially added the exponent parameters to avoid imposing any unnecessary restrictions on optimization.
# The fact that the exponents, "2.083", "1.030", and "0.994" are so close to whole numbers however suggests to us 
# some deeper dimensional significance which we have yet to determine. Coupled with our desire to reduce 
# the dimensionality of our problem space, we set these values to the nearest whole numbers and 
# ran the optimization again:



# We also noted the estimate overpredicted for the P<5 range, 
# so we tried adding data below P<5 for T=1 where the relationship was most clearly linear
# Here's the result with the exponents:
#
# Aguess = np.array([0.111,1.379,2.359,1.051,0.977])
# exp(-1.379**2/T**2.359) + 0.111/T**1.051 * p**0.977
# max error:  0.06661600088642097 
# mean error: 0.030296935164155915 
#
# and without:
#
#

# Iguess = np.array([0.289,0.855,1.212,19.057,58.592,3.030])
# [0.111,1.380,0.289,0.855,1.212,2.375]
# exp(-p**19.057/(T**58.592 + 0.289) - 0.855*p/T ) + (1-exp(-p**19.057/(T**58.592 + 0.289) - 0.855*p/T )) * ( exp(-1.380**2/T**2.375) + 0.111 * p**0.977/T**1.053 )
# max error:  0.1390981406346533 
# mean error: 0.05022545857049926 
# cost:  0.13917318918068738
