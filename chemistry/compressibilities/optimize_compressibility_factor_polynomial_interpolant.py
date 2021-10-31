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
    Fout = np.array([Z     for (p,T), Z in observations if p >= 5 or (p>=1.2 and T<1.05) ])
    Fin  = np.array([(p,T) for (p,T), Z in observations if p >= 5 or (p>=1.2 and T<1.05) ])
    Zout = np.array([Z     for (p,T), Z in observations])
    Zin  = np.array([(p, T, 1, (p/T), (1/T), (p/T)*(1/T), (p/T)**2, (1/T)**2) for (p,T), Z in observations])

def max_absolute_error(estimated, observed):
    return np.max(np.abs(observed-estimated))

def max_percent_absolute_error(estimated, observed):
    return np.max(np.abs(observed-estimated/observed))

def mean_percent_absolute_error(estimated, observed):
    return np.mean(np.abs(observed-estimated/observed))

def mean_absolute_error(estimated, observed):
    return np.mean(np.abs(observed-estimated))

def F(Fparams, Fin):
    p = Fin[:,0]
    T = Fin[:,1]
    a = (Fparams[0])
    b = (Fparams[1])
    g = (Fparams[2])
    h = (Fparams[3])
    k = (Fparams[4])
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

def Fcost1(Fparams):
    return max_absolute_error(F(Fparams, Fin), Fout)

def Fcost2(Fparams):
    return mean_absolute_error(F(Fparams, Fin), Fout)

def Fcode(Fparams):
    a = (Fparams[0])
    b = (Fparams[1])
    g = (Fparams[2])
    h = (Fparams[3])
    k = (Fparams[4])
    return f'{a:.3f} {b:+.3f}*p/T + {g:+.3f}/T'
# exp(-{b:.3f}**2/T**{g:.3f}) + {a:.3f} * p**{k:.3f}/T**{h:.3f}
# exp(-{b:.3f}**2/T**2) + {a:.3f}/T**2 * p**1

def Ftext(Fparams):
    arraytext = ','.join(f'{Fparams[i]:.3f}' for i in range(len(Fparams)))
    return( f'''#
# Fguess = np.array([{arraytext}])
# max error:  {Fcost1(Fparams)} 
# {Fcode(Fparams)}
# mean error: {Fcost2(Fparams)} ''')

Fguess = np.array([1.110,0.105,-0.942,-0.000,-4.432])
Fsolutions = [Fguess + np.array([random.gauss(0,0.1) for j in range(len(Fguess))]) for i in range(50000)]
Fsolutions = genetic_algorithm([Fcost1], Ftext, Fsolutions, survival_rate=0.8, mutant_deviation=0.1)



def G(Gparams, Gin):
    return np.exp(np.matmul(Gin,  Gparams))

def Gcode(Gparams):
    a0 = (Gparams[0])
    a1 = (Gparams[1])
    a2 = (Gparams[2])
    a3 = (Gparams[3])
    a4 = (Gparams[4])
    a5 = (Gparams[5])
    return f'exp({a0:.3f} {a1:+.3f}*(p/T) {a2:+.3f}*(1/T) {a3:+.3f}*(p/T)*(1/T) {a4:+.3f}*(p/T)**2  {a5:+.3f}*(1/T)**2 )'

def I(Iparams, Iin):
    return np.matmul(Iin,  Iparams[:6]) /  np.matmul(Iin,  Iparams[6:])

def Z(Zparams, Zin):
    # F2 = F(Zparams[:5], Zin[:,:2])
    F2 = F(Fguess, Zin[:,:2]) # use the best parameterization, parameterize separately
    G2 = G(Zparams[5:5+6], Zin[:,2:])
    return G2 + (1-G2)*F2

def Zcost1(Zparams):
    return max_absolute_error(Z(Zparams,Zin), Zout)

def Zcost2(Zparams):
    return mean_absolute_error(Z(Zparams,Zin), Zout)

def Zcode(Zparams):
    # Fcodetext = Fcode(Zparams[:5])  
    Fcodetext = Fcode(Fguess)   # use the best parameterization, parameterize separately
    Gcodetext = Gcode(Zparams[5:5+6])
    return f'({Gcodetext}) + (1-{Gcodetext})*({Fcodetext})'

def Ztext(Zparams):
    arraytext = ','.join(f'{Zparams[i]:.3f}' for i in range(len(Zparams)))
    return( f'''#
# Zguess = np.array([{arraytext}])
# {Zcode(Zparams)}
# max error:  {Zcost1(Zparams)} 
# mean error: {Zcost2(Zparams)} ''')

Zguess = np.array([1.104, 0.101, -0.924, 13.249, -4.144,      1.151, -0.477, -0.309, 0.194, -0.219, 0.123])
# Zguess = np.array([0.103,1.245,2.083,1.030,0.994]) # best found for the other model
Zsolutions = [Zguess + np.random.normal(0, 0.3, len(Zguess)) for i in range(1000000)]
Zsolutions = [x for x in Zsolutions if not isnan(Zcost1(x))]
Zsolutions = sorted(Zsolutions, key=Zcost1)[0:30000]
Zsolutions = genetic_algorithm([Zcost1], Ztext, Zsolutions, survival_rate=0.8, mutant_deviation=1)

