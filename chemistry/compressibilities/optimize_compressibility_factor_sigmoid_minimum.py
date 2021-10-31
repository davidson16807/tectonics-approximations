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
    Lout = np.array([Z     for (p,T), Z in observations if p >= 5 or (p>=1.2 and T<1.05) ])
    Lin  = np.array([(p,T) for (p,T), Z in observations if p >= 5 or (p>=1.2 and T<1.05) ])
    Zout = np.array([Z     for (p,T), Z in observations])
    Zin  = np.array([(p, T) for (p,T), Z in observations])

def max_absolute_error(estimated, observed):
    return np.max(np.abs(observed-estimated))

def max_percent_absolute_error(estimated, observed):
    return np.max(np.abs(observed-estimated/observed))

def mean_percent_absolute_error(estimated, observed):
    return np.mean(np.abs(observed-estimated/observed))

def mean_absolute_error(estimated, observed):
    return np.mean(np.abs(observed-estimated))

def L(Lparams, Lin):
    p = Lin[:,0]
    T = Lin[:,1]
    V = p/T
    T1 = 1/T
    a0 = (Lparams[0])
    a1 = (Lparams[1])
    a2 = (Lparams[2])
    a3 = (Lparams[3])
    a4 = (Lparams[4])
    return a0 + a1*V**a4 + a2*T1**a3
    # return a0 + a1*p/T + a2/T + a3*(p/T)*(1/T) + a4*(p/T)**2 + a5*(1/T)**2

def Lcost1(Lparams):
    return max_absolute_error(L(Lparams, Lin), Lout)

def Lcost2(Lparams):
    return mean_absolute_error(L(Lparams, Lin), Lout)

def Lcode(Lparams):
    a0 = (Lparams[0])
    a1 = (Lparams[1])
    a2 = (Lparams[2])
    a3 = (Lparams[3])
    a4 = (Lparams[4])
    return f'{a0:.3f} {a1:+.3f}*(p/T)**{a4:+.3f} {a2:+.3f}/T**{a3:+.3f}'
    # return f'{a0:.3f} {a1:+.3f}*p/T {a2:+.3f}/T {a3:+.3f}*(p/T)*(1/T) {a4:+.3f}*(p/T)**2 {a5:+.3f}*(1/T)**2'

def Ltext(Lparams):
    arraytext = ','.join(f'{Lparams[i]:.3f}' for i in range(len(Lparams)))
    return( f'''#
# Lguess = np.array([{arraytext}])
# max error:  {Lcost1(Lparams)} 
# {Lcode(Lparams)}
# mean error: {Lcost2(Lparams)} ''')

# Lguess = np.array([1.098,0.118,-0.946,0.981,0.954])
Lguess = np.array([1.104, 0.101, -0.924, 1,1]) # best found where exponents are 1
Lsolutions = [Lguess + np.array([random.gauss(0,0.1) for j in range(len(Lguess))]) for i in range(1000000)]
Lsolutions = sorted(Lsolutions, key=Lcost1)[0:50000]
Lsolutions = genetic_algorithm([Lcost1], Ltext, Lsolutions, survival_rate=0.8, mutant_deviation=0.3)



def S(Sparams, Sin):
    p = Sin[:,0]
    T = Sin[:,1]
    V = p/T
    T1 = 1/T
    a0 = (Sparams[0])
    a1 = (Sparams[1])
    return 1/(1+np.exp(a0*(T1-a1)))

def Scode(Sparams):
    a0 = (Sparams[0])
    a1 = (Sparams[1])
    return f' 1/(1+exp({a0:.3f}*(T1-{a1:.3f})))'



def I(Iparams, Iin):
    p = Iin[:,0]
    T = Iin[:,1]
    V = p/T
    T1 = 1/T
    a0 = (Iparams[0])
    a1 = (Iparams[1])
    Lvalue = L(Iparams[2:2+5], Iin)
    Svalue = S(Iparams[2+5:2+5+2], Iin)
    return 1/(1+V*a0*np.exp((Lvalue-Svalue)*a1))

def Icode(Iparams):
    a0 = (Iparams[0])
    a1 = (Iparams[1])
    Lcodetext = Lcode(Iparams[2:2+5])
    Scodetext = Scode(Iparams[2+5:2+5+2])
    return f'1/(1+V*{a0:.3f}*np.exp(({Lcodetext}-{Scodetext})*{a1:.3f}))'



def Z(Zparams, Zin):
    Ivalue = I(Zparams, Zin)
    Lvalue = L(Zparams[2:2+5], Zin)
    return Ivalue + (1-Ivalue)*Lvalue

def Zcost1(Zparams):
    return max_absolute_error(Z(Zparams,Zin), Zout)

def Zcost2(Zparams):
    return mean_absolute_error(Z(Zparams,Zin), Zout)

def Zcode(Zparams):
    Icodetext = Icode(Zparams)
    Lcodetext = Lcode(Zparams)
    return f'({Icodetext}) + (1-{Icodetext})*({Lcodetext})'

def Ztext(Zparams):
    arraytext = ','.join(f'{Zparams[i]:.3f}' for i in range(len(Zparams)))
    return( f'''#
# Zguess = np.array([{arraytext}])
# {Zcode(Zparams)}
# max error:  {Zcost1(Zparams)} 
# mean error: {Zcost2(Zparams)} ''')

Zguess = np.array([3,3,      1.12, 0.101, -0.928, 1,1,     7.7, -0.84])
# Zguess = np.array([1.098,0.118,-0.946,0.981,0.954,      18.033,-7.974,-24.599,3.465,0.116,9.261])
# Zguess = np.array([0.103,1.245,2.083,1.030,0.994]) # best found for the other model
Zsolutions = [Zguess]+[Zguess + np.random.normal(0, 0.3, len(Zguess)) for i in range(100000)]
Zsolutions = [x for x in Zsolutions if not isnan(Zcost1(x))]
Zsolutions = sorted(Zsolutions, key=Zcost1)[0:50000]
Zsolutions = genetic_algorithm([Zcost1], Ztext, Zsolutions, survival_rate=0.8, mutant_deviation=1)

