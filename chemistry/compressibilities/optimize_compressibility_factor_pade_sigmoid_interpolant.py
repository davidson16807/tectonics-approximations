from math import *
import csv
import random
import numpy as np
import plotly.graph_objects as go

from optimize import genetic_algorithm, gradient_descent_algorithm

def Gobserved(Z, p, T):
    try:
        F = 1.104 + 0.101*p/T - 0.924/T
        return (Z-F)/(-F+1)
        # return -log(1/I+1)
    except ValueError:
        return None

with open('pTZ.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(csvreader, None) # skip header
    observations = [( np.array([float(p),float(T)]), float(Z))
                      for (i,p,Z,T) in csvreader ]

Fout = np.array([Z     for (p,T), Z in observations if p >= 5 or (p>=1.2 and T<1.05) ])
Fin  = np.array([(p,T) for (p,T), Z in observations if p >= 5 or (p>=1.2 and T<1.05) ])
Zout = np.array([Z     for (p,T), Z in observations])
Zin  = np.array([(p, T, 1, (p/T), (1/T), (p/T)*(1/T), (p/T)**2, (1/T)**2) for (p,T), Z in observations])

Gout = np.array([ Gobserved(Z, p, T)
    for (p,T), Z in observations if p<5 and T<2 and Gobserved(Z, p, T) is not None])
Gin  = np.array([(1, (p/T), (1/T), (p/T)*(1/T), (p/T)**2, (1/T)**2)
    for (p,T), Z in observations if p<5 and T<2 and Gobserved(Z, p, T) is not None])

Gp = Gin[:,1] / Gin[:,2]
GT = 1 / Gin[:,2]

Gref = [go.Scatter3d(x=Gp**2/(GT**6 + 1/GT**2), y=GT, z=-np.log(Gout), mode='markers')]
f = go.Figure(data= Gref)
f.show()

Gref = [go.Scatter3d(x=Gp/GT, y=GT, z=(-np.log(Gout)), mode='markers')]
f = go.Figure(data= Gref)
f.show()


Gref = [go.Scatter3d(x=Gp/GT, y=GT, z=-np.log(Gout)-0.3*Gp/GT, mode='markers')]
Gfit1 = [go.Scatter3d(x=Gp/GT, y=GT, z=0.3*Gp/GT, mode='markers')]
Gfit2 = [go.Scatter3d(x=Gp/GT, y=GT, z=GT/Gp, mode='markers')]
f = go.Figure(data= Gref)
f.show()


Gref = [go.Scatter3d(x=Gp, y=GT, z=(Gout), mode='markers')]
Gfit = [go.Scatter3d(x=Gp, y=GT, z=np.exp(-(1/(1.1*GT**6/Gp**2 + 0.9*Gp**1/GT**1))), mode='markers')]
f = go.Figure(data= Gref+Gfit)
f.show()



# Gfit =[go.Scatter3d(visible=False, x=((Gp)/GT**3)**i, y=GT, z=-np.log(Gout), mode='markers') for i in np.arange(0, 2, 0.1)]
# steps = [dict(method="update", args=[{"visible": [j==i for j in range(len(f.data))]}],) for i in range(len(f.data))]
# f.update_layout(
#     sliders=[dict(
#         active=10,
#         currentvalue={"prefix": "Frequency: "},
#         pad={"t": 50},
#         steps=steps
#     )]
# )
# f.show()

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
    a0 = (Fparams[0])
    a1 = (Fparams[1])
    a2 = (Fparams[2])
    a3 = (Fparams[3])
    a4 = (Fparams[4])
    return a0 + a1*(p/T)**a4 + a2/T**a3
    # return a0 + a1*p/T + a2/T + a3*(p/T)*(1/T) + a4*(p/T)**2 + a5*(1/T)**2

def Fcost1(Fparams):
    return max_absolute_error(F(Fparams, Fin), Fout)

def Fcost2(Fparams):
    return mean_absolute_error(F(Fparams, Fin), Fout)

def Fcode(Fparams):
    a0 = (Fparams[0])
    a1 = (Fparams[1])
    a2 = (Fparams[2])
    a3 = (Fparams[3])
    a4 = (Fparams[4])
    return f'{a0:.3f} {a1:+.3f}*(p/T)**{a4:+.3f} {a2:+.3f}/T**{a3:+.3f}'
    # return f'{a0:.3f} {a1:+.3f}*p/T {a2:+.3f}/T {a3:+.3f}*(p/T)*(1/T) {a4:+.3f}*(p/T)**2 {a5:+.3f}*(1/T)**2'

def Ftext(Fparams):
    arraytext = ','.join(f'{Fparams[i]:.3f}' for i in range(len(Fparams)))
    return( f'''#
# Fguess = np.array([{arraytext}])
# max error:  {Fcost1(Fparams)} 
# {Fcode(Fparams)}
# mean error: {Fcost2(Fparams)} ''')

Fguess = np.array([1.098,0.118,-0.946,0.981,0.954])
# Fguess = np.array([1.104, 0.101, -0.924, 1,1]) # best found where exponents are 1
Fsolutions = [Fguess + np.array([random.gauss(0,0.1) for j in range(len(Fguess))]) for i in range(1000000)]
Fsolutions = sorted(Fsolutions, key=Fcost1)[0:50000]
Fsolutions = genetic_algorithm([Fcost1], Ftext, Fsolutions, survival_rate=0.8, mutant_deviation=0.3)



def G(Gparams, Gin):
    return 1/(1+np.exp(-np.matmul(Gin,  Gparams[:6]) / np.matmul(Gin,  Gparams[6:])))

def Gcode(Gparams):
    a0 = (Gparams[0])
    a1 = (Gparams[1])
    a2 = (Gparams[2])
    a3 = (Gparams[3])
    a4 = (Gparams[4])
    a5 = (Gparams[5])
    b0 = (Gparams[6])
    b1 = (Gparams[7])
    b2 = (Gparams[8])
    b3 = (Gparams[9])
    b4 = (Gparams[10])
    b5 = (Gparams[11])
    return f'''1/
(1+exp(-({a0:.3f} {a1:+.3f}*(p/T) {a2:+.3f}*(1/T) {a3:+.3f}*(p/T)*(1/T) {a4:+.3f}*(p/T)**2  {a5:+.3f}*(1/T)**2)  /  
        ({b0:.3f} {b1:+.3f}*(p/T) {b2:+.3f}*(1/T) {b3:+.3f}*(p/T)*(1/T) {b4:+.3f}*(p/T)**2  {b5:+.3f}*(1/T)**2) ))
'''

def Gcost1(Gparams):
    return max_absolute_error(G(Gparams,Gin), Gout)

def Gcost2(Gparams):
    return mean_absolute_error(G(Gparams,Gin), Gout)

def Gtext(Gparams):
    arraytext = ','.join(f'{Gparams[i]:.3f}' for i in range(len(Gparams)))
    return( f'''#
# Gguess = np.array([{arraytext}])
# I2={Gcode(Gparams)}
# max error:  {Gcost1(Gparams)} 
# mean error: {Gcost2(Gparams)} ''')

Gguess = np.array([
    17.341,-7.712,-23.699,3.872,-0.404,9.205, 
    1, 0,0,0,0,0
])

Ginitial_count = 10**8
Gsample_count = 10**5
Gsample = [np.random.normal(0, 10, len(Gguess)) for i in range(Gsample_count)]
representative_count = int(len(Gsample) * len(Gsample)/Ginitial_count)
Gcost_threshold = Gcost1(sorted(Gsample, key=Gcost1)[representative_count])
Gsolutions = []
for i in range(Ginitial_count):
    if i % 10**6 == 0:
        print(i)
    Gsolution = np.random.normal(0, 1, len(Gguess))
    cost = Gcost1(Gsolution)
    if cost is not None and cost < Gcost_threshold:
        Gsolutions.append(Gsolution)


Gsolutions_cost = [Gcost1(x) for x in Gsolutions]
    Gsolutions, Gsolutions_cost = best_half(Gsolutions, Gsolutions_cost)



Gsolutions = genetic_algorithm([Gcost1], Gtext, Gsolutions, survival_rate=0.8, mutant_deviation=1, scale_mutation_if_slow = 0.999)


def Z(Zparams, Zin):
    # F2 = F(Zparams[:5], Zin[:,:2])
    F2 = F(Fguess, Zin[:,:2]) # use the best parameterization, parameterize separately
    G2 = G(Zparams[5:], Zin[:,2:])
    return G2 + (1-G2)*F2

def Zcost1(Zparams):
    return max_absolute_error(Z(Zparams,Zin), Zout)

def Zcost2(Zparams):
    return mean_absolute_error(Z(Zparams,Zin), Zout)

def Zcode(Zparams):
    # Fcodetext = Fcode(Zparams[:5])  
    Fcodetext = Fcode(Fguess)   # use the best parameterization, parameterize separately
    Gcodetext = Gcode(Zparams[5:])
    return f'({Gcodetext}) + (1-{Gcodetext})*({Fcodetext})'

def Ztext(Zparams):
    arraytext = ','.join(f'{Zparams[i]:.3f}' for i in range(len(Zparams)))
    return( f'''#
# Zguess = np.array([{arraytext}])
# {Zcode(Zparams)}
# max error:  {Zcost1(Zparams)} 
# mean error: {Zcost2(Zparams)} ''')

Zguess = np.array([
    1.037,0.106,-0.919,1.542,0.985,
    17.341,-7.712,-23.699,3.872,-0.404,9.205, 
    1, 0,0,0,0,0])
# Zguess = np.array([1.098,0.118,-0.946,0.981,0.954,      18.033,-7.974,-24.599,3.465,0.116,9.261])
# Zguess = np.array([0.103,1.245,2.083,1.030,0.994]) # best found for the other model
Zsolutions = [Zguess]+[Zguess + np.random.normal(0, 0.3, len(Zguess)) for i in range(1000000)]
Zsolutions = [x for x in Zsolutions if not isnan(Zcost1(x))]
Zsolutions = sorted(Zsolutions, key=Zcost1)[0:50000]
Zsolutions = genetic_algorithm([Zcost1], Ztext, Zsolutions, survival_rate=0.8, mutant_deviation=1, scale_mutation_if_slow = 0.999)

