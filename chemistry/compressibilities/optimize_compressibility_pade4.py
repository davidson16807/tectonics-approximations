from math import *
import csv
import random
import numpy as np

from optimize import OptimizationBasis, genetic_algorithm

numerator_basis = [
    '1', 
    'p', 
    'T', 
    'p**2', 
    'p*T', 
    'T**2', 
    'p**3', 
    'p**2*T',
    'p*T**2',
    'T**3',
    'p**4',
    'p**3*T',
    'p**2*T**2',
    'p*T**3',
    'T**4',
    # 'sqrt(p*p+T*T)',
    #'sqrt(p)',
    #'sqrt(T)',
]
denominator_basis = [
    '1', 
    'p', 
    'T', 
    'p**2', 
    'p*T', 
    'T**2', 
    'p**3', 
    'p**2*T',
    'p*T**2',
    'T**3',
    'p**4',
    'p**3*T',
    'p**2*T**2',
    'p*T**3',
    'T**4',
    #'sqrt(p*p+T*T)',
    #'sqrt(p)',
    #'sqrt(T)',
]

numerator_mask   = np.array([i <  len(numerator_basis) for i in range(parameter_count)])
denominator_mask = np.array([i >= len(numerator_basis) for i in range(parameter_count)])

with open('pTZ.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(csvreader, None) # skip header
    observations = [( [eval(x, {}, {'p':float(p),'T':float(T)}) for x in (numerator_basis + denominator_basis)], float(Z))
                      for (i,p,Z,T) in csvreader ]
    output = np.array([Z     for basis, Z in observations])
    input_ = np.array([basis for basis, Z in observations])

def pade_cost(coefficients):
    pade = np.dot(input_*coefficients, numerator_mask) / np.dot(input_*(coefficients), denominator_mask)
    error = output - pade
    # return np.sum(np.abs(error) / output) / len(output)
    return np.max(np.abs(error))

def pade_mean(coefficients):
    pade = np.dot(input_*coefficients, numerator_mask) / np.dot(input_*(coefficients), denominator_mask)
    error = output - pade
    return np.sum(np.abs(error) / output) / len(output)
    # return np.max(np.abs(error) / output)

def code(coefficients):
    a = coefficients[0:len(numerator_basis)]
    b = coefficients[len(numerator_basis):(len(numerator_basis)+len(denominator_basis))]
    numerator = '+'.join(f'{a[i]:.3f}*{numerator_basis[i]}' for i in range(len(a)))
    denominator = '+'.join(f'{(b[i]):.3f}*{denominator_basis[i]}' for i in range(len(b)))
    return(f'({numerator}) / ({denominator})')

def pade_basis_text(basis):
    return basis.replace('**3','³').replace('**2','²').replace('*','').replace('4','⁴').replace('sqrt(','√').replace(')','')

def pade_text(coefficients):
    a = coefficients[0:len(numerator_basis)]
    b = coefficients[len(numerator_basis):(len(numerator_basis)+len(denominator_basis))]
    numerator = ' '.join(f'{a[i]:+.2f}{pade_basis_text(numerator_basis[i])}' for i in range(len(a)))
    denominator = ' '.join(f'{(b[i]):+.2f}{pade_basis_text(denominator_basis[i])}' for i in range(len(b)))
    paramtext = ','.join(f'{coefficients[i]:.3f}' for i in range(len(coefficients)))
    return(
f"""# {numerator} 
# {'-'*max(len(numerator),len(denominator))}
# {denominator}
# COST:   {pade_cost(coefficients)}
# MEAN E: {pade_mean(coefficients)}
# PARAMS: np.array([{paramtext}])
# CODE:   {code(coefficients)}
#""")

basis = OptimizationBasis(pade_cost, pade_text)

guess = np.array([2.543,-0.868,1.356,3.688,-3.646,-3.515,-2.956,-2.788,-1.917,-0.857,0.391,2.000,2.869,-0.857,4.861,3.178,-0.144,0.083,1.266,-2.868,-2.666,0.020,-3.560,-1.470,-1.129,0.024,2.721,0.316,0.452,4.699])
solutions = [guess]+[guess+np.array([random.gauss(0,1) for j in range(len(numerator_basis)+len(denominator_basis))]) for i in range(30000)]
solutions = genetic_algorithm(basis, solutions, survival_rate=0.8, mutation_offset=1)



# BEST:
# +2.541 -0.87p +1.36T +3.69p² -3.65pT -3.52T² -2.96p³ -2.79p²T -1.92pT² -0.86T³ +0.39p⁴ +2.00p³T +2.87p²T² -0.86pT³ +4.86T⁴ 
# --------------------------------------------------------------------------------------------------------------------------
# +3.181 -0.14p +0.08T +1.27p² -2.87pT -2.67T² +0.02p³ -3.56p²T -1.47pT² -1.13T³ +0.02p⁴ +2.72p³T +0.32p²T² +0.45pT³ +4.70T⁴
# COST:   0.08238938225311054
# MEAN E: 0.03814189990107604
# PARAMS: np.array([2.543,-0.868,1.356,3.688,-3.646,-3.515,-2.956,-2.788,-1.917,-0.857,0.391,2.000,2.869,-0.857,4.861,3.178,-0.144,0.083,1.266,-2.868,-2.666,0.020,-3.560,-1.470,-1.129,0.024,2.721,0.316,0.452,4.699])
# CODE:   (2.543*1+-0.868*p+1.356*T+3.688*p**2+-3.646*p*T+-3.515*T**2+-2.956*p**3+-2.788*p**2*T+-1.917*p*T**2+-0.857*T**3+0.391*p**4+2.000*p**3*T+2.869*p**2*T**2+-0.857*p*T**3+4.861*T**4) / (3.178*1+-0.144*p+0.083*T+1.266*p**2+-2.868*p*T+-2.666*T**2+0.020*p**3+-3.560*p**2*T+-1.470*p*T**2+-1.129*T**3+0.024*p**4+2.721*p**3*T+0.316*p**2*T**2+0.452*p*T**3+4.699*T**4)
#
