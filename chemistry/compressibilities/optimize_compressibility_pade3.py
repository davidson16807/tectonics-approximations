from math import *
import csv
import random
import numpy as np

from optimize import OptimizationBasis, genetic_algorithm

def get_basis(p,T):
    return np.array([
            #NUMERATOR:
            1, 
            p, 
            T, 
            p**2, 
            p*T, 
            T**2, 
            p**3, 
            p**2*T,
            p*T**2,
            T**3,
            # p**4,
            # p**3*T,
            # p**2*T**2,
            # p*T**3,
            # T**4,
            #sqrt(p*p+T*T),
            #sqrt(p),
            #sqrt(T),

            #DENOMINATOR:
            1, 
            p, 
            T, 
            p**2, 
            p*T, 
            T**2, 
            p**3, 
            p**2*T,
            p*T**2,
            T**3,
            # p**4,
            # p**3*T,
            # p**2*T**2,
            # p*T**3,
            # T**4,
            #sqrt(p*p+T*T),
            #sqrt(p),
            #sqrt(T),
        ])

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
    # 'p**4',
    # 'p**3*T',
    # 'p**2*T**2',
    # 'p*T**3',
    # 'T**4',
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
    # 'p**4',
    # 'p**3*T',
    # 'p**2*T**2',
    # 'p*T**3',
    # 'T**4',
    #'sqrt(p*p+T*T)',
    #'sqrt(p)',
    #'sqrt(T)',
]

parameter_count = len(numerator_basis)+len(denominator_basis)

numerator_mask   = np.array([i <  len(numerator_basis) for i in range(parameter_count)])
denominator_mask = np.array([i >= len(numerator_basis) for i in range(parameter_count)])

# np.row_stack((a,b))

with open('pTZ.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(csvreader, None) # skip header
    observations = [( get_basis(float(p),float(T)), float(Z))
                      for (i,p,Z,T) in csvreader ]
    output = np.array([Z     for basis, Z in observations])
    input_ = np.array([basis for basis, Z in observations])

def pade_cost(coefficients):
    terms = coefficients*input_
    pade = np.dot(terms,numerator_mask) / np.dot(terms,denominator_mask)
    error = output - pade
    return np.sum(np.abs(error) / output) / len(output)

def code(params):
    a = params[0:len(numerator_basis)]
    b = params[len(numerator_basis):(len(numerator_basis)+len(denominator_basis))]
    numerator = '+'.join(f'{a[i]:.3f}*{numerator_basis[i]}' for i in range(len(a)))
    denominator = '+'.join(f'{exp(b[i]):.3f}*{denominator_basis[i]}' for i in range(len(b)))
    return(f'({numerator}) / ({denominator})')

def pade_basis_text(basis):
    return basis.replace('**3','³').replace('**2','²').replace('*','').replace('4','⁴').replace('sqrt(','√').replace(')','')

def pade_text(params):
    a = params[0:len(numerator_basis)]
    b = params[len(numerator_basis):(len(numerator_basis)+len(denominator_basis))]
    numerator = ' '.join(f'{a[i]:+.2f}{pade_basis_text(numerator_basis[i])}' for i in range(len(a)))
    denominator = ' '.join(f'{exp(b[i]):+.2f}{pade_basis_text(denominator_basis[i])}' for i in range(len(b)))
    return(
f"""# {numerator} 
# {'-'*max(len(numerator),len(denominator))}
# {denominator}
# cost: {pade_cost(params)}
#""")


basis = OptimizationBasis(pade_cost, pade_text)


solutions = [np.array([random.gauss(0,1) for j in range(len(numerator_basis)+len(denominator_basis))]) for i in range(20000)]
solutions = genetic_algorithm(basis, solutions, survival_rate=0.8, mutation_offset=0.3)










# BEST:
# -5.591 +0.96p +0.59T +5.20p² +4.47pT +0.70T² -0.53p³ -5.96p²T -0.08pT² -4.60T³ 
# ----------------------------------------------------------------------------------
# +0.011 +0.12p +0.02T +1.65p² +2004.64pT +49.31T² +0.99p³ +0.00p²T +1.54pT² +0.00T³
# cost: 0.05457421662157773
# params: np.array([-5.589,0.963,0.588,5.197,4.468,0.702,-0.534,-5.960,-0.084,-4.602,-4.377,-2.086,-3.846,0.501,7.603,3.898,-0.007,-5.621,0.433,-5.356])
# code: (-5.589*1+0.963*p+0.588*T+5.197*p**2+4.468*p*T+0.702*T**2+-0.534*p**3+-5.960*p**2*T+-0.084*p*T**2+-4.602*T**3) / (0.013*1+0.124*p+0.021*T+1.651*p**2+2004.641*p*T+49.314*T**2+0.993*p**3+0.004*p**2*T+1.541*p*T**2+0.005*T**3)
