import random
import numpy as np


def gradient_descent_algorithm(cost_getter, text_getter, x, mobility=3e-5, delta=1e-5):
    def gradient(a, cost, delta=1e-5):
        return np.array([ (cost(x+dx_i)-cost(x-dx_i)) / (2.*np.linalg.norm(dx_i)) for dx_i in delta*np.identity(len(x)) ])
    def gradient_descent(x, cost, mobility=3e-6, delta=1e-8):
        return x + mobility * gradient(x, cost, delta)
    best, best_cost = x, float('inf')
    grad_cost = 0.*x
    try:
        while True:
            # grad_cost = gradient(x, cost_getter, delta)
            cost_x = cost_getter(x)
            # if cost_x < best_cost: print('NEW PERSONAL BEST!\a')
            best, best_cost = (x, cost_x) if cost_x < best_cost else (best, best_cost)
            x -= mobility * gradient(x, cost_getter, delta)
            print(text_getter(x))
            print('cost: ', cost_x)
            # print('∇cost:', grad_cost)
    except KeyboardInterrupt as e:
        return best

def genetic_algorithm(
        cost_getters, 
        text_getter, 
        solutions, 
        iterations=float('inf'), 
        survival_rate=0.5, 
        mutant_deviation=0.3, 
        mutant_deviations=None, 
        scale_mutation_if_slow = 0.99):
    mutant_deviations = (0*solutions[0]) + mutant_deviation if mutant_deviations is None else mutant_deviations
    def mate(a,b):
        return np.array([random.choice([ai,bi]) for ai,bi in zip(a,b)])
    def mutate(a):
        return np.array([ai + random.gauss(0,mutant_deviations[i]) for i,ai in enumerate(a)])
    def select(solutions):
        return solutions[int(random.paretovariate(1)) % len(solutions)]
    best_costs = [float('inf') for j in cost_getters]
    i = 0
    try:
        while i < iterations:
            i+=1
            for j, cost_getter in enumerate(cost_getters):
                solutions = sorted(solutions, key=cost_getter)
                cutoff = int(survival_rate*len(solutions))
                solutions[cutoff:len(solutions)] = [
                    mutate(mate(select(solutions), select(solutions)))
                    for k in solutions[cutoff:len(solutions)]
                ]
                new_best_cost = cost_getter(solutions[0])
                if new_best_cost < best_costs[j] and iterations>1:
                    best_costs[j] = new_best_cost
                    print('NEW PERSONAL BEST!\a')
                    mutant_deviations /= scale_mutation_if_slow**2
                    print(f'mutant deviations: {mutant_deviations}')
                else:
                    mutant_deviations *= scale_mutation_if_slow
                    print(f'mutant deviations: {mutant_deviations}')
                print(text_getter(solutions[0]))
    except KeyboardInterrupt as e:
        return solutions
    return solutions

