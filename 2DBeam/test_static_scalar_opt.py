import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import minimize, minimize_scalar
from functools import partial

from solve_functions import eigenvalue_solve, static_solve
from system_functions import build_system_matrix

k_full_target, m_full_target = build_system_matrix([1.0,1.0,1.0])
eigenmodes_target, eigenfrequencies_target = eigenvalue_solve(k_full_target, m_full_target)
static_deformation_target = static_solve(k_full_target)

k_full_initial, m_full_initial = build_system_matrix([0.001,0.0012,0.0014])
#k_full_initial, m_full_initial = build_system_matrix([0.999,0.97,0.99])
eigenmodes_initial, eigenfrequencies_initial = eigenvalue_solve(k_full_initial, m_full_initial)
static_deformation_initial = static_solve(k_full_initial)

x = [0.0, 1.0, 2.0, 3.0]


def evaluate_residual(a_cur, a_tar):
    return np.linalg.norm(np.subtract(a_cur, a_tar))/np.amax(np.absolute(a_tar))

def objective_function(static_deformation_target, sys_param):

    k_full_cur, m_full_cur = build_system_matrix([sys_param, 1.0, 1.0])
    static_deformation_cur = static_solve(k_full_cur)

    f1 = evaluate_residual(static_deformation_cur['y'],static_deformation_target['y'])
    f2 = evaluate_residual(static_deformation_cur['g'],static_deformation_target['g'])
    f = 0.67*f1**2+0.33*f2**2

    print('F: ', str(f))

    return f

optimizable_function = partial(objective_function, static_deformation_target)

res_scalar = minimize_scalar(optimizable_function, options={'gtol': 1e-6, 'disp': True})
print(res_scalar.x)