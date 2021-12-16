import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import minimize, minimize_scalar
from functools import partial

from solve_functions import eigenvalue_solve, static_solve
from system_functions import build_system_matrix


def evaluate_residual(a_cur, a_tar):
    return np.linalg.norm(np.subtract(a_cur, a_tar))/np.amax(np.absolute(a_tar))


def check_and_flip_sign(eigenmodes):

    for idx in range(len(eigenmodes['y'])):
        if eigenmodes['y'][idx][1] < 0:
            eigenmodes['y'][idx] = np.negative(eigenmodes['y'][idx])
            eigenmodes['g'][idx] = np.negative(eigenmodes['g'][idx])

    return eigenmodes


k_full_target, m_full_target = build_system_matrix([1.0,1.0,1.0])
eigenmodes_target, eigenfrequencies_target = eigenvalue_solve(k_full_target, m_full_target)
eigenmodes_target = check_and_flip_sign(eigenmodes_target)
static_deformation_target = static_solve(k_full_target)

x = [0.0, 1.0, 2.0, 3.0]

# seems to work robustly for modes 0 and 1
consider_mode = 0

def objective_function(mode, eigenmodes_target, eigenfrequencies_target, sys_params):

    k_full_cur, m_full_cur = build_system_matrix(sys_params)
    eigenmodes_cur, eigenfrequencies_cur = eigenvalue_solve(k_full_cur, m_full_cur)

    eigenmodes_cur = check_and_flip_sign(eigenmodes_cur)

    f1 = evaluate_residual(eigenmodes_cur['y'][mode], eigenmodes_target['y'][mode])
    f2 = evaluate_residual(eigenmodes_cur['g'][mode], eigenmodes_target['g'][mode])
    f3 = evaluate_residual([eigenfrequencies_cur[mode]], [eigenfrequencies_target[mode]])

    # deformation and frequency relatively more important, than rotation
    # weights = [0.4, 0.2, 0.4]

    # deformation, rotation, frequency relatiive similar importance
    weights = [0.33, 0.33, 0.33]

    gamma = 2
    components = [weights[0]*f1**gamma, weights[1]*f2**gamma, weights[2]*f3**gamma]
    f = sum(components)

    print('Sys params: ', ', '.join([str(val) for val in sys_params]))
    print('Components: ', ', '.join([str(val) for val in components]))
    print('Objective funtion: ', str(f))

    return f

optimizable_function = partial(objective_function, consider_mode, eigenmodes_target, eigenfrequencies_target)

# defining bounds
bnds = ((0.001, 100),(0.001, 100),(0.001, 100))

# alternatively inequality constraints
cnstrts = [{'type': 'ineq', 'fun': lambda x: 100 - x[0]},
           {'type': 'ineq', 'fun': lambda x: 100 - x[1]},
           {'type': 'ineq', 'fun': lambda x: 100 - x[2]},
           {'type': 'ineq', 'fun': lambda x: x[0] - 0.001},
           {'type': 'ineq', 'fun': lambda x: x[1] - 0.001},
           {'type': 'ineq', 'fun': lambda x: x[2] - 0.001}]

# SLSQP works with bounds
res_scalar = minimize(optimizable_function, [0.12, 0.15, 0.17], method='SLSQP', bounds=bnds, tol=1e-3, options={'gtol': 1e-3, 'ftol': 1e-3, 'disp': True})

# trust-constr runs, but not ideal
# res_scalar = minimize(optimizable_function, [0.12, 0.15, 0.17], method='trust-constr', bounds=bnds, tol=1e-3, options={'gtol': 1e-3, 'disp': True})

# Nelder-Mead, BFGS does not work
# res_scalar = minimize(optimizable_function, [0.12, 0.15, 0.17], method='Nelder-Mead', tol=1e-3, options={'gtol': 1e-3, 'ftol': 1e-3, 'disp': True})
# res_scalar = minimize(optimizable_function, [0.12, 0.15, 0.17], method='BFGS', tol=1e-3, options={'gtol': 1e-3, 'ftol': 1e-3, 'disp': True})

# TNC, L-BFGS-B, Powell does not work
# res_scalar = minimize(optimizable_function, [0.12, 0.15, 0.17], method='L-BFGS-B', bounds=bnds, tol=1e-2, options={'gtol': 1e-3, 'disp': True})

# COBYLA does not work
# res_scalar = minimize(optimizable_function, [0.12, 0.15, 0.17], method='COBYLA', constraints=cnstrts, tol=1e-2, options={'gtol': 1e-3, 'disp': True})

# SLSQP works with constraints as well
# res_scalar = minimize(optimizable_function, [0.12, 0.15, 0.17], method='SLSQP', constraints=cnstrts, tol=1e-3, options={'gtol': 1e-3, 'ftol': 1e-3, 'disp': True})

# trust-constr does not work
# res_scalar = minimize(optimizable_function, [0.12, 0.15, 0.17], method='trust-constr', constraints=cnstrts, tol=1e-3, options={'gtol': 1e-3, 'disp': True})

print(res_scalar.x)

k_full_final, m_full_final = build_system_matrix(res_scalar.x)
eigenmodes_final, eigenfrequencies_final = eigenvalue_solve(k_full_final, m_full_final)
eigenmodes_final = check_and_flip_sign(eigenmodes_final)
static_deformation_final = static_solve(k_full_final)

print('Eigenvalue solve results')
for idx, freq in enumerate(eigenfrequencies_target):

    print('   Mode: ', str(idx+1))
    print('      Eigenfeq: ')
    print('         Target: ', str(eigenfrequencies_target[idx]))
    print('         Final: ', str(eigenfrequencies_final[idx]))

    print('         Res y: ', str(evaluate_residual(eigenmodes_final['y'][idx], eigenmodes_target['y'][idx])))
    print('         Res g: ', str(evaluate_residual(eigenmodes_final['g'][idx], eigenmodes_target['g'][idx])))

    if idx == consider_mode:
        fig, axs = plt.subplots(2)
        fig.suptitle('Eigenvalue solve')
        axs[0].plot(x, eigenmodes_target['y'][idx],'k--', label='target')
        axs[0].plot(x, eigenmodes_final['y'][idx],'r-', label='final')
        axs[1].plot(x, eigenmodes_target['g'][idx], 'k--', label='target')
        axs[1].plot(x, eigenmodes_final['g'][idx], 'r-', label='final')
        plt.legend()
        plt.show()


print('Static solve results')

print('Res y: ', str(evaluate_residual(static_deformation_final['y'],static_deformation_target['y'])))
print('Res g: ', str(evaluate_residual(static_deformation_final['g'],static_deformation_target['g'])))

fig, axs = plt.subplots(2)
fig.suptitle('Static solve')
axs[0].plot(x, static_deformation_target['y'],'k--', label='target')
axs[0].plot(x, static_deformation_final['y'],'r-', label='final')
axs[1].plot(x, static_deformation_target['g'], 'k--', label='target')
axs[1].plot(x, static_deformation_final['g'], 'r-', label='final')
plt.legend()
plt.show()