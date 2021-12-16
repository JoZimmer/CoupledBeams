import matplotlib.pyplot as plt
import numpy as np

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


print('Eigenvalue solve results')
for idx, freq in enumerate(eigenfrequencies_target):
    print('   Mode: ', str(idx+1))
    print('      Eigenfeq: ')
    print('         Target: ', str(eigenfrequencies_target[idx]))
    print('         Initial: ', str(eigenfrequencies_initial[idx]))

    print('         Res y: ', str(evaluate_residual(eigenmodes_initial['y'][idx], eigenmodes_target['y'][idx])))
    print('         Res g: ', str(evaluate_residual(eigenmodes_initial['g'][idx], eigenmodes_target['g'][idx])))

    if idx == 0:
        fig, axs = plt.subplots(2)
        fig.suptitle('Eigenvalue solve')
        axs[0].plot(x, eigenmodes_target['y'][idx],'k--', label='target')
        axs[0].plot(x, eigenmodes_initial['y'][idx],'r-', label='initial')
        axs[1].plot(x, eigenmodes_target['g'][idx], 'k--', label='target')
        axs[1].plot(x, eigenmodes_initial['g'][idx], 'r-', label='initial')
        plt.legend()
        plt.show()


print('Static solve results')

print('Res y: ', str(evaluate_residual(static_deformation_initial['y'],static_deformation_target['y'])))
print('Res g: ', str(evaluate_residual(static_deformation_initial['g'],static_deformation_target['g'])))

fig, axs = plt.subplots(2)
fig.suptitle('Static solve')
axs[0].plot(x, static_deformation_target['y'],'k--', label='target')
axs[0].plot(x, static_deformation_initial['y'],'r-', label='initial')
axs[1].plot(x, static_deformation_target['g'], 'k--', label='target')
axs[1].plot(x, static_deformation_initial['g'], 'r-', label='initial')
plt.legend()
plt.show()