from scipy import linalg
import numpy as np

from system_functions import apply_bc_by_reduction, recuperate_bc_by_extension

def eigenvalue_solve(k,m, n_el=3, dof_labels=['y','g']):

    n_nodes = n_el + 1

    eig_values_raw, eigen_modes_raw = linalg.eigh(apply_bc_by_reduction(k),
                                                  apply_bc_by_reduction(m))

    eig_values = np.sqrt(np.real(eig_values_raw))

    eigenfrequencies = eig_values / 2. / np.pi #rad/s
    eig_periods = 1 / eigenfrequencies

    gen_mass = np.matmul(np.matmul(np.transpose(eigen_modes_raw), apply_bc_by_reduction(m)), eigen_modes_raw)

    is_identiy = np.allclose(gen_mass, np.eye(gen_mass.shape[0]))
    if not is_identiy:
        raise Exception('generalized mass is not identiy')

    eigenmodes = {}
    for dof in dof_labels:
        eigenmodes[dof] = []

    for i in range(len(eigenfrequencies)):
        for j, dof in enumerate(dof_labels):
            eigenmodes[dof].append(np.zeros(n_nodes))
            eigenmodes[dof][i][1:] = eigen_modes_raw[j:,i][::len(dof_labels)]

    return eigenmodes, eigenfrequencies


def static_solve(k, force_vector=np.array([750.0, 25.0, 1250.0, 35.0, 1500.0, 45.0]), n_el=3, dof_labels=['y','g']):
    n_nodes = n_el + 1

    static_result = np.linalg.solve(apply_bc_by_reduction(k), force_vector)

    static_deformation = {}
    for dof in dof_labels:
        static_deformation[dof] = np.zeros(n_nodes)

    for i, label in enumerate(dof_labels):
        static_deformation[label][1:] = static_result[i::len(dof_labels)]

    return static_deformation


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    from system_functions import build_system_matrix


    k_full, m_full = build_system_matrix([1.0,1.0,1.0])

    x = [0.0, 1.0, 2.0, 3.0]

    eigenmodes, eigenfrequencies = eigenvalue_solve(k_full, m_full)
    print('Eigenvalue solve results')
    for idx, freq in enumerate(eigenfrequencies):
        print('   Mode: ', str(idx+1))
        print('      Eigenfeq: ', str(freq))
        print('      y-entries: ', eigenmodes['y'][idx])
        print('      g-entries: ', eigenmodes['g'][idx])

        if idx == 0:
            fig, axs = plt.subplots(2)
            fig.suptitle('Eigenvalue solve')
            axs[0].plot(x, eigenmodes['y'][idx])
            axs[1].plot(x, eigenmodes['g'][idx])
            plt.show()

    static_deformation = static_solve(k_full)
    print('Static solve results')
    print('   y-entries: ', static_deformation['y'])
    print('   g-entries: ', static_deformation['g'])

    fig, axs = plt.subplots(2)
    fig.suptitle('Static solve')
    axs[0].plot(x, static_deformation['y'])
    axs[1].plot(x, static_deformation['g'])
    plt.show()