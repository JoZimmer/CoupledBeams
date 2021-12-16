import numpy as np

from element_functions import get_stiffness_matrix_var, get_stiffness_matrix_tar, get_mass_matrix_var, get_mass_matrix_tar

def build_system_matrix(var_handles, n_el=3, n_dofs_per_node=2):

    n_nodes = n_el + 1
    nodes_per_elem = 2

    k_sys = np.zeros((n_nodes * n_dofs_per_node,
                      n_nodes * n_dofs_per_node))

    m_sys = np.zeros((n_nodes * n_dofs_per_node,
                      n_nodes * n_dofs_per_node))

    for el_id in range(0,n_el,1):

        k_el = get_stiffness_matrix_var(var_handles[0])
        m_el = get_mass_matrix_var(var_handles[1],var_handles[2])

        start = n_dofs_per_node * el_id
        end = start + n_dofs_per_node * nodes_per_elem

        k_sys[start: end, start: end] += k_el
        m_sys[start: end, start: end] += m_el

    return k_sys, m_sys

def apply_bc_by_reduction(matrix, dofs_of_bc=[0,1], n_el=3, n_dofs_per_node=2, axis='both'):

    n_nodes = n_el + 1
    nodes_per_elem = 2

    n_dofs_total = np.arange(n_nodes * n_dofs_per_node)
    dofs_to_keep = list(set(n_dofs_total) - set(dofs_of_bc))

    if axis == 'both':
        ixgrid = np.ix_(dofs_to_keep, dofs_to_keep)
    # for a force vector
    elif axis == 'row_vector':
        ixgrid = np.ix_(dofs_to_keep, [0])
        matrix = matrix.reshape([len(matrix), 1])

    return matrix[ixgrid]

def recuperate_bc_by_extension(matrix, dofs_of_bc=[0,1], n_el=3, n_dofs_per_node=2, axis='both'):

    n_nodes = n_el + 1
    nodes_per_elem = 2

    n_dofs_total = np.arange(n_nodes * n_dofs_per_node)
    dofs_to_keep = list(set(n_dofs_total) - set(dofs_of_bc))

    if axis == 'both':
        rows = len(n_dofs_total)
        cols = rows
        ixgrid = np.ix_(dofs_to_keep,dofs_to_keep)
        extended_matrix = np.zeros((rows,cols))
    elif axis == 'row':
        rows = len(n_dofs_total)
        cols = matrix.shape[1]
        # make a grid of indices on interest
        ixgrid = np.ix_(dofs_to_keep, np.arange(matrix.shape[1]))
        extended_matrix = np.zeros((rows, cols))
    elif axis == 'column':
        rows = matrix.shape[0]
        cols = len(n_dofs_total)
        # make a grid of indices on interest
        ixgrid = np.ix_(np.arange(matrix.shape[0]), dofs_to_keep)
        extended_matrix = np.zeros((rows, cols))
    elif axis == 'row_vector':
        rows = len(n_dofs_total)
        cols = 1
        ixgrid = np.ix_(dofs_to_keep, [0])
        matrix = matrix.reshape([len(matrix), 1])
        extended_matrix = np.zeros((rows, cols))
    elif axis == 'column_vector':
        rows = len(n_dofs_total)
        cols = 1
        ixgrid = np.ix_(dofs_to_keep)
        extended_matrix = np.zeros((rows,))

    extended_matrix[ixgrid] = matrix

    return extended_matrix

if __name__ == '__main__':
    print('Full system stiffness matrix: ')
    k_full, m_full = build_system_matrix([1.0,1.0,1.0])
    print(k_full)
    print(np.shape(k_full))

    print('Reduced system stiffness matrix: ')
    k_red = apply_bc_by_reduction(k_full)
    print(k_red)
    print(np.shape(k_red))
    print(np.array_equal(k_full[2:,2:],k_red))

    print('Extended (padded by zeros) stiffness matrix: ')
    k_ext = recuperate_bc_by_extension(k_red)
    print(k_ext)
    print(np.array_equal(np.array(np.shape(k_full)),np.array(np.shape(k_ext))))


    print('Full system mass matrix: ')
    print(m_full)
    print(np.shape(m_full))

    print('Reduced system stiffness matrix: ')
    m_red = apply_bc_by_reduction(m_full)
    print(m_red)
    print(np.shape(m_red))
    print(np.array_equal(m_full[2:,2:],m_red))

    print('Extended (padded by zeros) mass matrix: ')
    m_ext = recuperate_bc_by_extension(m_red)
    print(m_ext)
    print(np.array_equal(np.array(np.shape(m_full)),np.array(np.shape(m_ext))))