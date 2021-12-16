import numpy as np


def get_stiffness_matrix_var(alpha, L=15.0, E=8750.0,Iy=17500.0):

    k_yy_11 = 12.0 * E * Iy / L**3
    k_yy_12 = -k_yy_11

    k_gg_11 = 4.0 * E * Iy / L
    k_gg_12 = 2.0 * E * Iy / L

    k_yg = 6.0 * E * Iy / L**2
    akyg = alpha * k_yg

    k = np.array([[ k_yy_11, -akyg,    k_yy_12, -akyg],
                    [-akyg,     k_gg_11, akyg,     k_gg_12],
                    [ k_yy_12,  akyg,    k_yy_11,  akyg],
                    [-akyg,     k_gg_12, akyg,     k_gg_11]])
    return k


def get_stiffness_matrix_tar():
    return get_stiffness_matrix_var(alpha=1.0)


def get_mass_matrix_var(beta1, beta2, rho=150.0, L=15.0):

    m_yy_11 = rho * L * 13./35.
    m_yy_12 = rho * L * 9./70.

    m_gg_11 = rho * L**3 /105.
    m_gg_12 = rho * L**3 /140.

    m_yg_11 = rho * (L**2) * 11./210.
    m_yg_12 = rho * (L**2) * 13./420.

    b1myg = beta1 * m_yg_11
    b2myg = beta2 * m_yg_12

    m = np.array([[ m_yy_11, -b1myg,    m_yy_12,  b2myg],
                  [-b1myg,    m_gg_11, -m_yg_12, -m_gg_12],
                  [ m_yy_12, -m_yg_12,  m_yy_11,  b1myg],
                  [ b2myg,   -m_gg_12,  b1myg,    m_gg_11 ]])
    return m


def get_mass_matrix_tar():
    return get_mass_matrix_var(beta1=1.0,beta2=1.0)


if __name__ == '__main__':
    print('Stiffness target: ')
    print(get_stiffness_matrix_tar())
    print('Stiffness initial: ')
    print(get_stiffness_matrix_var(alpha=0.00005))

    print('Mass target: ')
    print(get_mass_matrix_tar())
    print('Mass initial: ')
    print(get_mass_matrix_var(beta1=0.00025,beta2=0.00075))