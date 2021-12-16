import numpy as np 
from source.utilities import global_definitions as GD

class BernoulliElement(object):

    def __init__(self, parameters, elem_length ,elem_id):

        self.parameters = parameters
        self.index = elem_id

        self.E = parameters['E_Modul']
        self.A = parameters['cross_section_area']
        self.rho = parameters['material_density']
        self.nu = parameters['nu']
        self.L = elem_length

        self.Iy = parameters['Iy']
        self.Iz = parameters['Iz']
        self.I_param = parameters['I_param']
        self.It = parameters['It']
        self.number_of_nodes = parameters['nodes_per_elem']

        self.dim = parameters['dimension']

        self.evaluate_torsional_inertia() # -> gives Ip
        self.evaluate_relative_importance_of_shear() # gives G, Py, Pz = 0 if Bernoulli here 

    def get_stiffness_matrix_var(self, alpha = 1.0, omega = 0.0, omega1 = 0.0):
        ''' 
        alpha: parameter for the coupling entry of y - g
        omega: parameter for the coupling entry of y - a (torsion)
        omega1: coupling parameter g - a 
        ''' 
        EA = self.E * self.A
        EIy = self.E * self.Iy
        EIz = self.E * self.Iz

        # bending around z - displacement in y roation gamma around z

        k_yy_11 = 12.0 * EIz / self.L**3
        k_yy_12 = -k_yy_11

        k_gg_11 = 4.0 * EIz / self.L
        k_gg_12 = 2.0 * EIz / self.L

        k_yg = 6.0 * EIz / self.L**2
        akyg = alpha * k_yg

        k_el_yg = np.array([[ k_yy_11, -akyg,    k_yy_12, -akyg],
                            [-akyg,     k_gg_11, akyg,     k_gg_12],
                            [ k_yy_12,  akyg,    k_yy_11,  akyg],
                            [-akyg,     k_gg_12, akyg,     k_gg_11]])

        if self.dim == '3D':
            
            # axial stiffness - along axis x - here marked as x
        
            k_xx_11 = 1.0
            k_xx_12 = -1.0
            k_el_x = EA/self.L * np.array([[k_xx_11, k_xx_12],
                                        [k_xx_12, k_xx_11]])

            # bending around y - displacement in z roation beta around y
            k_zz_11 = 12.0 * EIy / self.L**3
            k_zz_12 = -k_zz_11

            k_bb_11 = 4.0 * EIy / self.L
            k_bb_12 = 2.0 * EIy / self.L

            k_zb = 6.0 * EIy / self.L**2

            k_el_zb = np.array([[ k_zz_11, -k_zb,    k_zz_12, -k_zb],
                                [-k_zb,     k_bb_11, k_zb,     k_bb_12],
                                [ k_zz_12,  k_zb,    k_zz_11,  k_zb],
                                [-k_zb,     k_bb_12, k_zb,     k_bb_11]])

        # ===========================
        # TORSION - COUPLING
        # ===========================

            k_aa = self.G * self.It / self.L 
            k_ya = omega *  self.E * self.I_param / self.L**3 # from eccentricity approach to have reasonable scales of start e * k_yy_1112.0 *

            k_aa_11 = k_aa
            

            k_ga = omega1 * 6.0 * self.E * self.I_param  / self.L**2 

        if self.dim == '2D':
            # NOTE this is without axial stuff
            k_el = k_el_yg
        
        elif self.dim == '3D':
            k_el = np.array([[k_el_x[0][0], 0., 0., 0., 0., 0., k_el_x[0][1], 0., 0., 0., 0., 0.],
                             [0., k_el_yg[0][0], 0., k_ya, 0., k_el_yg[0][1], 0., k_el_yg[0][2], 0., -k_ya, 0., k_el_yg[0][3]],
                             [0., 0., k_el_zb[0][0], 0., k_el_zb[0][1], 0., 0., 0., k_el_zb[0][2], 0., k_el_zb[0][3],0.],
                             [0., k_ya, 0., k_aa_11, 0., -k_ga, 0., -k_ya, 0., -k_aa_11, 0., k_ga],
                             [0., 0., k_el_zb[0][1], 0., k_el_zb[1][1], 0., 0., 0., k_el_zb[1][2], 0., k_el_zb[1][3],0.],
                             [0., k_el_yg[0][1], 0., -k_ga, 0., k_el_yg[1][1], 0., k_el_yg[1][2], 0., k_ga, 0.,k_el_yg[1][3]],

                             [k_el_x[1][0], 0., 0., 0., 0., 0., k_el_x[1][1], 0., 0., 0., 0., 0.],
                             [0., k_el_yg[0][2], 0., -k_ya, 0., k_el_yg[1][2], 0., k_el_yg[2][2], 0., k_ya, 0., k_el_yg[2][3]],
                             [0., 0., k_el_zb[0][2],  0., k_el_zb[1][2], 0., 0., 0., k_el_zb[2][2], 0., k_el_zb[2][3],0.],
                             [0., -k_ya, 0., -k_aa_11, 0., k_ga, 0., k_ya,  0., k_aa_11, 0., -k_ga],
                             [0., 0., k_el_zb[0][3], 0., k_el_zb[1][3], 0., 0., 0., k_el_zb[2][3], 0., k_el_zb[3][3],0.],
                             [0., k_el_yg[0][3], 0., k_ga, 0., k_el_yg[1][3], 0., k_el_yg[2][3], 0., -k_ga, 0.,k_el_yg[3][3]]])

        return k_el

    def get_stiffness_matrix_tar(self):
        '''
        return the correct stiffness matrix with parameter alpha = 1.0
        '''
        return get_stiffness_matrix_var(alpha=1.0, omega=0.0)


    def get_mass_matrix_var(self, beta1 = 1.0, beta2 = 1.0, psi1 = 0.0, psi2 = 0.0, psi3 = 0.0):
        ''' 
        in Euler - Bernoulli Theroy the rotary part of the element mass matrix is neglected
        from Appendix A - Straight Beam Element Matrices - page 228
        https://link.springer.com/content/pdf/bbm%3A978-3-319-56493-7%2F1.pdf

        BENDING:
            beta1: coupling parameter y - g at same node
            beta2: coupling parameter y - g at opposing nodes
        TORSION:
            psi1:  coupling parameter y - a at same node
            psi2:  coupling parameter y - a at opposing nodes
            psi3:  coupling parameter g - a 
        '''

        if self.dim == '2D':
            m_yy_11 = self.rho * self.A * self.L* 13./35.
            m_yy_12 = self.rho * self.A * self.L* 9./70.

            m_gg_11 = self.rho * self.A * self.L**3 /105.
            m_gg_12 = self.rho * self.A * self.L**3 /140.

            m_yg_11 = self.rho * self.A * (self.L**2) * 11./210.
            m_yg_12 = self.rho * self.A * (self.L**2) * 13./420.

            b1myg = beta1 * m_yg_11
            b2myg = beta2 * m_yg_12

            m_el = np.array([[ m_yy_11, -b1myg,    m_yy_12,  b2myg],
                             [-b1myg,    m_gg_11, -m_yg_12, -m_gg_12],
                             [ m_yy_12, -m_yg_12,  m_yy_11,  b1myg],
                             [ b2myg,   -m_gg_12,  b1myg,    m_gg_11 ]])
            return m_el

        elif self.dim == '3D':
            m_const = self.rho * self.A * self.L

            # axial inertia - along axis x - here marked as x
            m_x = m_const / 6.0
            m_x_11 = 2.
            m_x_12 = 1.
            m_el_x = m_x * np.array([[m_x_11, m_x_12],
                                    [m_x_12, m_x_11]])

            # bending - inertia along axis y, rotations around axis z - here marked as gamma - g
            # translation
            Py = self.Py
            m_yg = m_const / 210 / (1 + Py) ** 2
            #
            m_yg_11 = 70 * Py ** 2 + 147 * Py + 78
            m_yg_12 = (35 * Py ** 2 + 77 * Py + 44) * self.L / 4
            m_yg_13 = 35 * Py ** 2 + 63 * Py + 27
            m_yg_14 = -(35 * Py ** 2 + 63 * Py + 26) * self.L / 4
            #
            m_yg_22 = (7 * Py ** 2 + 14 * Py + 8) * self.L ** 2 / 4
            m_yg_23 = - m_yg_14
            m_yg_24 = -(7 * Py ** 2 + 14 * Py + 6) * self.L ** 2 / 4
            #
            m_yg_33 = m_yg_11
            m_yg_34 = -m_yg_12
            #
            m_yg_44 = m_yg_22
            #
            m_el_yg_trans = m_yg * np.array([[m_yg_11, m_yg_12, m_yg_13, m_yg_14],
                                            [m_yg_12, m_yg_22, m_yg_23, m_yg_24],
                                            [m_yg_13, m_yg_23, m_yg_33, m_yg_34],
                                            [m_yg_14, m_yg_24, m_yg_34, m_yg_44]])

            # rotation
            m_yg = self.rho * self.Iz / 30 / (1 + Py) ** 2 / self.L
            #
            m_yg_11 = 36
            m_yg_12 = -(15 * Py - 3) * self.L
            m_yg_13 = -m_yg_11
            m_yg_14 = m_yg_12
            #
            m_yg_22 = (10 * Py ** 2 + 5 * Py + 4) * self.L ** 2
            m_yg_23 = - m_yg_12
            m_yg_24 = (5 * Py ** 2 - 5 * Py - 1) * self.L ** 2
            #
            m_yg_33 = m_yg_11
            m_yg_34 = - m_yg_12
            #
            m_yg_44 = m_yg_22
            #
            m_el_yg_rot = m_yg * np.array([[m_yg_11, m_yg_12, m_yg_13, m_yg_14],
                                        [m_yg_12, m_yg_22, m_yg_23, m_yg_24],
                                        [m_yg_13, m_yg_23, m_yg_33, m_yg_34],
                                        [m_yg_14, m_yg_24, m_yg_34, m_yg_44]])

            # sum up translation and rotation
            #m_el_yg = m_el_yg_trans + m_el_yg_rot
            m_el_yg = m_el_yg_trans

            #if self.domain_size == '3D':
            # bending - inertia along axis z, rotations around axis y - here marked as beta - b
            # translation
            Pz = self.Pz
            m_zb = m_const / 210 / (1 + Pz) ** 2
            #
            m_zb_11 = 70 * Pz ** 2 + 147 * Pz + 78
            m_zb_12 = -(35 * Pz ** 2 + 77 * Pz + 44) * self.L / 4.
            m_zb_13 = 35 * Pz ** 2 + 63 * Pz + 27
            m_zb_14 = (35 * Pz ** 2 + 63 * Pz + 26) * self.L / 4.
            #
            m_zb_22 = (7 * Pz ** 2 + 14 * Pz + 8) * self.L ** 2 / 4.
            m_zb_23 = -m_zb_14
            m_zb_24 = -(7 * Pz ** 2 + 14 * Pz + 6) * self.L ** 2 / 4.
            #
            m_zb_33 = m_zb_11
            m_zb_34 = - m_zb_12
            #
            m_zb_44 = m_zb_22
            #
            m_el_zb_trans = m_zb * np.array([[m_zb_11, m_zb_12, m_zb_13, m_zb_14],
                                                [m_zb_12, m_zb_22, m_zb_23, m_zb_24],
                                                [m_zb_13, m_zb_23, m_zb_33, m_zb_34],
                                                [m_zb_14, m_zb_24, m_zb_34, m_zb_44]])
            # rotation
            m_zb = self.rho * self.Iy / 30. / (1 + Pz) ** 2 / self.L

            m_zb_11 = 36.
            m_zb_12 = (15. * Pz - 3) * self.L
            m_zb_13 = -m_zb_11
            m_zb_14 = m_zb_12

            m_zb_22 = (10 * Pz ** 2 + 5 * Pz + 4) * self.L ** 2
            m_zb_23 = -m_zb_12
            m_zb_24 = (5 * Pz ** 2 - 5 * Pz - 1) * self.L ** 2

            m_zb_33 = m_zb_11
            m_zb_34 = -m_zb_12

            m_zb_44 = m_zb_22

            m_el_zb_rot = m_zb * np.array([[m_zb_11, m_zb_12, m_zb_13, m_zb_14],
                                            [m_zb_12, m_zb_22, m_zb_23, m_zb_24],
                                            [m_zb_13, m_zb_23, m_zb_33, m_zb_34],
                                            [m_zb_14, m_zb_24, m_zb_34, m_zb_44]])

            # sum up translation and rotation
            # m_el_zb = m_el_zb_trans + m_el_zb_rot
            m_el_zb = m_el_zb_trans

            # torsion inertia - around axis x - here marked as alpha - a
            m_a = m_const * self.Ip / self.A / 6.0
            m_a_11 = 2
            m_a_12 = 1
            m_el_a = m_a * np.array([[m_a_11, m_a_12],
                                     [m_a_12, m_a_11]])

            # ===================================================
            # TORSION COUPLING - BENDING AROUND Z (Y-DSIP; G-ROT)
            # ===================================================
            
            # coupling of displacement y with alpha a

            m_ya = m_const  / 10 
            
            m_ya_11 = 7/20 * m_ya * psi1 # see reason for that choice taken form frias but basically to have a reaonable scale 
            m_ya_12 = 3/20 * m_ya * psi2 # seems to be smaler than m_ya_11 maybe this could be taken into accoutn
           
            # coupling gamma - alpha (also from Frias with ratio of m_ya - m_ga)
            m_ga = m_const * self.L 
            m_ga_11 = 1/20 * m_ga *  psi3
            m_ga_12 = 1/30 * m_ga *  psi3
            

            # ===================================================
            # TORSION COUPLING - BENDING AROUND Y (Z-DSIP; B-ROT)
            # ===================================================

            m_za = -m_const * 0

            m_za_11 = 7/20 * m_za * psi1 
            m_za_12 = 3/20 * m_za * psi2

            m_ba = m_const * self.L

            m_ba_11 = 1/20 * m_ba * psi3
            m_ba_12 = 1/30 * m_ba * psi3


            # assemble all components
            
            m_el = np.array([[m_el_x[0][0], 0., 0., 0., 0., 0., m_el_x[0][1], 0., 0., 0., 0., 0.],
                                [0., m_el_yg[0][0], 0., m_ya_11, 0., m_el_yg[0][1], 0., m_el_yg[0][2], 0., m_ya_12, 0., m_el_yg[0][3]],
                                [0., 0., m_el_zb[0][0], -m_za_11, m_el_zb[0][1], 0., 0., 0., m_el_zb[0][2], -m_za_12, m_el_zb[0][3],0.],
                                [0., m_ya_11, -m_za_11, m_el_a[0][0], m_ba_11, m_ga_11, 0., m_ya_12, -m_za_12, m_el_a[0][1], m_ba_12, m_ga_12],
                                [0., 0., m_el_zb[0][1], m_ba_11, m_el_zb[1][1], 0., 0., 0., m_el_zb[1][2], -m_ba_12, m_el_zb[1][3], 0.],
                                [0., m_el_yg[0][1], 0., m_ga_11, 0., m_el_yg[1][1], 0., m_el_yg[1][2], 0., -m_ga_12, 0., m_el_yg[1][3]],

                                [m_el_x[1][0], 0., 0., 0., 0., 0., m_el_x[1][1], 0., 0., 0., 0., 0.],
                                [0., m_el_yg[0][2], 0., m_ya_12, 0., m_el_yg[1][2], 0., m_el_yg[2][2], 0., m_ya_11, 0., m_el_yg[2][3]],
                                [0., 0., m_el_zb[0][2], -m_za_12, m_el_zb[1][2], 0., 0., 0., m_el_zb[2][2], -m_za_11, m_el_zb[2][3], 0.],
                                [0., m_ya_12, -m_za_12, m_el_a[1][0], m_ba_12, m_ga_12, 0., m_ya_11, -m_za_11, m_el_a[1][1], -m_ba_11, -m_ga_11],
                                [0., 0., m_el_zb[0][3], -m_ba_12, m_el_zb[1][3], 0., 0., 0., m_el_zb[2][3], -m_ba_11, m_el_zb[3][3], 0.],
                                [0., m_el_yg[0][3], 0., -m_ga_12, 0., m_el_yg[1][3], 0., m_el_yg[2][3], 0., -m_ga_11, 0., m_el_yg[3][3]]])

            if self.dim == '2D':
                m_el_1 = np.array([
                            [m_el_yg[0][0], m_el_yg[0][1], m_el_yg[0][2], m_el_yg[0][3]],
                            [m_el_yg[0][1], m_el_yg[1][1], m_el_yg[1][2], m_el_yg[1][3]],

                            [m_el_yg[0][2], m_el_yg[1][2], m_el_yg[2][2], m_el_yg[2][3]],
                            [m_el_yg[0][3], m_el_yg[1][3], m_el_yg[2][3], m_el_yg[3][3]]])
            return m_el
      
    def get_mass_matrix_tar(self):
        return get_mass_matrix_var(beta1=1.0,beta2=1.0,psi1=0.0, psi2=0.0)


    def evaluate_relative_importance_of_shear(self):
        '''
        for Bernoulli Py and Pz are just 0, this is the Timoshenko element with no shear influence
        ''' 
        self.G = self.E / 2 / (1 + self.nu)
        # relative importance of the shear deformation to the bending one
        self.Py = 0.
        self.Pz = 0.

    def evaluate_torsional_inertia(self):
        # polar moment of inertia
        # assuming equivalency with circle
        self.Ip = self.Iy + self.Iz

    