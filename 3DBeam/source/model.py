from copy import deepcopy
import numpy as np 
import matplotlib.pyplot as plt
from scipy import linalg
from scipy.optimize import minimize, minimize_scalar
from functools import partial
from source.bernoulli_element import BernoulliElement
from source.utilities import utilities as utils
from source.utilities import global_definitions as GD

num_zero = 1e-15

class BeamModel(object):

    def __init__(self, parameters, adjust_mass_density_for_total=False, optimize_frequencies_init = True, apply_k_geo = False):

 
        # MATERIAL; GEOMETRIC AND ELEMENT INFORMATION
        self.parameters = parameters
        self.nacelle_mass = self.parameters['nacelle_mass']
        self.dim = self.parameters['dimension']
        self.n_dofs_node = GD.DOFS_PER_NODE[self.dim]
        self.dof_labels = GD.DOF_LABELS[self.dim]
        self.dofs_of_bc = self.parameters['dofs_of_bc']
        if self.parameters['type_of_bc'] == 'spring':
            self.dofs_of_bc = [0] # nur in längsrichtung einspannen 
        
        self.n_elems = parameters['n_elements']
        self.n_nodes = self.n_elems + 1
        self.nodes_per_elem = self.parameters['nodes_per_elem']

        self.parameters['intervals'] = []
        for idx, val in enumerate(parameters["defined_on_intervals"]):
            self.parameters["intervals"].append({
                'bounds': val['interval_bounds'],
                # further quantities defined by polynomial coefficient as a function of running coord x
                'c_a': val["area"],
                'c_iz': val["Iz"],
                'c_d':val["D"]
            })

        self.initialize_elements() 

        self.apply_k_geo = apply_k_geo
        self.build_system_matricies()
        self.calculate_total_mass_and_volume()
        self.calculate_eigengewicht_knoten_kraft()
        
        self.eigenvalue_solve()

        if adjust_mass_density_for_total:
            from source.optimizations import Optimizations

            mass_opt = Optimizations(self)

            mass_opt.adjust_density_for_target_total_mass(target_total_mass=self.parameters['total_mass_tower'],
                                                          print_to_console=True)

        if optimize_frequencies_init:
            from source.optimizations import Optimizations

            self.init_opt = Optimizations(self)

            target_freqs = self.parameters['eigen_freqs_target']

            self.init_opt.adjust_sway_z_stiffness_for_target_eigenfreq(target_freqs[0], 
                                                                        target_mode = 0,
                                                                        print_to_console=True)

            self.init_opt.adjust_mass_density_for_target_eigenfreq(target_freqs[0], 
                                                                    target_mode = 0,
                                                                    print_to_console=True)

            if self.dim == '3D':
                self.init_opt.adjust_sway_y_stiffness_for_target_eigenfreq(target_freqs[1], 
                                                                        target_mode = 1,
                                                                        print_to_console=False)


                self.init_opt.adjust_torsional_stiffness_for_target_eigenfreq(target_freqs[2],
                                                                            target_mode = 2,
                                                                            print_to_console=False)


# # ELEMENT INITIALIZATION AND GLOBAL MATRIX ASSAMBLAGE

    def initialize_elements(self):
        '''
        eine Listen von Bernoulli elementen erstellen -> Steifigkeitsmatrizen werden in seperater Funktion erstellt
        '''
        self.nodal_coordinates = {}
        self.elements = []

        lx_i = self.parameters['lx_total_beam'] / self.n_elems
        self.nodal_coordinates['x0'] = np.zeros(self.n_nodes)
        self.nodal_coordinates['y0'] = np.zeros(self.n_nodes)
        for i in range(self.n_nodes):
            self.nodal_coordinates['x0'][i] = i * lx_i

        if not self.parameters['intervals']:
            for i in range(self.n_elems):
                e = BernoulliElement(self.parameters, elem_length = lx_i , elem_id = i)
                self.elements.append(e)
        else:
            # MIT INTERVAL WEISE DEFINIERTEN CROSS SECTION PARAMETERN
            self.relative_length_ratios = [1.0] * self.n_elems
            param_elem_length_sum = sum(self.relative_length_ratios)
            param_elem_length_cumul = [0.0]
            for idx, el_length in enumerate(self.relative_length_ratios):
                param_elem_length_cumul.append(sum(self.relative_length_ratios[:idx+1]))
            param_elem_length_cumul_norm = [
                x/param_elem_length_sum for x in param_elem_length_cumul]

            self.parameters['x'] = [x * self.parameters['lx_total_beam'] for x in param_elem_length_cumul_norm]
            self.parameters['x_mid'] = [(a+b)/2 for a, b in zip(self.parameters['x'][:-1], self.parameters['x'][1:])]            

            for i in range(self.n_elems):
                element_params = self.initialize_element_geometric_parameters(i)
                self.parameters.update(element_params)
                coord = np.array([[self.parameters['x'][i], 0.0, 0.0],
                                [self.parameters['x'][i + 1], 0.0, 0.0]])
                e = BernoulliElement(self.parameters, elem_length = lx_i , elem_id = i)
                self.elements.append(e)

    def build_system_matricies(self, params_yg = [1.0,1.0,1.0], EIz_param = 1.0, params_k_ya = [0.0,0.0], params_m_ya = [0.0,0.0,0.0]):
        ''' 
        assigns K_el and M_el to the object. BC is applied already
        The coupling prameters are introduced
        params_yg:
            0: alpha - stiffness coupling yg
            1: beta1 - mass coupling yg1
            2: beta2 - mass coupling yg2
        EIz_param:
            Faktor für EIz 
        params_k_ya:
            0: omega - stiffness coupling ya
            1: omega1 - stiffness coupling gamma - a
        params_m_ya:
            1: psi1: - mass couling ya_11
            2: psi2: - mass coupling ya_12
            3: psi3: - mass coupling ga_11/12 
        '''
        self.k = np.zeros((self.n_nodes * self.n_dofs_node,
                            self.n_nodes * self.n_dofs_node))

        self.m = np.zeros((self.n_nodes * self.n_dofs_node,
                            self.n_nodes * self.n_dofs_node))
        
        for element in self.elements:

            k_el = element.get_stiffness_matrix_var()
            m_el = element.get_mass_matrix_var()
            
            if self.apply_k_geo:
                # TODO pretension sollte auch als Vektor über die Höhe gegeben werdne können -> Staffelung
                pretension = self.parameters['vorspannkraft']
                k_el_geo = element.get_stiffness_matrix_geometry(axial_force = pretension)

                k_el += k_el_geo

            start = self.n_dofs_node * element.index
            end = start + self.n_dofs_node * self.nodes_per_elem

            # assemble to global 
            self.k[start: end, start: end] += k_el
            self.m[start: end, start: end] += m_el

        self.comp_k = self.apply_bc_by_reduction(self.k)
        self.comp_m = self.apply_bc_by_reduction(self.m)

        # Muss nach der reduction der restlichen BCs gemacht werden
        if self.parameters['type_of_bc'] == 'spring':
            self.comp_k = self.apply_spring_bc(self.comp_k)
        
        # Aus Paper: Beton Turm Strukturoptimierung - masse am letzen knoten hinzufügen
        self.comp_m[-2,-2] += self.nacelle_mass
        self.comp_m[-3,-3] += self.nacelle_mass
        
        self.b = self.get_damping()
        self.comp_b = self.apply_bc_by_reduction(self.b)
        
    def initialize_element_geometric_parameters(self, i):
        element_params = {}
        # element properties
        x = self.parameters['x_mid'][i]
        # area
        element_params['A'] = self.evaluate_characteristic_on_interval(x, 'c_a')

        # second moment of inertia
        element_params['Iz'] = self.evaluate_characteristic_on_interval(x, 'c_iz')

        # Durchmesser
        element_params['D'] = self.evaluate_characteristic_on_interval(x, 'c_d')
       
        return element_params

    def evaluate_characteristic_on_interval(self, running_coord, characteristic_identifier):
        '''
        NOTE: continous polynomial defined within interval
        starting from the local coordinate 0.0
        so a shift is needed
        see shifted coordinate defined as: running_coord-val['bounds'][0]

        TODO: might not be robust enough with end-check -> add some tolerance
        '''
        from numpy.polynomial import Polynomial
        for val in self.parameters['intervals']:
            if "End" not in val['bounds']:
                if val['bounds'][0] <= running_coord < val['bounds'][1]:
                    polynom = Polynomial(val[characteristic_identifier])
                    polynom_x = polynom (running_coord - val['bounds'][0])
                    return polynom_x

                    # Alternative, da die werte schon als sectional means kommen
                    # return val[characteristic_identifier][0]

            elif "End" in val['bounds']:
                if val['bounds'][0] <= running_coord <= self.parameters['lx']:
                    polynom = Polynomial(val[characteristic_identifier])
                    polynom_x = polynom(running_coord - val['bounds'][0])
                    return polynom_x

                    # Alternative, da die werte schon als sectional means kommen
                    # return val[characteristic_identifier][0]

    def update_equivalent_nodal_mass(self):
        '''
        The mass matrix for the beam element can be
        lumped or consistent

        For the equivalent nodal mass a lumped distribution
        is assumned

        Here only element and NO point mass contribution
        '''
        self.parameters['nodal_mass'] = [0 for val in self.parameters['x']]
        self.volume = 0

        for idx, element in enumerate(self.elements):
            self.parameters['nodal_mass'][idx] += 0.5 * element.A * element.rho * element.L
            self.parameters['nodal_mass'][idx+1] += 0.5 * element.A * element.rho * element.L
            self.volume += element.V

    def calculate_total_mass_and_volume(self, print_to_console=False):
        # update influencing parameters
        self.update_equivalent_nodal_mass()

        self.parameters['m_tot'] = sum(self.parameters['nodal_mass'])
        
    def calculate_eigengewicht_knoten_kraft(self):
        '''
        'nodal_mass' muss mit self.update_equivalent_nodal_mass() schon ausgeführt worden sein
        Gewichtskraft in [N]
        Negativ da Druckkraft
        '''
        self.eigengewicht = np.zeros(self.n_nodes) 
        volume_total = 0
        for i, element in enumerate(self.elements):
            
            self.eigengewicht[i] = element.V * element.rho * GD.GRAVITY
            volume_total += element.V
            #self.eigengewicht[-2-i] += self.eigengewicht[-1-i] + element.V * element.rho * GD.GRAVITY

        # In Druckkraft umwandeln
        self.eigengewicht *= -1


# # BOUNDARY CONDITIONS

    def apply_bc_by_reduction(self, matrix, axis = 'both'):
        n_dofs_total = np.arange(self.n_nodes * self.n_dofs_node)
        dofs_to_keep = list(set(n_dofs_total) - set(self.dofs_of_bc))
        
        if axis == 'both':
            ixgrid = np.ix_(dofs_to_keep, dofs_to_keep)
        # for a force vector
        elif axis == 'row_vector':
            ixgrid = np.ix_(dofs_to_keep, [0])
            matrix = matrix.reshape([len(matrix), 1])
        elif axis == 'row':
            # make a grid of indices on interest
            ixgrid = np.ix_(dofs_to_keep, np.arange(matrix.shape[1]))

        return matrix[ixgrid]

    def recuperate_bc_by_extension(self, matrix, axis = 'both'):
        n_dofs_total = np.arange(self.n_nodes * self.n_dofs_node)
        dofs_to_keep = list(set(n_dofs_total) - set(self.dofs_of_bc))

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

    def apply_spring_bc(self, k):
        '''
        den einträgen werden die in den parameters['spring_stiffness'] gegebenen Steifigkeiten zugeordnet
        bisher nur 2D u und gamma als spring
        '''
        if self.dim == '3D':
            raise Exception('Spring stiffness Boundary condition not impelemnted for 3D')
         
        k[0,0] = self.parameters['spring_stiffness'][0]
        k[1,1] = self.parameters['spring_stiffness'][1]

        return k

# # SOLVES

    def eigenvalue_solve(self):
        '''
        fills: 
        self.eigenfrequencies as a list 
        self.eigenmodes as a dictionary with the dof as key and list as values 
        list[i] -> mode id 
        '''
        try:
            self.eigen_values_raw, self.eigen_modes_raw = linalg.eigh(self.comp_k, self.comp_m)
        except np.linalg.LinAlgError:
            # Not handled since this error migh abort an optimization
            return 0

        self.eig_values = np.sqrt(np.real(self.eigen_values_raw))
        
        self.eigenfrequencies = self.eig_values / 2. / np.pi #rad/s
        self.eig_periods = 1 / self.eigenfrequencies

        # generalized mass stuff -> nicht so wichtig 
        gen_mass = np.matmul(np.matmul(np.transpose(self.eigen_modes_raw), self.comp_m), self.eigen_modes_raw)

        is_identiy = np.allclose(gen_mass, np.eye(gen_mass.shape[0]), rtol=1e-05)
        if not is_identiy:
            print ('generalized mass matrix:')
            for i in range(gen_mass.shape[0]):
                for j in range(gen_mass.shape[1]):
                    if gen_mass[i][j] < 1e-3:
                        gen_mass[i][j] = 0
            print (gen_mass)
            raise Exception('generalized mass is not identiy')
    
        if not is_identiy:
            return 0

        # eigenmodes nach Richtung sortieren
        self.eigenmodes = {}
        for dof in self.dof_labels:
            self.eigenmodes[dof] = []
        # NOTE: her only fixed free boundary implemented 
        # could also be done with apply BC and recuperate BC to make it generic   
        for i in range(len(self.eigenfrequencies)):
            for j, dof in enumerate(self.dof_labels):
                self.eigenmodes[dof].append(np.zeros(self.n_nodes))
                dof_and_mode_specific = self.eigen_modes_raw[j:,i][::len(self.dof_labels)]
                start = self.n_nodes - dof_and_mode_specific.shape[0]
                self.eigenmodes[dof][i][start:] = self.eigen_modes_raw[j:,i][::len(self.dof_labels)]

        self.eigenmodes = utils.check_and_flip_sign_dict(self.eigenmodes)

        self.eig_freqs_sorted_indices = np.argsort(self.eigenfrequencies)
        
    def static_analysis_solve(self, load_vector_file = None, apply_mean_dynamic = False, directions = 'y', 
                                    return_result = False, add_eigengewicht=True, add_imperfektion = True):
        ''' 
        pass_load_vector_file (directions ist dann unrelevant)
        if apply_mean_dynamic: the dynamic load file is taken and at each point the mean magnitude
        direction: if 'all' all load directions are applied
        ''' 
        if load_vector_file != None:
            load_vector = np.load(load_vector_file)
        elif apply_mean_dynamic:
            dyn_load = np.load(self.parameters['dynamic_load_file'])
            if directions == 'all':
                load_vector = np.apply_along_axis(np.mean, 1, dyn_load)
            else:
                for direction in directions:
                    dir_id = GD.DOF_LABELS['3D'].index(direction)
                    mean_load = np.apply_along_axis(np.mean, 1, dyn_load)[dir_id::GD.n_dofs_node['3D']]
                    load_vector[dir_id::GD.n_dofs_node['3D']] = mean_load
        else:
            raise Exception('No load vector file provided')

        if add_eigengewicht:
            load_vector[0::GD.DOFS_PER_NODE[self.dim]] += self.eigengewicht.reshape(self.n_nodes,1)
        if add_imperfektion:
            if add_eigengewicht == False:
                print ('WARINING! Einfluss aus Imperfektion ohne Eigengewicht berechnet!!')
            theta_x = self.parameters['imperfektion']
            # horizontale ersatzkraft = Schiefstellung an höhe x * Vertikalkraft an dieser stelle?!
            Fv_id = GD.DOF_LABELS[self.dim].index('x')
            Fh_id = GD.DOF_LABELS[self.dim].index('y')
            step = GD.DOFS_PER_NODE[self.dim]
            # -1 da Vertikalkraft nach unten wirkt die äquivalente Horizontalkraft aber in richtung der bereits wirkenden horizontalkraft sein sollte
            self.Fh_imp_x = theta_x * load_vector[Fv_id::step] * -1 
            load_vector[Fh_id::step] += self.Fh_imp_x

        self.static_deformation = {}
        self.reaction = {}

        load = self.apply_bc_by_reduction(load_vector, axis='row_vector')

        #missing ground node -> get bc by extension
        static_result = np.linalg.solve(self.comp_k, load)
        static_result = self.recuperate_bc_by_extension(static_result, 'row_vector')
        f = np.dot(self.k, static_result)
        self.resisting_force = load_vector - f

        for i, label in enumerate(self.dof_labels):
            self.static_deformation[label] = static_result[i::self.n_dofs_node]
            self.reaction[label] = self.resisting_force[i::self.n_dofs_node][:, 0]

        # internal forces - Schnittgrößen müssen an jedem Element gemacht werden 
        internal_forces_list = []
        self.internal_forces = {}
        for i, element in enumerate(self.elements):
            # schnittgröße am ground node ist gleich der resisting force
            if i == 0:
                internal_forces_list.append(self.resisting_force[:self.n_dofs_node][:,0])

            start = i*self.n_dofs_node
            stop = start + 2*self.n_dofs_node

            u_i = static_result[start:stop]
            f_i = np.dot(element.k_element, u_i)
            
            # negatives Schnittufer sind die hinteren einträge diese müssen mit den äußeren Lasten im Gleichgewicht stehen
            internal_forces_list.append(f_i[self.n_dofs_node:][:,0])
            # if i == self.n_elems-1:
            #     internal_forces_list.append(f_i[self.n_dofs_node:self.n_dofs_node*2][:,0])#*-1
        
        internal_forces_arr = np.array(internal_forces_list)
        for i, label in enumerate(self.dof_labels):
            self.internal_forces[label] = internal_forces_arr[:,i]

        if return_result:
            return self.internal_forces

# # RETURN OPTIMIZED PARAMETERS

    def update_optimized_parameters(self, optimization_results = None):
        ''' 
        NOTE: as it is implemented only a homogenous cross sections are updated 
        ''' 

        if self.optimize_frequencies_init:
            print ('\nupdating:')
            for param, val in self.init_opt.opt_geometric_props.items():
                if param == 'Ip':
                    init = self.elements[0].Ip
                else:
                    init = self.parameters[param]

                print ('    ',param, 'from', round(init), 'to',round(val[0]), 'increasing by:', utils.increasing_by(init, val[0]), '%')
                self.parameters[param] = val[0]
                
        if optimization_results:
            for param, val in optimization_results.items():
                init = self.parameters[param]
                print ('    ',param, 'from', init, 'to',val, '(k_ya, kga)')#, 'increasing by:', utils.increasing_by(init, val[0]), '%')
                self.parameters[param] = list(val)

        print('\nupdated optimized parameters')
        return self

# # EXTRA FUNCTIONS

    def get_transformation_matrix(self,ey,ez):
        ''' 
        transform the element matrix from the shear center to the coordinate center (geometric center)
            ey: eccentricity in y direction -> coupling z displacement with torsion
            ez: eccentricity in z direction -> coupling y displacement with torsion
        NOTE not used 
        '''
         
        # eigentlich auch noch ez, ey aus stockwerk oben drüber --> node 2
        T = np.identity(GD.n_dofs_node[self.dim]*2)
        T_transform_part_node1 =  np.array([[1., 0., ez],
                                            [0., 1., -ey],
                                            [0., 0., 1.]])
        T_transform_part_node2 =  np.array([[1., 0., ez],
                                            [0., 1., -ey],
                                            [0., 0., 1.]])
        
        start = GD.dof_lables[self.dim].index('y')
        end = GD.dof_lables[self.dim].index('a')+1
        T[start:end, start:end] = T_transform_part_node1

        start = GD.n_dofs_node[self.dim] + GD.dof_lables[self.dim].index('y')
        end = GD.n_dofs_node[self.dim] + GD.dof_lables[self.dim].index('a')+1
        T[start:end, start:end] = T_transform_part_node2
                                    
        return T

    def get_damping(self):
        """
        Calculate damping b based upon the Rayleigh assumption
        using the first 2 eigemodes - here generically i and j
        """

        mode_i = 0
        mode_j = 1
        zeta_i = self.parameters['damping_coeff']
        zeta_j = zeta_i

        self.eigenvalue_solve()

        self.rayleigh_coefficients = np.linalg.solve(0.5 *
                                                     np.array(
                                                         [[1 / self.eig_values[self.eig_freqs_sorted_indices[mode_i]],
                                                           self.eig_values[self.eig_freqs_sorted_indices[mode_i]]],
                                                          [1 / self.eig_values[self.eig_freqs_sorted_indices[mode_j]],
                                                           self.eig_values[self.eig_freqs_sorted_indices[mode_j]]]]),
                                                     [zeta_i, zeta_j])

        # return back the whole matrix - without BCs applied
        return self.rayleigh_coefficients[0] * self.m + self.rayleigh_coefficients[1] * self.k

