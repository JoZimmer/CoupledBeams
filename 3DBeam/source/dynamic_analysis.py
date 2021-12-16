import numpy as np 
import matplotlib.pyplot as plt 

from source.utilities import global_definitions as GD 
from source.utilities import utilities as utils

''' 
This class is entirely copied from ParOptBeam
''' 

def transform_into_modal_coordinates(modal_transform_matrix, matrix, modes_considered):
    modal_transform_matrix_red = modal_transform_matrix[:,:modes_considered]
    matrix_transformed = np.matmul(np.matmul(np.transpose(modal_transform_matrix_red), matrix), modal_transform_matrix_red)
    matrix_transformed = matrix_transformed.round(decimals=8) # rounded off 10 **-8 to zero 
    if np.count_nonzero(matrix_transformed - np.diag(np.diagonal(matrix_transformed))) == 0: 
        matrix_as_array = np.diagonal(matrix_transformed)
    else:
        raise Exception('transformed matrix non-diagonal')
    return matrix_as_array

class DynamicAnalysis(object):
    """
    Derived class for the dynamic analysis of a given structure model

    """

    # using these as default or fallback settings
    DEFAULT_SETTINGS = {
        "type": "dynamic_analysis",
        "settings": {},
        "input": {},
        "output": {}}

    def __init__(self, structure_model, parameters):

        self.structure_model = structure_model
        self.structure_model.eigenvalue_solve()
        self.parameters = parameters
        self.damping_dummy = np.zeros(self.structure_model.m.shape)
        # time parameters
        time_integration_scheme = self.parameters['settings']['time']['integration_scheme']
        start = self.parameters['settings']['time']['start']
        stop = self.parameters['settings']['time']['end']
        # TODO check if this is the correct way
        self.dt = self.parameters['settings']['time']['step']
        steps = int((self.parameters['settings']['time']['end']
                     - self.parameters['settings']['time']['start']) / self.dt) + 1
        self.array_time = np.linspace(start, stop, steps)

        force = np.load(self.parameters['input']['file_path'])

        self.force = force

        # check dimensionality
        # of time and force
        len_time_gen_array = len(self.array_time)
        len_time_force = len(self.force[0])
        if len_time_gen_array != len_time_force:
            err_msg = "The length " + \
                str(len_time_gen_array) + \
                " of the time array generated based upon parameters\n"
            err_msg += "specified in \"runs\" -> for \"type\":\"dynamic_analysis\" -> \"settings\" -> \"time\"\n"
            err_msg += "does not match the time series length " + \
                str(len_time_force) + " of the load time history\n"
            err_msg += "specified in \"runs\" -> for \"type\":\"dynamic_analysis\" -> \"input\" -> \"file_path\"!\n"
            raise Exception(err_msg)

        # of nodes-dofs
        n_dofs_model = structure_model.n_nodes * GD.n_dofs_node['3D']
        n_dofs_force = len(self.force)

        if n_dofs_model != n_dofs_force:
            err_msg = "The number of the degrees of freedom " + \
                str(n_dofs_model) + " of the structural model\n"
            err_msg += "does not match the degrees of freedom " + \
                str(n_dofs_force) + " of the load time history\n"
            err_msg += "specified in \"runs\" -> for \"type\":\"dynamic_analysis\" -> \"input\" -> \"file_path\"!\n"
            err_msg += "The structural model has:\n"
            err_msg += "   " + \
                str(structure_model.n_elems) + " number of elements\n"
            err_msg += "   " + \
                str(structure_model.n_nodes) + " number of nodes\n"
            err_msg += "   " + str(n_dofs_model) + " number of dofs.\n"
            err_msg += "The naming of the force time history should reflect the number of nodes\n"
            err_msg += "using the convention \"dynamic_force_<n_nodes>_nodes.npy\"\n"
            digits_in_filename = [
                s for s in self.parameters['input']['file_path'].split('_') if s.isdigit()]
            if len(digits_in_filename) == 1:
                err_msg += "where currently <n_nodes> = " + \
                    digits_in_filename[0] + " (separated by underscores)!"
            else:
                err_msg += "but found multiple digits: " + \
                    ', '.join(digits_in_filename) + \
                    " (separated by underscores)!"
            raise Exception(err_msg)

        rows = len(self.structure_model.apply_bc_by_reduction(self.structure_model.k))

        # initial condition of zero displacement and velocity used for the time being.

        u0 = np.zeros(rows)  # initial displacement
        v0 = np.zeros(rows)  # initial velocity
        a0 = np.zeros(rows)  # initial acceleration
        initial_conditions = np.array([u0, v0, a0])
     
        if 'run_in_modal_coordinates' in self.parameters['settings']:
            if self.parameters['settings']['run_in_modal_coordinates']:
                self.transform_into_modal = True
                num_of_modes_specified = self.parameters['settings']['number_of_modes_considered']
                min_number_of_modes = 1
                max_number_of_modes = rows 
                if num_of_modes_specified < min_number_of_modes:
                    err_msg = "specified number of modes is less than minimum required"
                    raise Exception(err_msg)
                elif num_of_modes_specified > max_number_of_modes: 
                    err_msg = "specified number of modes is more than maximum possible"
                    raise Exception(err_msg)
                else:
                    self.num_of_modes_considered = num_of_modes_specified
                    
                u0 = np.zeros(self.num_of_modes_considered)  # initial displacement
                v0 = np.zeros(self.num_of_modes_considered)  # initial velocity
                a0 = np.zeros(self.num_of_modes_considered)  # initial acceleration
                initial_conditions = np.array([u0, v0, a0])
            else:
                self.transform_into_modal = False
                pass
        else:
            self.transform_into_modal = False
            pass
        
        # TODO check if concept of comp - computational model is robust and generic enough
        self.comp_m = np.copy(self.structure_model.comp_m)
        self.comp_k = np.copy(self.structure_model.comp_k)
        #self.comp_b = np.zeros(self.comp_k.shape)
        self.comp_b = np.copy(self.structure_model.comp_b)
        # tranformation to the modal coordinates
        if self.transform_into_modal:
            self.comp_m = transform_into_modal_coordinates(
                self.structure_model.eigen_modes_raw, self.comp_m, self.num_of_modes_considered)
            self.comp_b = transform_into_modal_coordinates(
                self.structure_model.eigen_modes_raw, self.comp_b, self.num_of_modes_considered)
            self.comp_k = transform_into_modal_coordinates(
                self.structure_model.eigen_modes_raw, self.comp_k, self.num_of_modes_considered)

        if force.shape[1] != len(self.array_time):
            err_msg = "The time step for forces does not match the time step defined"
            raise Exception(err_msg)

        # external forces
        force = self.structure_model.apply_bc_by_reduction(self.force, 'row')

        if self.transform_into_modal:
            force = np.dot(np.transpose(
                self.structure_model.eigen_modes_raw[:,:self.num_of_modes_considered]), force)

        #print(self.parameters)
        if self.parameters["settings"]["solver_type"] == "Linear":
            from source.solving_strategies.strategies.linear_solver import LinearSolver
            self.solver = LinearSolver(self.array_time, time_integration_scheme, self.dt,
                                       [self.comp_m, self.comp_b,self.comp_k],
                                       initial_conditions, force,
                                       self.structure_model)
        elif self.parameters["settings"]["solver_type"] == "Picard":
            from source.solving_strategies.strategies.residual_based_picard_solver import ResidualBasedPicardSolver
            self.solver = ResidualBasedPicardSolver(self.array_time, time_integration_scheme, self.dt,
                                                    [self.comp_m, self.comp_b,self.comp_k],
                                                    initial_conditions, force,
                                                    self.structure_model)
        elif self.parameters["settings"]["solver_type"] == "NewtonRaphson":
            from source.solving_strategies.strategies.residual_based_newton_raphson_solver import ResidualBasedNewtonRaphsonSolver
            self.solver = ResidualBasedNewtonRaphsonSolver(self.array_time, time_integration_scheme, self.dt,
                                                           [self.comp_m, self.comp_b,
                                                               self.comp_k],
                                                           initial_conditions, force,
                                                           self.structure_model)
        else:
            err_msg = "The requested solver type \"" + \
                self.parameters["settings"]["solver_type"]
            err_msg += "\" is not available \n"
            err_msg += "Choose one of: \"Linear\", \"Picard\", \"NewtonRaphson\"\n"
            raise Exception(err_msg)

    def solve(self):

        #print("Solving the structure for dynamic loads \n")
        self.solver.solve()


        # transforming back to normal coordinates :
        if self.transform_into_modal:
            self.solver.displacement = np.matmul(
                self.structure_model.eigen_modes_raw[:,:self.num_of_modes_considered], self.solver.displacement)
            self.solver.velocity = np.matmul(
                self.structure_model.eigen_modes_raw[:,:self.num_of_modes_considered], self.solver.velocity)
            self.solver.acceleration = np.matmul(
                self.structure_model.eigen_modes_raw[:,:self.num_of_modes_considered], self.solver.acceleration)

        self.solver.displacement = self.structure_model.recuperate_bc_by_extension(
            self.solver.displacement, axis='row')
        self.solver.velocity = self.structure_model.recuperate_bc_by_extension(
            self.solver.velocity, axis='row')
        self.solver.acceleration = self.structure_model.recuperate_bc_by_extension(
            self.solver.acceleration, axis='row')
        # computing the reactions
        f1 = np.dot(self.structure_model.m, self.solver.acceleration)
        f2 = np.dot(self.structure_model.b, self.solver.velocity)
        f3 = np.dot(self.structure_model.k, self.solver.displacement)
        self.solver.dynamic_reaction = self.force - f1 - f2 - f3
        #TODO : elastic support reaction computation

    def output_kinetic_energy(self, total = False):
        ''' 
        - potentielle energie eines mechanischen Systems ist durch die lage gegeben: also Feder wird gespannt -> Epot ist die Fläche unter der Last- Verformungskurve 
        - kinetische energie wie üblich 1/2 m v² 
        ''' 
        print("\nCalculate modal kinetic energy...")       
        
        m = self.structure_model.m  # which mass matrix ? lumped masses ?
        k = self.structure_model.k  
        vel = self.solver.velocity
        disp = self.solver.displacement

        self.kin_energy = np.zeros(len(self.array_time))
        self.el_energy = np.zeros(len(self.array_time))

        pre = 'KINETIC '
        for i in range(0,len(self.array_time)):
            self.kin_energy[i] = 0.5 * np.dot(np.transpose(vel[:,i]),np.dot(m,vel[:,i]))
            if total:
                pre = 'TOTAL '
                self.el_energy[i] = 0.5 * np.dot(np.transpose(disp[:,i]),np.dot(k,disp[:,i]))


        self.sum_energy = self.kin_energy + self.el_energy

        self.sum_energy_over_time = np.sum(np.multiply(self.sum_energy, self.dt/self.array_time[-1]))

        ## EXPLIZIT WAY
            # vel = {}

            # for idx, label in zip(list(range(GD.DOFS_PER_NODE[self.structure_model.domain_size])),
            #                       GD.DOF_LABELS[self.structure_model.domain_size]):
            #     start = idx
            #     step = GD.DOFS_PER_NODE[self.structure_model.domain_size]
            #     stop = self.solver.displacement.shape[0] + idx - step
            #     vel[label] = self.solver.velocity[start:stop +1:step][:]

            # # TODO: make robust for 2d and 3d, here hard coded for 3d
            # # here taking sqrt
            # vel_magn = np.power(np.power(vel['y'],2) + np.power(vel['z'],2), 0.5)

            # # here power of previous sqrt
            # kin_energy = 0.5 * np.dot(self.structure_model.parameters['m'],np.power(vel_magn,2))

            # self.sum_energy = kin_energy #+ el_energy

            # # first here introducing dt (integral) 
            # # and division with total length of time to get a normed integral
            # # could be extended to possible 
            # sum_energy_over_time = np.sum(np.multiply(self.sum_energy, self.dt/self.array_time[-1]))

        result_data = self.sum_energy

        print ('\n  ' + pre + 'ENERGY sum over time:', round(self.sum_energy_over_time))

        #return self.sum_energy, self.sum_over_time

