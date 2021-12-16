import numpy as np
import matplotlib.pyplot as plt
from os.path import join as os_join

import scipy.optimize as spo
from scipy.optimize import minimize, minimize_scalar
from scipy import linalg
from functools import partial
import source.postprocess

from source.utilities import utilities as utils 

class Optimizations(object):

    def __init__(self, model, optimization_parameters=None):
        
        ''' 
        an inital model is given to this object
        '''
        self.model = model
        self.opt_geometric_props = {}
        if optimization_parameters:
        
            self.opt_params = optimization_parameters
            self.consider_mode = optimization_parameters['consider_mode']
            self.method = optimization_parameters['method']
            self.weights = optimization_parameters['weights']
        else:
            self.opt_params = False


# FOR EIGENFREQUENCIES
    ''' 
    These are functions taken from ParOptBeam
    ''' 

    def adjust_sway_y_stiffness_for_target_eigenfreq(self, target_freq, target_mode, print_to_console=False):
        '''
        displacement in z direction -> sway_y = schwingung um y - Achse
        '''
        
        initial_iy = list(e.Iy for e in self.model.elements)

        # using partial to fix some parameters for the
        self.optimizable_function = partial(self.bending_y_geometric_stiffness_objective_function,
                                            target_freq,
                                            target_mode,
                                            initial_iy)
        init_guess = 1.0

        upper_bnd = self.model.elements[0].Iz / self.model.elements[0].Iy

        bnds_iy = (0.001, upper_bnd)#100)  # (1/8,8)

        # minimization_result = minimize(self.optimizable_function,
        #                                init_guess,
        #                                method='L-BFGS-B',  # 'SLSQP',#
        #                                bounds=(bnds_iy, bnds_a_sz))

        min_res = minimize_scalar(self.optimizable_function, tol=1e-06)#, options={'disp':True})

        # returning only one value!
        opt_fctr = min_res.x

        # NOTE this is only for constant Iy over the height
        self.opt_geometric_props['Iy'] = [min_res.x * iy_i for iy_i in initial_iy]

        if print_to_console:
            print('INITIAL Iy:', ', '.join([str(val) for val in initial_iy]))
            print()
            print('OPTIMIZED Iy: ', ', '.join([str(opt_fctr * val) for val in initial_iy]))
            
            print()
            print('FACTOR: ', opt_fctr)
            print()

    def bending_y_geometric_stiffness_objective_function(self, target_freq, target_mode, initial_iy, multiplier_fctr):

        for e in self.model.elements:
            e.Iy = multiplier_fctr * initial_iy[e.index]
            # assuming a linear dependency of shear areas
            # NOTE: do not forget to update further dependencies
            e.evaluate_relative_importance_of_shear()
            e.evaluate_torsional_inertia()

        # re-evaluate
        self.model.build_system_matricies(self.model.parameters['inital_params_yg'], 
                                          self.model.parameters['params_k_ya'], 
                                          self.model.parameters['params_m_ya'])

        self.model.eigenvalue_solve()

        eig_freq_cur = self.model.eigenfrequencies[target_mode]
        

        return (eig_freq_cur - target_freq)**2 / target_freq**2

    def adjust_sway_z_stiffness_for_target_eigenfreq(self, target_freq, target_mode, print_to_console=False):
        '''
        sway_z = schwingung in y richtung, um z Achse,
        '''
        initial_iz = list(e.Iz for e in self.model.elements)

        # using partial to fix some parameters for the
        self.optimizable_function = partial(self.bending_z_geometric_stiffness_objective_function,
                                            target_freq,
                                            target_mode,
                                            initial_iz)
        initi_guess = 1.0

        # NOTE this is correct only for homogenous cross section along length

        upper_bnd = self.model.elements[0].Iy / self.model.elements[0].Iz

        bnds_iz = (0.001, upper_bnd)#(0.001, 100)  # (1/8,8)

        # minimization_result = minimize(self.optimizable_function,
        #                                initi_guess,
        #                                method ='L-BFGS-B',
        #                                bounds = bnds_iz)

        min_res = minimize_scalar(self.optimizable_function, method='Bounded', tol=1e-06, bounds=bnds_iz)#, options={'disp':True})

        # returning only one value!
        #opt_iz_fctr = minimization_result.x
        opt_iz_fctr = min_res.x

        self.opt_geometric_props['Iz'] = [min_res.x * iz_i for iz_i in initial_iz]
        if print_to_console:
            print('  INITIAL iz:', ', '.join([str(val) for val in initial_iz]))
            print()
            print('  OPTIMIZED iz: ', ', '.join(
                [str(opt_iz_fctr * val) for val in initial_iz]))
            print()
            print('  FACTOR: ', opt_iz_fctr)
            print ('  Final Func:', min_res.fun)
            print()

    def bending_z_geometric_stiffness_objective_function(self, target_freq, target_mode, initial_iz, multiplier_fctr):
        
        for e in self.model.elements:
            e.Iz = multiplier_fctr * initial_iz[e.index]

            # NOTE: do not forget to update further dependencies
            e.evaluate_relative_importance_of_shear()
            e.evaluate_torsional_inertia()

        # re-evaluate
        self.model.build_system_matricies(self.model.parameters['inital_params_yg'], 
                                          self.model.parameters['params_k_ya'], 
                                          self.model.parameters['params_m_ya'])

        self.model.eigenvalue_solve()

        eig_freq_cur = self.model.eigenfrequencies[target_mode]        # mode_type_results is an ordered list

        result = (eig_freq_cur - target_freq)**2 / target_freq**2

        return result

    def adjust_torsional_stiffness_for_target_eigenfreq(self, target_freq, target_mode, print_to_console=False):
        initial_it = list(e.It for e in self.model.elements)
        initial_ip = list(e.Ip for e in self.model.elements)

        # NOTE: single parameter optimization seems not to be enough

        # using partial to fix some parameters for the
        self.optimizable_function = partial(self.torsional_geometric_stiffness_objective_function,
                                            target_freq,
                                            target_mode,
                                            initial_it,
                                            initial_ip)

        self.weights = [0.0,0.0,0.0]

        # NOTE: some additional reduction factor so that ip gets changes less

        init_guess = (1.0, 1.0)
        
        # NOTE: this seems not to be enough
        # bnds_it = (1/OptimizableStraightBeam.OPT_FCTR, OptimizableStraightBeam.OPT_FCTR)
        # bnds_ip = (1/OptimizableStraightBeam.OPT_FCTR, OptimizableStraightBeam.OPT_FCTR)

        # NOTE: seems that the stiffness contribution takes lower bound, the inertia one the upper bound
        bnds_it = (1/100, 10)
        bnds_ip = (1/11, 20)

        if self.opt_params:
            init_guess = self.opt_params['init_guess']
            bnds_it = self.opt_params['bounds']
            bnds_ip = self.opt_params['bounds']

        # NOTE: TNC, SLSQP, L-BFGS-B seems to work with bounds correctly, COBYLA not
        min_res = minimize(self.optimizable_function,
                                       init_guess,
                                       method='L-BFGS-B',
                                       bounds=(bnds_it, bnds_ip),
                                       options={'disp':False})

        # returning only one value!
        opt_fctr = min_res.x
        self.opt_geometric_props['It'] = [min_res.x[0] * it_i for it_i in initial_it]
        self.opt_geometric_props['Ip'] = [min_res.x[1] * ip_i for ip_i in initial_ip]
        if print_to_console:
            print('\nFACTORS It, Ip: ', ', '.join([str(val) for val in opt_fctr]))
            print ('final frequency: ', self.model.eigenfrequencies[target_mode])
            print()

    def torsional_geometric_stiffness_objective_function(self, target_freq, target_mode, initial_it, initial_ip, multiplier_fctr):

        for e in self.model.elements:
            e.It = multiplier_fctr[0] * initial_it[e.index]
            e.Ip = multiplier_fctr[1] * initial_ip[e.index]

        # re-evaluate
        self.model.build_system_matricies(self.model.parameters['inital_params_yg'], 
                                          self.model.parameters['params_k_ya'], 
                                          self.model.parameters['params_m_ya'])

        self.model.eigenvalue_solve()
        weights = [0]

        eig_freq_cur = self.model.eigenfrequencies[target_mode]

        return (eig_freq_cur - target_freq)**2 *100# / target_freq**2

# TORSION COUPLING OPTIMIZATIONS
    ''' 
    Coupling with one design variable.
    Either y-a or g-a
    ''' 
    
    def eigen_ya_stiffness_opt(self, which = 'kya'):
        ''' 
        Optimizes EITHER 'kya' or 'kga' to couple the y-displacement (gamma-rotation) to the torsional twist.
        The optimization target is hard coded in here -> see eigenmodes_target_*y*a
        The eigenfrequnecy_target is mostly not used (uncomment it in the objective function to see what happens).
        which: 'kya' or 'kga'
        ''' 
        if self.model.parameters['params_k_ya'] != [0.0,0.0]:
            raise Exception('inital parameters of ya are not 0 - check if sensible')

        eigenmodes_target_y = self.model.eigenmodes['y'][self.consider_mode]*0.9 # an assumption: y gets less if a is also deforming
        eigenmodes_target_a = np.linspace(0, eigenmodes_target_y[-1] * 0.012, eigenmodes_target_y.shape[0]) # 0.12 is the ratio of caarc tip a / tip y 1st mode
        eigenfreq_target = self.model.eigenfrequencies[self.consider_mode]

        self.inital = {'y':self.model.eigenmodes['y'][self.consider_mode],'a':self.model.eigenmodes['a'][self.consider_mode]}
        self.targets = {'y':eigenmodes_target_y, 'a':eigenmodes_target_a}

        self.optimizable_function = partial(self.obj_func_eigen_ya_stiffnes, self.consider_mode, 
                                            eigenmodes_target_y, eigenmodes_target_a, eigenfreq_target,
                                            which)
        ids = ['kya','kga']
        bounds = None
        #bounds = self.opt_params['bounds'][ids.index(which)]
        method_scalar = 'brent'
        #bounds = (0.001, 100)#,(0.001, 100))#,(0.001, 100))
        if bounds:
            method_scalar = 'bounded'

        res_scalar = minimize_scalar(self.optimizable_function, method=method_scalar, bounds= bounds, tol=1e-6)
        # SLSQP works with bounds
        #res_scalar = minimize(self.optimizable_function, x0= 0.0, method=self.method, bounds=bounds, tol=1e-6, options={'disp': True})

        # SLSQP works with constraints as well
        #res_scalar = minimize(self.optimizable_function, x0 = init_guess, method='SLSQP', constraints=cnstrts, tol=1e-3, options={'gtol': 1e-3, 'ftol': 1e-3, 'disp': True})

        #print( 'final F: ', str(self.optimizable_function))

        #self.optimized_design_params = res_scalar.x
        if which == 'kya':
            self.optimized_design_params = {'params_k_ya':[res_scalar.x, 0.0]}
        elif which == 'kga':
            self.optimized_design_params = {'params_k_ya':[0.0, res_scalar.x]}

        print('\noptimization result for design variable k'+which+':', res_scalar.x)

    def obj_func_eigen_ya_stiffnes(self, mode_id, eigenmodes_target_y, eigenmodes_target_a, eigenfreq_target, which, design_param):
        ''' 
        Objective function for one design variable (either kya or kga). 
        ''' 
        
        if isinstance(design_param, np.ndarray):
            if design_param.size == 2:
                if design_param[0] == design_param[1]:
                    design_param = design_param[0]
                else:
                    raise Exception('design parameter has 2 variables that differ')
            else:
                design_param = design_param[0]
        if which == 'kya':
            self.model.build_system_matricies(params_k_ya=[design_param, 0.0]) 
        elif which == 'kga':
            self.model.build_system_matricies(params_k_ya=[0.0, design_param])

        self.model.eigenvalue_solve()

        eigenmodes_cur = self.model.eigenmodes
        eigenfreq_cur = self.model.eigenfrequencies[self.consider_mode]

        f1 = utils.evaluate_residual(eigenmodes_cur['y'][mode_id], eigenmodes_target_y)
        f2 = utils.evaluate_residual(eigenmodes_cur['a'][mode_id], eigenmodes_target_a)
        #f3 = utils.evaluate_residual([eigenfreq_cur], [eigenfreq_target])

        weights = self.weights

        f = weights[0]*f1**2 + weights[1]*f2**2 # + weights[2] * f3**2

        return f

    ''' 
    Coupling with two design variables.
    Either y-a and g-a
    ''' 
    def eigen_vectorial_ya_opt(self, target_to_use = 'custom'):
        '''
        optimizing BOTH the stiffness coupling entries
            K_ya
            K_ga
        and the mass coupling entries
            mostly mass couling is not necessary or sensible - see also Thesis JZ ch. 3.3.3
            in optimization_parameters a boolean for turning this option on and off is used
            M_ya, M_yg (both with the same parameter )
        target_to_use:
            - 'custom': 0.9 times the initial lateral displacement & ratio alpha/disp = 0.012; a displacement is assumed linear
            - 'realistic': taking values from full 3D FE simulation of exccentirc building (ARiedls work)
            - 'semi_realistic': uses values from the optimization_params: 'ratio_a_y_tar', 'factor_y'; a twist displacement is amplified -> original shape is taken
        ''' 

        include_mass = self.opt_params['include_mass']

        # defining bounds
        # NOTE: k_ya takes lower bounds than 0.1
        bnds = self.opt_params['bounds']
        init_guess = self.opt_params['init_guess']#,1.0]#[0.0, 0.0,0.0]#[0.12, 0.15, 0.17] 
        self.n_iter = 0
        self.optimization_history = {'iter':[0],'func':[], 'k_ya':[init_guess[0]], 'k_ga':[init_guess[1]]}
        if include_mass:
            self.optimization_history['m_ya_ga'] = [init_guess[2]]
        def get_callback(x):
            # not really working
            self.n_iter += 1
            #self.optimization_history['func'].append(self.optimizable_function(x))
            self.optimization_history['k_ya'].append(x[0])
            self.optimization_history['k_ga'].append(x[1])
            self.optimization_history['iter'].append(self.n_iter)
            if include_mass:
                self.optimization_history['m_ya_ga'].append(x[2])
        def print_callback(x):
            print (x[0], x[1], x[2], self.optimizable_function(x))

        if self.model.parameters['params_k_ya'] != [0.0,0.0]:
            raise Exception('inital parameters of ya are not 0 - check if the targets are still sensible')

        if target_to_use == 'custom':
            eigenmodes_target_y = self.model.eigenmodes['y'][self.consider_mode]*0.9
            eigenmodes_target_a = np.linspace(0, eigenmodes_target_y[-1] * self.opt_params['ratio_a_y_tar'], eigenmodes_target_y.shape[0]) # 0.012 is the ratio of caarc tip a / tip y 1st mode
            eigenfreq_target = self.model.eigenfrequencies[self.consider_mode]

        elif target_to_use == 'realistic':
            modi = np.load(os_join(*['inputs', 'eigenvectors', 'EigenvectorsGid.npy']))
            z_coords = np.load(os_join(*['inputs','eigenvectors', 'z_coords_gid_45.npy']))
            # is only available with 45 nodes but is fitted if the current model has a different number of nodes
            if self.model.nodal_coordinates['x0'].size == 46:
                eigenmodes_target_y = modi[self.consider_mode][:,4]
                eigenmodes_target_a = modi[self.consider_mode][:,2] # here the ratio is 0.00373
            else:
                modi_fitted = utils.get_eigenform_polyfit(modi[self.consider_mode], z_coords, self.model.nodal_coordinates['x0'], plot_compare=False)
                eigenmodes_target_y = modi_fitted['eigenmodes']['y']
                eigenmodes_target_a = -1*modi_fitted['eigenmodes']['a']
            eigenfreq_target = self.opt_params['eigen_freqs_tar'] #self.model.eigenfrequencies[self.consider_mode]

        elif target_to_use == 'semi_realistic':
            ''' 
            assumes a reduction of y displacement by a custom factor < 1
            uses shape of a initial with a factor to get a-y tip ratio as specified 
            -> reason see inital shapes uncoupled max normed: a has the typical torsion shape -> needs amplification
            ''' 
            ratio_a_y = self.opt_params['ratio_a_y_tar']
            factor_y = self.opt_params['factor_y']
            eigenmodes_target_y = self.model.eigenmodes['y'][self.consider_mode]*factor_y

            a_factor = ratio_a_y * max(eigenmodes_target_y)/max(self.model.eigenmodes['a'][self.consider_mode])
            eigenmodes_target_a = self.model.eigenmodes['a'][self.consider_mode] * a_factor

            eigenfreq_target = self.opt_params['eigen_freqs_tar'] #self.model.eigenfrequencies[self.consider_mode]



        self.inital = {'y':self.model.eigenmodes['y'][self.consider_mode],'a':self.model.eigenmodes['a'][self.consider_mode]}
        self.targets = {'y':eigenmodes_target_y, 'a':eigenmodes_target_a}

        
        self.optimizable_function = partial(self.obj_func_eigen_vectorial_k_ya, self.consider_mode, eigenmodes_target_y, eigenmodes_target_a, eigenfreq_target, include_mass)
        self.optimization_history['func'].append(self.optimizable_function(init_guess))
        if not include_mass:
            print ('\nnot optimizing the mass entries, thus...')
            if len(bnds) != 2:
                bnds = bnds[:2]
                print ('  ...dropping the 3rd bound given')
            if len(init_guess) != 2:
                init_guess = init_guess[:2]
                print ('  ...dropping the 3rd initial guess given\n')

        # alternatively inequality constraints
        cnstrts = [{'type': 'ineq', 'fun': lambda x: 100 - x[0]},
                    {'type': 'ineq', 'fun': lambda x: 100 - x[1]},
                    {'type': 'ineq', 'fun': lambda x: 100 - x[2]},
                    {'type': 'ineq', 'fun': lambda x: x[0] - 0.001},
                    {'type': 'ineq', 'fun': lambda x: x[1] - 0.001},
                    {'type': 'ineq', 'fun': lambda x: x[2] - 0.001}]

        # SLSQP works with bounds
        res_scalar = minimize(self.optimizable_function,
                              x0 = init_guess,
                              method=self.method,
                              bounds=bnds, 
                              callback=get_callback,
                              options={'ftol': 1e-6, 'disp': True})

        evals = [0,10,10]
        #print ('func with manual opt params: ', self.optimizable_function(evals))

        self.optimized_design_params = {'params_k_ya':res_scalar.x[:2]}

        if include_mass:
            self.optimized_design_params['params_m_ya'] = [res_scalar.x[-1],res_scalar.x[-1],0.0]

        self.optimization_history['k_ya'].append(res_scalar.x[0])
        self.optimization_history['k_ga'].append(res_scalar.x[1])
        self.optimization_history['iter'].append(self.n_iter+1)

        digits = 5
        # SLSQP works with constraints as well
        # res_scalar = minimize(self.optimizable_function, x0 = init_guess, method='SLSQP', constraints=cnstrts, tol=1e-3, options={'gtol': 1e-3, 'ftol': 1e-3, 'disp': True})
        print()
        print('optimized parameters:')
        print ('  k_ya:', round(res_scalar.x[0],digits), 'absolute:', round(self.model.comp_k[1][3]))
        print ('  k_ga:', round(res_scalar.x[1],digits), 'absolute:', round(self.model.comp_k[3][5]))
        if include_mass:
            print ('  m_ya:', round(res_scalar.x[2],digits+4), 'absolute m_ya_11:', round(self.model.comp_m[1][3]), 'absolute m_ya_12:', round(self.model.comp_m[1][9]))

    def obj_func_eigen_vectorial_k_ya(self, mode_id, eigenmodes_target_y, eigenmodes_target_a, eigenfreq_target, include_mass, design_params):
        ''' 
        Objective function for more than one design variable (kya and kga optional mass entries).
        ''' 
        if include_mass:
            self.model.build_system_matricies(params_k_ya = design_params[:2], params_m_ya=[design_params[-1], design_params[-1],0])
        else:
            self.model.build_system_matricies(params_k_ya = design_params)

        self.model.eigenvalue_solve()

        eigenmodes_cur = self.model.eigenmodes
        eigenfreq_cur = self.model.eigenfrequencies[mode_id]

        f1 = utils.evaluate_residual(eigenmodes_cur['y'][mode_id], eigenmodes_target_y)
        f2 = utils.evaluate_residual(eigenmodes_cur['a'][mode_id], eigenmodes_target_a)
        f3 = utils.evaluate_residual([eigenfreq_cur], [eigenfreq_target])

        weights = self.weights

        gamma = 2
        components = [weights[0]*f1**gamma, weights[1]*f2**gamma, weights[2]*f3**gamma]
        f = sum(components)

        return f


# MASS MATRIX OPTIMIZATIONS
    
    def mass_entries_opt_ya(self):

        target = np.eye(self.model.n_dofs_node * self.model.n_elems)
        
        self.optimizable_function = partial(self.obj_func_gen_mass, target)

        bounds = self.opt_params['bounds']#,(0.001, 100))
        init_guess = self.opt_params['init_guess']#,1.0]#[0.0, 0.0,0.0]#[0.12, 0.15, 0.17] 


        #res_scalar = minimize_scalar(self.optimizable_function, method=method, bounds= bounds, options={'gtol': 1e-6, 'disp': True})
        # SLSQP works with bounds
        res_scalar = minimize(self.optimizable_function, x0= init_guess, method=self.method, bounds=bounds, options={'ftol': 1e-5, 'disp': True})

        print ('optimizaion result:', res_scalar.x)

    def obj_func_gen_mass(self, target, design_params):
        '''
        1. design_params are psi1, psi2 -> only ya entries the rest 0
        ''' 
        self.model.build_system_matricies(params_m_ya=[design_params[0],design_params[1], 0.0])

        eig_values_raw, eigen_modes_raw = linalg.eigh(self.model.comp_k, self.model.comp_m)
        
        gen_mass_cur = np.matmul(np.matmul(np.transpose(eigen_modes_raw), self.model.comp_m), eigen_modes_raw)
        
        f1 = utils.evaluate_residual(gen_mass_cur, target)

        return f1**2

