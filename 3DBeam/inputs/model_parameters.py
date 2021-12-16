from os.path import join as os_join
from os.path import sep as os_sep

''' 
a collection of parameters to be importet in the run module
Axes and Dimensons:

   CAARC has one long and one short side. The dimensons are defined such that for a angle of attack of 0° the long side is perpendicular to the wind.
   The wind direction is then y, 
   -> ly = 30m, lz = 45m
   weak bending is thus Mz with sway_z = 0.2 Hz
   and strong bending My with sway_y = 0.23 Hz
   (both for B uncoupled)
'''


# # AVAILABLE OPTIMIZATION OPTIONS AND INPUT PARAMETER SETTINGS
opt_methods = ['Nelder-Mead', 'SLSQP', 'Powell', 'BFGS', 'L-BFGS-B','TNC', 'CG']

optimization_variables = ['kya','kga','both']

weights = [[0.33,0.33,.33],[0.5,0.5,0.],[.4,.4,.2],[.45,.45,.1],[.8,.1,.1],[.1,.8,.1],[.2,.4,.4],[0.1,0.1,0.8]]

parameters = {'A':
               {  
                'dimension': '3D',
                'n_elements': 10,
                'lx_total_beam': 240,#180.0, #
                'material_density': 160,#150, 
                'E_Modul': 286100000.0,#8750.0,
                'nu': 0.1, #3D only
                'damping_coeff': 0.025,#0.001,#
                'nodes_per_elem': 2,
                'cross_section_area': 72*24,#1350,
                'B':72,
                'D':24, 
                'Iy': 746496.0,#105241.5,#1650.0,
                'Iz': 82944.0,#141085.0,#1750.0, 
                'I_param':100000.0, #NOTE if kga only 10000 else 100000 used for the coupling eintries to make it independet of the others -> initialy in the scale of Iz or Iy 
                'It': 829440.0,#270281.1,#3400.0,
                'modes_to_consider': 15,
                'static_load_magnitude': -20.0,
                'dynamic_load_file': os_join(*["inputs","forces","dynamic_force_11_nodes.npy"]),
                'inital_params_yg': [1.0,1.0,1.0],#[0.001,0.0012,0.0014]
                'params_k_ya': [0,0],#[1.0, 1.0, 1.0, 1.0] # omega: stiffness coupling y - a omega1: y -g
                'params_m_ya': [0.0,0.0,0.0],#[1.0, 1.0, 1.0] # omega: stiffness coupling, psi1, psi2: mass y-a, psi3: mass g-a
                'eigen_freqs_tar':[0.231,0.429,0.536],
                'eigen_freqs_orig':[0.231,0.429,0.536],
                
             },

            'B':{  
                'dimension': '3D',
                'n_elements': 10,
                'lx_total_beam': 180,
                'material_density': 100, #NOTE diminished from 160 to 100, this makes the inital frequency tunig in the correct direction BUT sense is questionable
                'E_Modul': 286100000.0,#2.861e8,
                'nu': 0.1, 
                'damping_coeff': 0.025,
                'nodes_per_elem': 2,
                'cross_section_area': 30*45,
                'B':45,
                'D':30,
                'Iy': 227812.5,
                'Iz': 101250.0,
                'I_param':10000.0, # 3 elems: 120000.0
                'It': 50000.0,#229062.0, #NOTE often set to Iy + Iz here diminsihed it by 100 000 such that the torsion optimization worked. 200 000 makes the coupling not worling 
                'modes_to_consider': 15,
                'static_load_magnitude': -20.0,
                'dynamic_load_file': os_join(*["inputs","forces","dynamic_force_11_nodes.npy"]),
                'inital_params_yg': [1.0,1.0,1.0],#[0.001,0.0012,0.0014]
                'params_k_ya': [0,0],#[1.0, 1.0, 1.0, 1.0] # omega: stiffness coupling y - a omega1: y -g
                'params_m_ya': [0.0,0.0,0.0],#[1.0, 1.0, 1.0] # omega: stiffness coupling, psi1, psi2: mass y-a, psi3: mass g-a
                'eigen_freqs_orig':[0.20,0.23,0.4], #[0.23,0.4,0.54]#
                'eigen_freqs_tar':[0.2591, 0.3249, 1.3555],#[1.62822, 2.04137, 8.517]  
                 #Andi
             }}
optimization_parameters = {'A':{
                            'eigen_freqs_tar': [0.231,0.429,0.536],# 'caarc_freqs_A'
                            'ratio_a_y_tar': 0.012,
                            'factor_y':0.9,
                            'coupling_target':'custom',
                            'readjust':{'sway_z':True,'sway_y':False,'torsion':False},
                            'consider_mode':0,
                            'var_to_optimize' :optimization_variables[2],#'ya'#'ga' both
                            'include_mass':False,
                            'method': opt_methods[1],
                            'init_guess':[0.,10.,0.0],#[-2.,5., 0.],#[0.1, 3] # k_ya, k_ga, m_cross
                            'bounds':((0.01, 100),(0.01, 100),(.0000001, 100)), #NOTE m seems to be corrected down always 
                            'weights':weights[0],
                            'save_optimized_parameters':False
                            },
                            'B':
                            {
                            'caarc_freqs_orig': [0.2,0.23,0.5], 
                            'eigen_freqs_tar':[0.2591, 0.3249, 1.3555],#[1.62822, 2.04137, 8.517],
                            'ratio_a_y_tar': 0.00372,#0.012,#
                            'factor_y':0.9,
                            'coupling_target':'semi_realistic', # custom, semi_realistic, realistic
                            'readjust':{'sway_z':True,'sway_y':False,'torsion':False},
                            'consider_mode':0,
                            'var_to_optimize' :optimization_variables[2],#'ya'#'ga' both
                            'include_mass':False,
                            'method': opt_methods[1],
                            'init_guess':[0.,10.,0],#3 elems: 0,10,0 # k_ya, k_ga, m_cross
                            'bounds':((0.01, 100),(0.01, 100),(.0000001, 100)), #NOTE m seems to be corrected down always 
                            'weights':weights[1],
                            'save_optimized_parameters':False
                            }}

dynamic_analysis_parameters = { "settings": 
                           {
                           "solver_type": "Linear",
                           "run_in_modal_coordinates": False,
                           "time":{
                                    "integration_scheme": "GenAlpha",
                                    "start": 0.0,
                                    "end": 595.0,
                                    "step" : 0.02},
                           "intial_conditions": {
                                    "displacement": None,
                                    "velocity": None,
                                    "acceleration" : None}
                           },
                        "input": 
                           {
                           "file_path": os_join(*["inputs","forces","dynamic_force_11_nodes.npy"])
                           }
                        }