import matplotlib.pyplot as plt
import copy
from os.path import join as os_join

from source.model import BeamModel
from source.optimizations import Optimizations
from source.postprocess import Postprocess
from source.utilities import utilities as utils
from plot_settings import plot_settings
from source.dynamic_analysis import DynamicAnalysis
from inputs import model_parameters

'''
coordinate system: x -> longitudinal axis, y -> perpendicular
Description:
A 3D Euler-Bernoulli Beam is extended to a torsion-bending coupled beam.
This is done by adjusting severeal Material parameters and introducing a coupling term into the stiffness matrix. 
3 Steps are done:
   1. Material parameters are optimized such that target Eigenfrequencies are reached. 
   2. The coupling paramters are optimized such the the eigen deformations are reached.
   3. The Material parameters are readjusted such that the Eigenfreuqencies are reached again

For the coupling optimization a target is required. Two different targets are implemented:
   - target based on a reference Eigenvalue analysis of a coupled model including the Eigenfrequencies and eigen deformations 
   - target based on the assumption that in a coupled case the bedning deformation is reduced by some factor whereas the torsional
     deformation increases in the same mode

'''
plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]

width = utils.cm2inch(5)
height = utils.cm2inch(8)

latex = False
savefig = False

plot_params = plot_settings.get_params(width =width, height=height, usetex=latex, minor_ticks=False)

#plt.rcParams.update({'figure.figsize': (width, height)})
plt.rcParams.update(plot_params)


# run_options = {
#                   'use_optimized_params':False,
#                   'plot_inital':True, # with ajdusted frequencies
#                   'plot_iteration':False,
#                   'plot_coupled':True,
#                   'plot_coupled_readjusted':False,
#                   'plot_obj_func':False, 
#                   'dynamic_analysis':[False,'b','reaction'],
#                   'static_analysis':False,
#                   'save_optimized_parameters':False,
#                  }

postprocess_parameters = {
                  'show_plots':True,
                  'savefig':False,
                  'optimization_procedure':{'intermediate_eigenmodes':True, 'objective_function':False},
                  'eigenmodes':{'do':False, 'n_modes':3, 'add_max_val':False, 'max_norm':False, 'rad_scale':True, 'caarc_A_only':False}, # with ajdusted frequencies
                  'static_analysis':{'do':False, 'load_directions':['y','a']},
                  'dynamic_analysis':{'do':True, 'time_hist':True, 'comp_stats':True,'comp_energy':True, 'normed':False,
                                                 'response':'total','response_type':'acceleration', 'node_id':10,'save_object':False, #response: use 'total' for acceleration total
                                                 'stats':['mean', 'std','max_est'],'unit':'MN', 'save_suffix':'_couple_uncouple_quer_1'},
                 }

model = 'B'

intermediate_results = {'eigenfrequencies':[],
                         'eigenmodes':[],
                         'dynamic_analysis':[],
                         'static_analysis':[]}

# # MODEL GEOMETRIC, MATERIAL AND ELEMENT PARAMETERS
parameters = model_parameters.parameters[model]
analysis_parameters = model_parameters.dynamic_analysis_parameters

# CREATE AN INITAL BEAM
# with initial frequency adjustments
beam = BeamModel(parameters, coupled=False, optimize_frequencies_init=True, use_translate_matrix=False)

# collect initial results
intermediate_results['eigenfrequencies'].append(beam.eigenfrequencies)
intermediate_results['eigenmodes'].append(beam.eigenmodes)

dynamic_analysis_init = DynamicAnalysis(beam, parameters = analysis_parameters)
dynamic_analysis_init.solve()
intermediate_results['dynamic_analysis'].append(dynamic_analysis_init)

beam.static_analysis_solve(apply_mean_dynamic=True, directions='all')
intermediate_results['static_analysis'].append(beam.static_deformation)

# # AVAILABLE OPTIMIZATION OPTIONS AND INPUT PARAMETER SETTINGS
                                                      
optimization_parameters = model_parameters.optimization_parameters[model]

# # TORSION COUPLING OPTIMIZATIONS
coupling_opt = Optimizations(beam, optimization_parameters)

var_to_optimize = optimization_parameters['var_to_optimize']
if var_to_optimize != 'both':
   # optimize onyl one stiffness variable
   coupling_opt.eigen_ya_stiffness_opt(which = var_to_optimize)   
else:
   # optimizing both stiffness variabels
   coupling_opt.eigen_vectorial_ya_opt(target_to_use=optimization_parameters['coupling_target'])

intermediate_results['coupling_optimization'] = coupling_opt

coupled_beam = beam.update_optimized_parameters(coupling_opt.optimized_design_params)
intermediate_results['eigenfrequencies'].append(coupled_beam.eigenfrequencies)
intermediate_results['eigenmodes'].append(coupled_beam.eigenmodes)

# # READJUSTMENT OF FREQUENCIES

opt_params_readjust = {'init_guess':[10,10], 'bounds':(0.001,100),'weights':None,'consider_mode':None,'method':None}
print ('\nOptimization of frequency after coupling...')

if optimization_parameters['readjust']['sway_z']: 
   freq_opt = Optimizations(model=coupled_beam, optimization_parameters=opt_params_readjust)
   print ('   ...readjusting stiffness for sway_z')
   freq_opt.adjust_sway_z_stiffness_for_target_eigenfreq(parameters['eigen_freqs_tar'][0], 
                                                         target_mode = 0,
                                                         print_to_console=True)
if optimization_parameters['readjust']['sway_y']: 
   freq_opt_y = Optimizations(model=coupled_beam, optimization_parameters=opt_params_readjust)
   print ('   ...readjusting stiffness for sway_y')
   freq_opt_y.adjust_sway_y_stiffness_for_target_eigenfreq(parameters['eigen_freqs_tar'][1], 
                                                         target_mode = 1,
                                                         print_to_console=True)

if optimization_parameters['readjust']['torsion']:
   freq_opt_torsion = Optimizations(model=coupled_beam, optimization_parameters=opt_params_readjust)                                                     
   print ('   ...readjusting stiffness for sway_torsion')
   freq_opt.adjust_torsional_stiffness_for_target_eigenfreq(parameters['eigen_freqs_tar'][2], 
                                                            target_mode = 2,
                                                            print_to_console=True)

intermediate_results['eigenfrequencies'].append(coupled_beam.eigenfrequencies)
intermediate_results['eigenmodes'].append(coupled_beam.eigenmodes)

dynamic_analysis_coupled = DynamicAnalysis(coupled_beam, parameters = analysis_parameters)
dynamic_analysis_coupled.solve()
intermediate_results['dynamic_analysis'].append(dynamic_analysis_coupled)

coupled_beam.static_analysis_solve(apply_mean_dynamic=True, directions='all')
intermediate_results['static_analysis'].append(coupled_beam.static_deformation)

# # POSTPROCESSING OBJECT
postprocess = Postprocess(postprocess_parameters, intermediate_results, beam, caarc_model=model)
            