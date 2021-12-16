import numpy as np
import matplotlib.pyplot as plt
import copy
import os
import pickle as pkl
from os.path import join as os_join
import json

from source.model import BeamModel
from source.optimizations import Optimizations
from source.postprocess import Postprocess
from source.utilities import utilities
from plot_settings import plot_settings
from source.dynamic_analysis import DynamicAnalysis
from inputs import model_parameters

plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]

width = utilities.cm2inch(6)
height = utilities.cm2inch(5)

mode = False #,1
pres = False

plot_params = plot_settings.get_params(width =width, height=height, usetex=mode, minor_ticks=False)

plt.rcParams.update({'axes.formatter.limits':(-3,3)}) 
plt.rcParams.update(plot_params)

# # INITIALIZING A POSTPROCESSING OBJECT
postprocess = Postprocess(show_plots = True, savefig = mode, savefig_latex = mode)

if postprocess.savefig or postprocess.savefig_latex:
   print ('\n****** FIGURE SAVINGS ACTIVE in CAARC B ******\n')

compare_options = {
                  'eigenmodes':{'do':False, 'n_modes':3, 'add_max_val':False, 'max_norm':False, 'rad_scale':True, 'caarc_A_only':False}, # with ajdusted frequencies
                  'static_analysis':{'do':False, 'load_directions':['y','a']},
                  'dynamic_analysis':{'do':True, 'time_hist':False, 'comp_stats':False,'comp_energy':True, 'normed':False,
                                                 'response':'total','response_type':'acceleration', 'node_id':10,'save_object':False, #response: use 'total' for acceleration total
                                                 'stats':['mean', 'std','max_est'],'unit':'MN', 'save_suffix':'_couple_uncouple_quer_1'},
                 }

parameters = model_parameters.parameters_B
n_nodes = parameters['n_elements'] + 1
parameters['dynamic_load_file'] = os_join(*['inputs','dynamic_force_'+str(int(n_nodes))+'_nodes.npy'])


# # CREATE THE UNCOUPLED BEAM
uncoupled_beam = BeamModel(parameters, coupled=False, optimize_frequencies_init=True, use_translate_matrix=False)

# # CREATE THE COUPLED BEAM
with open(os_join(*['optimized_parameters', 'coupled_beam_B.json']), 'r') as parameter_file:
    coupled_parameters = json.loads(parameter_file.read())
    coupled_beam = BeamModel(coupled_parameters, optimize_frequencies_init=False, use_translate_matrix=False)


if compare_options['eigenmodes']['do']:
    # figsize: w11, h8
    options = compare_options['eigenmodes']
    if not options['caarc_A_only']:
        postprocess.plot_eigenmodes_3D_compare(uncoupled_beam, coupled_beam, number_of_modes = 3, dofs_to_plot = ['y','a','z'], 
                                                add_max_deform= options['add_max_val'], max_normed=options['max_norm'],
                                                do_rad_scale =options['rad_scale'], filename_for_save = '0_no_name')
    else:
        postprocess.plot_eigenmodes_3D(uncoupled_beam, number_of_modes = 3, dofs_to_plot = ['y','a','z'], show_legend=False,
                                            add_max_deform= False, max_normed=False, caarc_A_only= options['caarc_A_only'],
                                            do_rad_scale = True, filename_for_save = 'eigenmodes_A')

if compare_options['static_analysis']['do']:
    # w11 h7
    options = compare_options['static_analysis']
    uncoupled_beam.static_analysis_solve(apply_mean_dynamic=True, directions=options['load_directions'])#all')
    coupled_beam.static_analysis_solve(apply_mean_dynamic=True, directions=options['load_directions'])#all')

    postprocess.plot_static_result(coupled_beam, init_deform = uncoupled_beam.static_deformation,
                                    load_type='mean', dofs_to_plot=options['load_directions'],#
                                    do_rad_scale=True, save_suffix='a_y_couple_uncouple', presentation=pres)

if compare_options['dynamic_analysis']['do']:
    options = compare_options['dynamic_analysis']
    analysis_parameters = model_parameters.dynamic_analysis_parameters

    if not os.path.isfile(os_join(*['output','dynamic_analysis_11_nodes_uncoupled.pkl'])):
        uncoupled_dynamic_analysis = DynamicAnalysis(uncoupled_beam, analysis_parameters)
        uncoupled_dynamic_analysis.solve()
    else:
        with open(os_join(*['output','dynamic_analysis_11_nodes_uncoupled.pkl']), 'rb') as dyn_input:
            uncoupled_dynamic_analysis = pkl.load(dyn_input)
        #options['save_object'] = False
    
    if not os.path.isfile(os_join(*['output','dynamic_analysis_11_nodes_coupled.pkl'])):
        coupled_dynamic_analysis = DynamicAnalysis(coupled_beam, analysis_parameters)
        coupled_dynamic_analysis.solve()
        
    else:
        with open(os_join(*['output','dynamic_analysis_11_nodes_coupled.pkl']), 'rb') as dyn_input:
            coupled_dynamic_analysis = pkl.load(dyn_input)
        options['save_object'] = False

    analyses = [uncoupled_dynamic_analysis, coupled_dynamic_analysis]
    dynamic_res_init = uncoupled_dynamic_analysis.solver

    if options['save_object']:
        
        filenames = ['dynamic_analysis_11_nodes_uncoupled.pkl','dynamic_analysis_11_nodes_coupled.pkl']
        for i,fname in enumerate(filenames):
            dest_path = ['output']
            dest_path.append(fname)
            with open(os_join(*dest_path),'wb') as dyn_output:
                pkl.dump(analyses[i], dyn_output)
            print ('\nsaved:', os_join(*dest_path))

    if options['time_hist']:
        # w7 h10 (time and freq übereinander) w8 h14.8 für quer format
        postprocess.plot_dynamic_results(coupled_dynamic_analysis, dof_label = options['response'], node_id = 10, 
                                            result_variable = options['response_type'], init_res = dynamic_res_init, 
                                            save_suffix = options['save_suffix'], include_fft=True, log=True, add_fft=False, unit=options['unit'])

    if options['comp_stats']:
        # w4.8 h3.5
        # analyses: list[uncoupled, coupled]
        postprocess.compare_stats(results = analyses, node_id=options['node_id'], response_label = options['response'], result_type = options['response_type'], 
                                  stats = options['stats'], unit=options['unit'], uncoupled_normed=options['normed'])

    if options['comp_energy']:
        uncoupled_dynamic_analysis.output_kinetic_energy(total=True)
        coupled_dynamic_analysis.output_kinetic_energy(total=True)
        postprocess.plot_compare_energies({'uncoupled': uncoupled_dynamic_analysis.sum_energy_over_time, 
                                           'coupled':coupled_dynamic_analysis.sum_energy_over_time})

