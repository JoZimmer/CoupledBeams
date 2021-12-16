import matplotlib.pyplot as plt 
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import cm
import copy
from os.path import join as os_join
from os.path import sep as os_sep
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

from source.utilities import utilities as utils
from source.utilities import global_definitions as GD
from source.utilities import CAARC_utilities as caarc_utils

'''
Legenden Außerhalb
unten 
axes[0,d_i].legend(bbox_to_anchor = (0.5, -0.2), loc ='upper center', ncol = 4)
rechts
axes[0,-1].legend(bbox_to_anchor = (1.04, 1), loc ='upper left')
'''

greek = {'y':'y','z':'z', 'x':'x','a':r'\alpha', 'b':r'\beta', 'g':r'\gamma'}

# custom rectangle size for figure layout
cust_rect = [0.05, 0.05, 0.95, 0.95]

COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

def cm2inch(value):
    return value/2.54

def convert_for_latex(string):
    l = list(string.split('_'))
    label = r'${}$'.format(l[0])
    for i in range(1,len(l)):
        if l[i] == 'alpha':
            l[i] = r'\alpha'
        label += r'$_{{{}}}$'.format(l[i])

    return label

def convert_load_for_latex(load_label,lower=False):
    ''' 
    use lower = True if line load and False if point load
    ''' 
    if lower:
        load = list(load_label)[0].lower()
    else:
        load = list(load_label)[0]
    direction = list(load_label)[1]

    return r'${}$'.format(load) + r'$_{{{}}}$'.format(GD.GREEK[direction])


class Postprocess(object):
    def __init__(self, parameters, intermediate_results, structure_model, caarc_model = 'B'):
        
        self.parameters = parameters
        self.show_plots = parameters['show_plots']
        self.savefig = parameters['savefig']
        self.nodal_coordinates = structure_model.nodal_coordinates
        self.rad_scale = np.sqrt(structure_model.parameters['cross_section_area'])

        #self.intermediate_results = intermediate_results

        dest_mode_results = os_join('plots',caarc_model,'eigenmode_results')
        dest_objective_func = os_join('plots',caarc_model,'objective_function')
        dest_1D_opt = os_join('plots',caarc_model,'ya_yg')
        dest_mass = os_join('plots',caarc_model,'mass_inclusion')
        dest_static = os_join('plots',caarc_model,'static_analysis')

        for result in intermediate_results:
            if result == 'eigenmodes':
                self.plot_eigenmodes_3D(intermediate_results['eigenfrequencies'], 
                                        intermediate_results['eigenmodes'],
                                        model=caarc_model,
                                        max_normed=False)
            if result == 'static_analysis':
                self.plot_static_result(static_deformation = intermediate_results[result][1], 
                                        init_deform = intermediate_results[result][0])
            if result == 'dynamic_analysis':
                options = self.parameters[result]
                if options['time_hist']:
                    # w7 h10 (time and freq übereinander) w8 h14.8 für quer format
                    self.plot_dynamic_results(intermediate_results[result][1], dof_label = options['response'], node_id = 10, 
                                                        result_variable = options['response_type'], init_res = intermediate_results[result][0].solver,
                                                        include_fft=True, log=True, add_fft=False, unit=options['unit'])

                if options['comp_stats']:
                    # w4.8 h3.5
                    # analyses: list[uncoupled, coupled]
                    self.compare_stats(results = intermediate_results[result], node_id=options['node_id'], response_label = options['response'], result_type = options['response_type'], 
                                            stats = options['stats'], unit=options['unit'], uncoupled_normed=options['normed'])

                if options['comp_energy']:
                    intermediate_results[result][0].output_kinetic_energy(total=True)
                    intermediate_results[result][1].output_kinetic_energy(total=True)
                    self.plot_compare_energies({'uncoupled': intermediate_results[result][0].sum_energy_over_time, 
                                                    'coupled':intermediate_results[result][1].sum_energy_over_time})

            if result == 'coupling_optimization':
                if parameters['optimization_procedure']['objective_function']:
                    self.plot_objective_function_3D(optimization_object= coupling_opt, evaluation_space_x=[-55, 55, 0.5], evaluation_space_y=[-30, 80, 0.5],
                                                include_opt_history=True, filename_for_save='objective_func_10_elems' + var_to_optimize, save_evaluation= True)
# # EIGENMODES

    # # 3D

    def plot_eigenmodes_3D(self,eigenfrequencies, eigenmodes, model = 'B',number_of_modes = 3, dofs_to_plot = ['y','z','a'],
                            max_normed = False, opt_targets = None, initial = None):

        ''' 
        Plotting eigendmodes of a 3D model.
            beam_model: a beam model object optimized or not 
            model: 'A' or 'B' for setting some plot style options
            opt_targets: the eigenmode target deformations dictionary with dofs or None
            initial: initial eigenmode deformations dictionary with dofs or None
            number_of_modes: ...to show 
            dofs_to_plot: ...
            max_normed: boolean -> if deformation should be normed with theri respective maximum
            rad_scale: scale for the rotational deformations (with dimensons from the model)
        ''' 
        steps = len(eigenfrequencies)

        if number_of_modes == 1:
            fig, ax = plt.subplots(nrows = steps, ncols = number_of_modes,  num='eigenmode results')#figsize=(2,3.5),
        
            ax = [ax]
        else:
            fig, ax = plt.subplots(nrows = steps,ncols = number_of_modes, sharey=True,  num='eigenmode results')#figsize=(5,4),

        norm, norm2, norm3 = 1, 1, 1
        if max_normed:
            self.rad_scale = False
        for step in range(steps):
            for i in range(number_of_modes):
                ax[step][i].plot( np.zeros(len(self.nodal_coordinates['x0'])),
                            self.nodal_coordinates['x0'],
                            #label = r'$structure$',
                            color = 'grey',
                            marker = 'o', 
                            linestyle = '--')
                ax[step][i].set_title(r'$mode$ ' + r'${}$'.format(str(i+1)) + '\n' +r'$f=$ ' + r'${}$'.format(str(round(eigenfrequencies[step][i],3))) +r' $Hz$')
                
                for d_i, dof in enumerate(dofs_to_plot):
                    scale=1.0
                    if self.rad_scale:
                        if dof in ['a','b','g']:
                            scale = self.rad_scale
                        
                    y = utils.check_and_flip_sign_array(eigenmodes[step][dof][i])

                    if i == 0:
                        if initial:
                            if dof in initial.keys():
                                y3 = initial[dof]
                                #print ('\ntarget max', dof, max(y3))
                                if max_normed:
                                    norm3 = 1/max(y3)
                                ax[step][i].plot(y3*norm3*scale,
                                            self.nodal_coordinates['x0'],
                                            label = r'${}$'.format(greek[dof]) + r'$_{initial}$',
                                            linestyle = ':',
                                            color = COLORS[d_i])
                        if opt_targets:
                            if dof in opt_targets.keys():
                                y2 = opt_targets[dof]
                                label_tar = r'${}$'.format(greek[dof]) + r'$_{target}$'
                                if max_normed:
                                    norm2 = 1/max(y2)
                                ax[step][i].plot(y2*norm2*scale,
                                            self.nodal_coordinates['x0'],
                                            label = label_tar,
                                            linestyle = '--',
                                            color = COLORS[d_i])

                    # ACTUAL MODE SHAPES
                    if max_normed:
                        norm = 1/max(y)
                    label = r'${}$'.format(greek[dof])
                    
                    ax[step][i].plot(y*norm*scale,
                                self.nodal_coordinates['x0'],
                                label = label,
                                linestyle = '-',
                                color = COLORS[d_i])
                
                ax[step][i].legend(loc = 'lower right')
                ax[step][i].grid()
                bins = {'A':6,'B':4}
                ax[step][i].locator_params(axis='y', nbins = bins[model])
                top_y = {'A':250,'B':200}
                ax[step][i].set_ylim(bottom=0,top=top_y[model])
                ax[step][i].set_xlabel(r'$defl. \, [m] $')
                ax[step][i].xaxis.set_label_coords(0.15, -0.1)
                ax[step][0].set_ylabel(r'height $[m]$') 

            ratio = max(utils.check_and_flip_sign_array(eigenmodes[step]['a'][0])) / max(utils.check_and_flip_sign_array(eigenmodes[step]['y'][0]))
            ratio_label = r'$\alpha_{max}/y_{max} = $' + str(round(ratio,3))

            props = dict(boxstyle='round', facecolor='white', edgecolor = 'lightgrey', alpha=0.8)

            # place a text box in upper left in axes coords
            ax[step][0].text(0.1, 0.97, ratio_label, transform=ax[step][0].transAxes, verticalalignment='top', bbox=props)

        opt_steps = ['step{}'.format(s_i) for s_i in range(1,steps+1)]
        pad = 5
        for ax_i, row in zip(ax[:,0], opt_steps):
            ax_i.annotate(row, xy=(0, 0.5), xytext=(-ax_i.yaxis.labelpad - pad, 0),
                        xycoords=ax_i.yaxis.label, textcoords='offset points',
                        size='large', ha='right', va='center')
                
        if self.savefig:
            plt.savefig(dest_mode_results + filename_for_save)
            print ('\nsaved: ', dest_mode_results + filename_for_save)
        if self.show_plots:
            plt.show()

# # OPTIMIZATIONS

    def plot_objective_function_2D(self,optimization_object, evaluation_space = [0,10, 0.01],design_var_label = 'design variable'):

        print ('\n EVALUATE AND PLOT OBJECTIVE FUNCTION\n')

        fig, ax = plt.subplots(figsize=(5,3), num='objective_func_1D')

        objective_function = optimization_object.optimizable_function
        x = np.arange(evaluation_space[0], evaluation_space[1], evaluation_space[2])
        result = np.zeros(len(x))
        for i, val in enumerate(x):
            result[i] = objective_function(val)
        
        ax.plot(x, result)
        opt_res = optimization_object.optimized_design_params
        ax.plot(opt_res, objective_function(opt_res),linestyle = 'None', marker='o',mfc='r',mec='k', ms=4, label='optimized variable ' + str(round(opt_res,4)))

        if optimization_object.optimized_design_params:
            ax.vlines(optimization_object.optimized_design_params, min(result), max(result), 
                        #label = 'optimized variable: ',# + str(round(optimization_object.final_design_variable, 2)), 
                        color = 'r', 
                        linestyle ='--')
        ax.set_title('optimizable function')
        ax.set_xlabel('values of ' + design_var_label )
        ax.set_ylabel(r'$ f = \sum w_{i} * e_{i} ^{2}$')
        ax.grid()
        ax.legend()
        # if self.savefig:
        #     plt.savefig(dest_mode_results)

        # if self.show_plots:
        #     plt.show()
        
        plt.close()

    def plot_objective_function_3D(self,optimization_object, evaluation_space_x = [-4,4,0.1], evaluation_space_y = None,
                                    include_opt_history = False, fig_label = '', filename_for_save ='0_no_name', add_3D_surf_plot = False, save_evaluation = True):
        '''
        deepcopy of the optimization is created to avoid undesired changes in the base objects of this class.
            evaluation_space[0,1]: space where the optimizable funciton shall be evaluated
            evaluation_space[2]: steps of evaluation in this space
            if space_y is not given the same space for x and y is used
        '''
        print ('\nEVALUATE AND PLOT OBJECTIVE FUNCTION...\n')
        print (filename_for_save + str(optimization_object.weights))
        optimization_obj_eval = copy.deepcopy(optimization_object)
        fig_title = 'objective_func_' + fig_label
        fig = plt.figure(num=fig_title)#, figsize=(5,3))
        if add_3D_surf_plot:
            ax1 = fig.add_subplot(121)
            ax2 = fig.add_subplot(122, projection='3d')
        else:
            ax1 = fig.add_subplot(111)
        if not self.savefig:
            ax1.set_title(r'objective function with weights: ' + r'${}$'.format(optimization_obj_eval.weights))

        objective_function = optimization_obj_eval.optimizable_function

        # get the optimization information before the function evaluation
        # needs explictly to be copied since the object changes during function evaluation
        if include_opt_history:
            k_ya_hist = np.copy(optimization_obj_eval.optimization_history['k_ya'])
            k_ga_hist = np.copy(optimization_obj_eval.optimization_history['k_ga'])
            func_hist = np.copy(optimization_obj_eval.optimization_history['func'])

        if not evaluation_space_y:
            evaluation_space_y = evaluation_space_x

        x = np.arange(evaluation_space_x[0],evaluation_space_x[1],evaluation_space_x[2])
        y = np.arange(evaluation_space_y[0],evaluation_space_y[1],evaluation_space_y[2])

        if x.shape != y.shape:
            raise Exception('shape of x and y input parameters must be the same for 3D plotting')
        x,y = np.meshgrid(x,y)
        z = np.zeros((len(x),len(y)))
        for i in range(z.shape[0]):
            for j in range(z.shape[1]):

                z[i][j] = objective_function((x[i][j], y[i][j]))

        if save_evaluation:
            fname = os_join(*['objective_functions',utils.join_whitespaced_string(filename_for_save)+'x.npy'])
            np.save(fname, x)
            fname = os_join(*['objective_functions',utils.join_whitespaced_string(filename_for_save)+'y.npy'])
            np.save(fname, y)
            fname = os_join(*['objective_functions',utils.join_whitespaced_string(filename_for_save)+'z.npy'])
            np.save(fname, z)

        levels = np.linspace(0.1, 0.75, 2000)

        cs = ax1.contourf(x,y,z, 50, cmap = 'viridis')
        l = list(cs.levels)
        level_lines = l[:2] + l[2:10:2]
        cs2 = ax1.contour(cs, levels = level_lines, colors = 'r', linewidths= 0.5)

        if include_opt_history:
            ax1.plot(k_ya_hist, 
                    k_ga_hist, 
                    linestyle = '-', color = 'grey', #linewidth= # line stuff
                    marker = 'o', mfc = 'grey', mec='k', ms = 3,  # marker stuff
                    label='iterations')
            ax1.scatter(k_ya_hist[0],k_ga_hist[0], marker='o',c='lime',edgecolors='k', s = 20,label='start',zorder=5)
            ax1.scatter(k_ya_hist[-1],k_ga_hist[-1], marker='o',c='red',edgecolors='k', s= 20, label='end',zorder=5)

        
        
        if add_3D_surf_plot:
            surf = ax2.plot_surface(x,y,z,
                                cmap= 'hsv',
                                rstride=1, cstride=1, 
                                linewidth=0,
                                antialiased=False,
                                vmin = z.min(), vmax = z.max())

            ax2.plot_wireframe(x,y,z)
            cbar = fig.colorbar(surf, shrink=0.5, aspect=20, ax = ax2)#,extend={'min','max'})

            
            ax2.set_xlabel('k_ya')
            ax2.set_ylabel('k_ga')
            ax2.set_zlabel(r'$ f = \sum^{3} w_{i} * e_{i} ^{2}$')
            ax2.grid()

        ax1.set_xlabel(r'$k_{ya}$')
        ax1.set_ylabel(r'$k_{ga}$')

        #ax1.set_xlim(-25,25)
        #ax1.margins(x = -0.3)
        #ax1.set_ylim(-10,65)
        #self.set_ax_size(utils.cm2inch(6), utils.cm2inch(4), ax=ax1)

        #cs_bar.ax.set_xlabel(r'$ f = (\mathbf{x})$')
        cs_bar = fig.colorbar(cs,  shrink=0.5, aspect=20, ax = ax1, pad=0.001)
        cs_bar.add_lines(cs2)
        cs_bar.set_label(r'$ f = (\mathbf{x})$')
        ax1.grid( linestyle='--')
        ax1.legend()

        if self.savefig:
            plt.savefig(dest_objective_func + filename_for_save)
            print ('\nsaved:',dest_objective_func + filename_for_save)

        if self.show_plots:
            plt.show()
        plt.close()
        
        del optimization_obj_eval

    def plot_optimization_history(self, optimization_object, include_func, norm_with_start = True):

        fig = plt.figure(num='opt_history')#, figsize=(cm2inch(7.3), cm2inch(4.8))

        for key, val in optimization_object.optimization_history.items():
            if not include_func and key == 'func':
                continue
            if key == 'iter':
                continue
            if norm_with_start:
            #     if val[0] == 0:
            #         val[0] = 0.01
                val_norm = [val_i - val[0] for val_i in val]
            else:
                val_norm = val
            plt.plot(np.arange(1,len(val)+1), val_norm, label=utils.prepare_string_for_latex(key))
            plt.xlabel('Iteration')
            plt.xlim(left=1)
            plt.ylabel(r'$x - x_{start}$')
            plt.grid()
        #plt.legend()
        plt.legend(bbox_to_anchor = (1.04, 1), loc ='upper left')

        if self.savefig:
            plt.savefig(dest_mass + 'mass_inc_iter')

        if self.show_plots:
            plt.show()

    def plot_multiple_result_vectors(self,beam_model, vector_list):
        '''
        the computed result of the respective targets(egenform or static disp.) are tracked during optimization
        here some of the are plotted to see the developments
        '''
        fig, ax = plt.subplots()

        if beam_model.opt_params['optimization_target'] == 'static_tip_disp':
            w_analytic = utils.analytic_function_static_disp(beam_model.parameters, np.arange(0,beam_model.parameters['lx_total_beam']+1))
            ax.plot(np.arange(0,len(w_analytic)),
                                w_analytic,
                                label = 'analytic',
                                color = 'r',
                                linestyle = '--',
                                marker = 'o')

        elif beam_model.opt_params['optimization_target'] == 'frequency':
            y = utils.analytic_eigenmode_shapes(beam_model)[0]
            ax.plot(beam_model.nodal_coordinates['x0'],
                                y,
                                label = 'correct 1st mode',
                                color = 'r',
                                linestyle = '--',
                                marker = 'o')
            
        ax.plot(beam_model.nodal_coordinates['x0'],
                beam_model.nodal_coordinates['y0'],
                label = 'structure',
                color = 'black',
                linestyle = '--')

        if len(vector_list) > 50:
            vectors_to_plot = vector_list[::int(len(vector_list)/10)]
        else:
            vectors_to_plot = vector_list
        for i, vector in enumerate(vectors_to_plot):
            if i == 0:
                ax.plot(beam_model.nodal_coordinates['x0'],
                        vector,
                        color = 'green',
                        linestyle= '--',
                        marker = 'o',
                        label = 'initial guess')
            elif i < len(vectors_to_plot)-2:
                ax.plot(beam_model.nodal_coordinates['x0'],
                        vector,
                        color = 'grey',
                        linestyle= '--')
            elif i < len(vectors_to_plot)-1:
                ax.plot(beam_model.nodal_coordinates['x0'],
                        vector,
                        color = 'grey',
                        linestyle= '--',
                        label = 'intermediate steps')
            else:
                ax.plot(beam_model.nodal_coordinates['x0'],
                        vector,
                        color = 'tab:blue',
                        linestyle= '-',
                        label = 'last iteration')

        ax.legend()
        ax.grid()
        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_title('target during optimization ' + beam_model.opt_params['method'] + ' init: '+ str(round(beam_model.init_guess,1)) + 
                        ' scaling exponent: ' + str(beam_model.opt_params['scaling_exponent']))
        plt.show()

# # STATIC ANALYSIS

    def plot_static_result(self,static_deformation, init_deform = None):

        title = 'static results '
        dofs_to_plot = self.parameters['static_analysis']['load_directions']

        fig, ax = plt.subplots(num=title )

        ax.plot(self.nodal_coordinates['y0'],
                self.nodal_coordinates['x0'],
                label = 'structure',
                marker = 'o',
                color = 'grey',
                linestyle = '--')

        ratio_text = ''
        next_line = '\n'
        dof_colors = {'y':'tab:blue','a':'tab:orange','z':'tab:green'}
        for d_i, dof in enumerate(dofs_to_plot):
            scale=1.0
            if self.rad_scale:
                if dof in ['a','b','g']:
                    scale = self.rad_scale
            label1 = 'coupled:     ' + r'${}$'.format(greek[dof]) + r'$_{max} =$'+ '{0:.2e}'.format(static_deformation[dof][-1][0]*scale)

            ax.plot(static_deformation[dof] * scale,
                    self.nodal_coordinates['x0'],
                    label = label1,
                    color = dof_colors[dof])
            if init_deform:
                label2 = 'uncoupled: ' + r'${}$'.format(greek[dof]) + r'$_{max} =$' + '{0:.2e}'.format(init_deform[dof][-1][0]*scale)
                ax.plot(init_deform[dof] * scale,
                        self.nodal_coordinates['x0'],
                        label = label2,
                        color = dof_colors[dof],
                        linestyle = '--')
            
            if d_i == len(dofs_to_plot)-1:
                next_line = ''

            ratio_un_c =  static_deformation[dof][-1][0]/init_deform[dof][-1][0]
            digits = 2#3 - len(str(int(ratio_un_c)))
            if abs(ratio_un_c) > 100:
                ratio_un_c = int(ratio_un_c)

            ratio_text += 'coupled '+ r'${}$'.format(greek[dof]) + r'$_{max}$' + r'$=$'+ r'${}$'.format(round(ratio_un_c,digits)) + ' of uncoupled' + next_line
            print ('  coupled deformation', greek[dof], 'is', round(ratio_un_c,2), 'of uncoupled.')
        
        txt_font = None

        print ('  w_max ist beam:     ', static_deformation['y'][-1][0])
        ax.legend(bbox_to_anchor = (1.04, 1), loc ='upper left')#, fontsize=txt_font)
        props = dict(boxstyle='round', facecolor='white', edgecolor = 'lightgrey', alpha=0.8)

        # place a text box in upper left in axes coords
        if txt_font == None:
            h = 0.4
        else:
            h = 0.6
        ax.text(1.09, h, ratio_text, transform=ax.transAxes, verticalalignment='top', bbox=props)#, fontsize=txt_font)

        ax.grid()
        ax.set_xlabel(r'defl. $[m]$')
        ax.set_ylabel(r'height $[m]$')
        ax.locator_params(axis='y', nbins = 4)
        ax.set_ylim(bottom=0,top=200)
        variant = ''

        if self.savefig:
            plt.savefig(dest_static + 'static_results')
            print ('\nsaved:', dest_static + 'static_results')
        if self.show_plots:
            plt.show()

# # DYNAMIC ANALYSIS

    def plot_dynamic_results(self, dynamic_analysis, dof_label, node_id, result_variable, init_res = None, include_fft = True,
                            log = False, add_fft = False, unit='N'):
        ''' 
        dynamic_analyis: analysis object
        dof_label: label of dof to plot
        node_id: node at whcih results to plot (starting from 0)
        result_variable: 'displacement','acceleration','reaction'
        ''' 
        dest = os_join(*['plots','CAARC_B','dynamic_results'])
        if result_variable != 'acceleration':
            dof = GD.DOF_LABELS['3D'].index(dof_label) + (node_id * GD.n_dofs_node['3D'])
        dof_fixed = dof_label
        unit_scale = 1
        unit_label = ''
        if result_variable == 'displacement':
            result_data = dynamic_analysis.solver.displacement[dof, :]
            if init_res:
                init_data = init_res.displacement[dof, :]
        elif result_variable == 'velocity':
            result_data = dynamic_analysis.solver.velocity[dof, :]
            if init_res:
                init_data = init_res.velocity[dof, :]
        elif result_variable == 'acceleration':
            # if dof_fixed == 'a':
            #     unit_label = r'[rad/s^{2}]'
            # else:
            unit_label = r'[m/s^{2}]'
            #dof_label = 'a_total'
            directions = {'total':{'dofs':['y','z','a'],'label':'a_total','rad':np.sqrt(15**2+22.5**2)},
                          'y':{'dofs':['y','a'],'label':'a_y','rad':22.5},'z':{'dofs':['z','a'],'label':'a_z','rad':15},
                          'a':{'dofs':['a'],'label':'a_alpha','rad':27}}

            dof_ids = [GD.DOF_LABELS['3D'].index(dof_i) + (node_id * GD.n_dofs_node['3D']) for dof_i in directions[dof_label]['dofs']]
            result_data, init_data = np.zeros(dynamic_analysis.array_time.size), np.zeros(dynamic_analysis.array_time.size)
            for i, dof_id in enumerate(dof_ids):
                if GD.DOF_LABELS['3D'][dof_id - node_id * GD.n_dofs_node['3D']] == 'a':
                    rad_scale =  directions[dof_label]['rad']
                else:
                    rad_scale = 1.0
                result_data += (dynamic_analysis.solver.acceleration[dof_id, :]*rad_scale)#**2
                init_data += (init_res.acceleration[dof_id, :]*rad_scale)#**2
            # result_data = np.sqrt(result_data)
            # init_data = np.sqrt(init_data)
            dof_label = directions[dof_label]['label']

        elif result_variable == 'action':
            result_data = dynamic_analysis.force[dof, :]
        elif result_variable == 'reaction':
            unit_label = GD.UNITS_POINT_LOAD_DIRECTION[dof_label]
            unit_label = unit_label.replace(unit_label[1], unit)
            unit_scale = GD.UNIT_SCALE[unit]
            if dof in dynamic_analysis.structure_model.dofs_of_bc:# or dof in dynamic_analysis.structure_model.elastic_bc_dofs:
                result_data = dynamic_analysis.solver.dynamic_reaction[dof, :] * unit_scale
            else:
                print ('\nReplacing the selected node by the ground node for reaction results')
                node_id = 0
                dof = GD.DOF_LABELS['3D'].index(dof_label)
                result_data = dynamic_analysis.solver.dynamic_reaction[dof, :]  * unit_scale
            if init_res:
                init_data = init_res.dynamic_reaction[dof, :] * unit_scale
            dof_label = GD.DIRECTION_RESPONSE_MAP[dof_label] 
            

        digits = 6
        mean = round(np.mean(result_data), digits) 
        std = round(np.std(result_data), digits) 
        if init_res:
            mean_init = round(np.mean(init_data), digits)
            std_init = round(np.std(init_data), digits) 


        plot_title = result_variable.capitalize() + ' at node ' + str(node_id) + ' in ' + dof_label + ' direction'

        if include_fft:
            fig, ax = plt.subplots(nrows=2, num = 'dyn_res')
        else:
            fig, ax = plt.subplots(num = 'dyn_res')
            ax = [ax]


        ax[0].set_xlabel(r'$t [s]$',fontsize= 10)
        ax[0].set_ylabel(convert_for_latex(dof_label) + r' ${}$'.format(unit_label), fontsize= 10)

        # UNCOUPLED
        res3 = ax[0].hlines(mean_init, dynamic_analysis.array_time[0], dynamic_analysis.array_time[-1], label='uncoupled mean', color = 'k', zorder=3)
        res4 = ax[0].hlines(mean_init + std_init, dynamic_analysis.array_time[0], dynamic_analysis.array_time[-1], label='uncoupled mean +/- std', color = 'k', 
                                                                                               linestyle= '--', zorder=3)
        ax[0].hlines(mean_init - std_init, dynamic_analysis.array_time[0], dynamic_analysis.array_time[-1], color = 'k', linestyle= '--', zorder=3)

        # COUPLED
        res5 = ax[0].hlines(mean, dynamic_analysis.array_time[0], dynamic_analysis.array_time[-1], label='coupled mean', color = 'r', zorder=3)
        res6 = ax[0].hlines(mean + std, dynamic_analysis.array_time[0], dynamic_analysis.array_time[-1], label='coupled mean +/- std', color = 'r', 
                                                                                               linestyle= '--', zorder=3)
        ax[0].hlines(mean - std, dynamic_analysis.array_time[0], dynamic_analysis.array_time[-1], color = 'r', linestyle= '--', zorder=3)
        
        label2 = None
        if init_res:
            res1, = ax[0].plot(dynamic_analysis.array_time, init_data, linewidth = 0.3, linestyle = '-', label = 'uncoupled result', color = 'gray')
            label2 = 'coupled result'

        res2, = ax[0].plot(dynamic_analysis.array_time, result_data, linewidth = 0.3,  label = label2)
        
        #ax[0].legend()
        if include_fft:
            leg_ax = 1
        else:
            leg_ax = 0
        
        #ax[leg_ax].legend(handles=[res1,res2,res3,res4,res5,res6],bbox_to_anchor = (0.5, -0.25), loc ='upper center', ncol = 2)#, fontsize=10)
        ax[0].set_xlim(0, dynamic_analysis.array_time[-1])
        if result_variable == 'acceleration':
            if dof_fixed == 'a':
                ax[0].set_ylim(bottom = -0.1,top= 0.1)
            else:
                ax[0].set_ylim(top= 0.8)

        f_orig = [0.2,0.23,0.5]
        if include_fft:
            eig_freqs = dynamic_analysis.structure_model.eigenfrequencies
            if result_variable != 'acceleration':
                sway_naming = GD.MODE_CATEGORIZATION_REVERSE['3D'][dof_fixed]
                f_id = GD.CAARC_MODE_DIRECTIONS['0_deg'][sway_naming]
                natural_f = eig_freqs[f_id]
                f_unc = f_orig[f_id]

            f_3 = eig_freqs[:3]
            sampling_freq = 1/dynamic_analysis.dt

            freq_half, series_fft = utils.get_fft(result_data, sampling_freq)
            freq_half_un, series_fft_un = utils.get_fft(init_data, sampling_freq)

            ax[1].plot(freq_half_un, series_fft_un, linewidth = 0.3, color='gray')
            ax[1].plot(freq_half, series_fft, linewidth = 0.3)
            
            if log:
                ax[1].set_xscale('log')
                ax[1].set_yscale('log')          
        
            if not log:
                ax[1].set_xlim(0.0,1.5)
                ax[1].set_ylim(bottom = 0, top = 0.06)

            ax[1].set_ylabel(r'$S($' + convert_for_latex(dof_label) + r'$)$',fontsize= 10)
            ax[1].set_xlabel(r'$f [Hz]$',fontsize= 10)
            text_pos_co={'Mz':(0.23,0.025),'My':(0.32,0.0035),'Mx':(0.85, 0.0000025),
                        'a_y':[0.25, 30e-8,1.36, 30e-8],'a_z':[0.34, 30e-8,1.36, 30e-8],'a_alpha':[0.13, 30e-8,0.7, 5e-8],'a_total':[0.25, 30e-8,1.36, 30e-8]}
            text_pos_un={'Mz':(0.12,0.05),'My':(0.14,0.05),'Mx':(0.3, 0.0005)}
            if result_variable == 'reaction':
                ax[1].axvline(natural_f, linestyle = '--', color = 'tab:blue', linewidth= 0.35,)
                ax[1].axvline(f_unc, linestyle = '--', color = 'gray', linewidth= 0.35,)
                # ax[1].axvline(eig_freqs[3], linestyle = '--', color = 'r', linewidth= 0.35,)
                # ax[1].axvline(eig_freqs[6], linestyle = '--', color = 'r', linewidth= 0.35,)
                        #label = r'$f$ ' + other_utils.prepare_string_for_latex(sway_naming) + r' = ${}$'.format(round(natural_f,2)) + r' $[Hz]$')
                ax[1].text(x = text_pos_co[dof_label][0],y = text_pos_co[dof_label][1], s= r'$f_{cou.}$ ' ,rotation = 90, color='tab:blue', fontsize=10)#
                ax[1].text(x = text_pos_un[dof_label][0],y = text_pos_un[dof_label][1], s =r'$f_{unc.}$ ',rotation = 90, color='gray',fontsize=10)# 
            
            if result_variable == 'acceleration':
                if dof_label == 'a_z':
                    fs = [f_3[1],f_3[2]]
                    lbl1 =  r'$f_{cou.,2}$ '
                else:
                    fs = [f_3[0],f_3[2]]
                    lbl1 = r'$f_{cou.,1}$ '
                
                for f_i in fs:
                    ax[1].axvline(f_i, linestyle = '--', color = 'tab:blue', linewidth= 0.35,)
                ax[1].text(text_pos_co[dof_label][0], text_pos_co[dof_label][1], lbl1 ,rotation = 90, color='tab:blue', fontsize=10) # M0.04
                ax[1].text(text_pos_co[dof_label][2], text_pos_co[dof_label][3], r'$f_{cou.,3}$ ' ,rotation = 90, color='tab:blue', fontsize=10)
                

            ax[0].locator_params(axis='y', nbins=5)
            yticks = {'Mz':[10e-3,1,10e1],'My':[10e-3,1,10e1],'Mx':[10e-5,10e-3,1,10e1],
                      'a_z':[10e-7,10e-5,10e-3,1],'a_y':[10e-7,10e-5,10e-3,1],'a_alpha':[10e-9,10e-7,10e-5,10e-3],'a_total':[10e-7,10e-5,10e-3,1]}
            ylim = {'Mz':[0.01,200],'My':[ 0.0001,200],'Mx':[1e-7,1],
                    'a_z':[1e-7,1],'a_y':[ 1e-7,1],'a_alpha':[1e-8,0.01],'a_total':[ 1e-7,1]}

            if not log:
                ax[1].locator_params(axis='y', nbins=4)
                ax[1].xaxis.set_major_locator(MultipleLocator(0.25)) #M:0.1
            
            else:
                if result_variable == 'acceleration':
                    ax[1].set_xlim(right = 50)# M: 10
                    ax[1].set_ylim(bottom=ylim[dof_label][0], top=ylim[dof_label][1])#
                    ax[1].set_yticks(yticks[dof_label])#,10e1
                else:
                    ax[1].set_xlim(right = 20)# M: 10
                    ax[1].set_ylim(bottom=ylim[dof_label][0], top=ylim[dof_label][1])#
                    ax[1].set_yticks(yticks[dof_label])#
            #ax[1].legend(loc= 'lower left')10e-5,,10e2
            # legend appears over the labels, thus adjusting the bottom of the plot box 
            #fig.subplots_adjust(bottom=0.2)
            #ax[1].minorticks_off()
            ax[1].yaxis.set_minor_locator(AutoMinorLocator())#,bottom=True)
            #ax[1].tick_params(axis='x', which='minor', right = False)#,bottom=True)
            #ax[1].xaxis.set_minor_locator(MultipleLocator(20))
            ax[1].grid()
        
        ax[0].grid()
        

        save_title =  dof_label + '_time_dyn_res_comp'
        if include_fft:
            save_title = dof_label + '_time_freq_dyn_res_comp'

        if self.savefig:
            plt.savefig(dest + os_sep + save_title)
            print('\nsaved:', dest + os_sep + save_title)

        if self.show_plots:
            plt.show()

        if add_fft:
            self.plot_fft(dof_label=dof_for_fft, dynamic_analysis=dynamic_analysis, init_dyn_res = init_res)    


    def plot_fft(self, dof_label, dynamic_analysis = None, given_series = None, sampling_freq = None, init_dyn_res = None):
        ''' 
        either give it:
            - a dynamic analysis object or 
            - directly a time series and the sample freqeuency
        
        dof_label: label of dof 
        ''' 
        fig = plt.figure(num='frequency_domain_result')#, figsize=(5,3)
        is_type = 'action '
        if dynamic_analysis:
            sampling_freq = 1/dynamic_analysis.dt
            dof = GD.DOF_LABELS['3D'].index(dof_label)
            time_series = dynamic_analysis.solver.dynamic_reaction[dof, :]
            given_series = time_series
            is_type = 'reaction '

        label2 = None
        if init_dyn_res:
            sampling_freq_init = 1/init_dyn_res.dt
            dof = GD.DOF_LABELS['3D'].index(dof_label)
            time_series_init = init_dyn_res.dynamic_reaction[dof, :]
            freq_half_init, series_fft_init = utils.get_fft(time_series_init, sampling_freq_init)
            plt.plot(freq_half_init, series_fft_init, label = 'uncoupled result', linestyle = '--', color = 'tab:orange')
            label2 = 'coupled result'


        freq_half, series_fft = utils.get_fft(given_series, sampling_freq)
        plt.plot(freq_half, series_fft, label = label2)


        plt.xlim(0.01,0.8)
        plt.ylabel('|Amplitude|')
        plt.xlabel('frequency [Hz]')
        plt.title(is_type + GD.DIRECTION_RESPONSE_MAP[dof_label] + ' in the frequency domain using FFT ')
        plt.grid(True)
        plt.legend()
        plt.show()

    def compare_stats(self, results, node_id, response_label, result_type, stats, unit, uncoupled_normed = False):
        ''' 
        results: list of dynamic_analysis objects 0: uncoupled, 1: coupled
        ''' 
        is_total = False
        if response_label == 'total':
            is_total = True

        dest = os_join(*['plots','CAARC_B','dynamic_results'])
        if result_type != 'acceleration':
            if response_label not in GD.DOF_LABELS['3D']:
                dof_label = GD.RESPONSE_DIRECTION_MAP[response_label]
            else:
                dof_label = response_label
                response_label = GD.DIRECTION_RESPONSE_MAP[response_label]
            dof_id = GD.DOF_LABELS['3D'].index(dof_label) + (node_id * GD.n_dofs_node['3D'])

        unit_scale = 1
        unit_label = ''
        if result_type == 'displacement':
            result_data = results[1].solver.displacement[dof_id, :]
            result_data_un = results[0].solver.displacement[dof_id, :]
        elif result_type == 'velocity':
            result_data = results[1].solver.velocity[dof_id, :]
            result_data_un = results[0].solver.velocity[dof_id, :]
        elif result_type == 'acceleration':
            unit_label = r'[m/s^{2}]'
            directions = {'total':{'dofs':['y','z','a'],'label':'a_total','rad':np.sqrt(15**2+22.5**2)},
                          'y':{'dofs':['y','a'],'label':'a_y','rad':22.5},'z':{'dofs':['z','a'],'label':'a_z','rad':15},
                          'a':{'dofs':['a'],'label':'a_alpha','rad':27}}
            dof_ids = [GD.DOF_LABELS['3D'].index(dof_i) + (node_id * GD.n_dofs_node['3D']) for dof_i in directions[response_label]['dofs']]
            result_data, result_data_un = np.zeros(results[1].array_time.size), np.zeros(results[1].array_time.size)
            for i, dof_id in enumerate(dof_ids):
                if GD.DOF_LABELS['3D'][dof_id - node_id * GD.n_dofs_node['3D']] == 'a':
                    rad_scale = directions[response_label]['rad']
                else:
                    rad_scale = 1.0
                result_data += (results[1].solver.acceleration[dof_id, :]*rad_scale)#**2
                result_data_un += (results[0].solver.acceleration[dof_id, :]*rad_scale)#**2
            # result_data = np.sqrt(result_data)
            # result_data_un = np.sqrt(result_data_un)
            response_label = directions[response_label]['label']
        elif result_type == 'action':
            result_data = results[1].force[dof_id, :]
            result_data_un = results[0].solver.force[dof_id, :]
        elif result_type == 'reaction':
            unit_label = GD.UNITS_POINT_LOAD_DIRECTION[dof_label]
            unit_label = unit_label.replace(unit_label[1], unit)
            unit_scale = GD.UNIT_SCALE[unit]


            stats = ['mean', 'std', 'max_est']
            # if dof_id in results[1].structure_model.dofs_of_bc:# or dof in results[1].structure_model.elastic_bc_dofs:
            #     result_data = results[1].solver.dynamic_reaction[dof_id, :] * unit_scale
            # else:
            print ('\nReplacing the selected node by the ground node for reaction results')
            #dof_id = 0
            dof_id = GD.DOF_LABELS['3D'].index(dof_label)
            result_data = results[1].solver.dynamic_reaction[dof_id, :]  * unit_scale
            
            result_data_un = results[0].solver.dynamic_reaction[dof_id, :] * unit_scale

        results_dict = {'uncoupled':{'time_hist':result_data_un,'stats':[]}, 'coupled':{'time_hist':result_data, 'stats':[]}}
        if result_type == 'acceleration' and not is_total:
            stats = ['std', 'max_est']
            if uncoupled_normed:
                stats = ['std', 'mean', response_label]

        for result in results_dict:
            time_series = results_dict[result]['time_hist']
            if result_type != 'acceleration':
                results_dict[result]['stats'].append(abs(np.mean(time_series)))
                print (np.mean(time_series))
            if result_type == 'acceleration' and is_total:#uncoupled_normed and
                results_dict[result]['stats'].append(abs(np.mean(time_series)))
            results_dict[result]['stats'].append(np.std(time_series))#+abs(np.mean(time_series)))
            results_dict[result]['stats'].append(utils.extreme_value_analysis_nist(time_series, 0.02)[1])

        directions_naming = {'Mz':'along wind', 'My':'cross wind','Mx':'torsion'}

        factor = 0.5
        colors = ['tab:blue', 'darkgray', 'gray',  'lightgray']
        fig, ax = plt.subplots(num='res compare')
        x = np.arange(0,len(results_dict['uncoupled']['stats']),1)* factor
        width = 0.15
        
        #if result_dict['norm_by'] == 'glob_max':
        norm_by = 1
        rot =30
        if uncoupled_normed:
            rot = 0
            norm_by = np.array([1/val for val in results_dict['uncoupled']['stats']])     #1/np.array(result_dict['dyn_max'])

        x_labels = [convert_for_latex(s) for s in stats]#results_dict['uncoupled']

        rects_un = ax.bar(x - width/2,  np.array(results_dict['uncoupled']['stats']) * norm_by, 
                        width, color = colors[1], label='uncoupled')
        rects_co = ax.bar(x + width/2, np.array(results_dict['coupled']['stats']) * norm_by, 
                        width, color = colors[0], label='coupled')
        dx = (x[1] - x[0])/2

        # ax.axvline(x[1]+dx, linestyle = '--', color = 'grey')
        # ax.axvline(x[3]+dx, linestyle = '--', color = 'grey')
        # ax.axhline(1,linestyle = '--', color = 'k',label= convert_for_latex('glob_max'))

        # for r_i, r in enumerate(['Mz','My','Mx']):
        #     ax.text(x[result_dict['labels'].index(r)],1.9, directions_naming[r])
        
        #ax.legend(bbox_to_anchor = (0.5, -0.4), loc ='upper center', ncol = 2)
        ax.legend()#loc= 'upper center', bbox_to_anchor=(0.5, -0.15), ncol= 3)
        if not uncoupled_normed:
            ax.set_ylabel(convert_for_latex(response_label) + r' ${}$'.format(unit_label))
        else:
            ax.set_ylabel(r'of uncoupled')
        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=rot)
        ax.set_ylim(0)
        if result_type == 'acceleration':
            if response_label == 'a_alpha':
                if not uncoupled_normed:
                    ax.yaxis.set_major_locator(MultipleLocator(0.03))#
            else:
                if not uncoupled_normed:
                    ax.yaxis.set_major_locator(MultipleLocator(0.5))
        else:
            ax.locator_params(axis='y', nbins=4)

        
        save_title = response_label + '_stats'

        if self.savefig:
            plt.savefig(dest + os_sep + save_title)
            print('\nsaved:', dest + os_sep + save_title)
        if self.show_plots:
            plt.show()

    def plot_compare_energies(self, values_dict):
        dest = os_join(*['plots','CAARC_B', 'dynamic_results'])
        fig = plt.figure(num='modal energy')
        norm = 1/values_dict['uncoupled']
        factor = 0.1
        hatches = ['//','o']
        colors = ['darkgray','tab:blue']
        for i, val in enumerate(values_dict.items()):
            #plt.bar(i*factor, val[1] * norm, label=val[0], width=0.1, color = 'w', hatch = hatches[i],edgecolor='tab:blue',) #
            plt.bar(i*factor, val[1] * norm, label=val[0], width=0.1, color = colors[i])#, hatch = hatches[i],edgecolor='tab:blue',) #
            #plt.bar(i*factor, val[1] * norm, label=val[0], width=0.1, color = 'w',  edgecolor = 'k') #hatches[i],edgecolor='tab:blue',
            plt.text(i*factor, val[1]*norm + 0.02,  r'${}$'.format(round(val[1]*norm,2)))#, fontsize =8)

        plt.xticks([])
        plt.yticks([])
        plt.legend(loc = 'lower center')
        label = r'$\sum^{T} E / \sum^{T} E_{uncoupled}$'#'sum of energy over time of uncoupled'
        plt.ylabel(label)

        save_tlt = 'sum_of_energy_new'
        if self.savefig:
            plt.savefig(dest + os_sep + save_tlt)
            print ('\nsaved:', dest + os_sep + save_tlt)
        if self.show_plots:
            plt.show()

    def set_ax_size(self, w,h, ax=None):
        """ w, h: width, height in inches """
        #if not ax: ax=plt.gca()
        l = ax.figure.subplotpars.left
        r = ax.figure.subplotpars.right
        t = ax.figure.subplotpars.top
        b = ax.figure.subplotpars.bottom
        figw = float(w)/(r-l)
        figh = float(h)/(t-b)
        ax.figure.set_size_inches(figw, figh)