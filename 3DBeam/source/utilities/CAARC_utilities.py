import numpy as np
from os.path import join as os_join
import matplotlib.pyplot as plt 
import source.postprocess as post
from source.utilities import utilities as utils
from source.utilities import global_definitions as GD

from plot_settings import plot_settings

COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

#params = plot_settings.get_params(w=7.3, h=5.2)


dest_folder = 'plots_new\\CAARC'

def get_CAARC_properties(src_file_name = 'CAARC_advanced_eigenmodes.txt', evaluate_at = None, interpolation_degree = 3):
    '''
    a dictionary is returned with information about the building A
        - storey (number), storey_level, mass, frequencies, eigenmodes
    the eigenmodes are appened to a list
    for each dof a list is created in a dictionary    
    '''
    src = os_join(*['inputs','eigenvectors', src_file_name])
    caarc = {}
    caarc['storey'] = np.flip(np.loadtxt(src, usecols = (0,))) # [-]
    caarc['storey_level'] = np.flip(np.loadtxt(src, usecols = (1,))) # [m]
    caarc['dimensons'] = {'x':240,'y':24, 'z':72}
    caarc['mass'] = 1231000.0
    caarc['frequencies'] = [0.231, 0.429, 0.536]
    caarc['eigenmodes'] = {'x':[],'y':[],'z':[],'a':[]}
    caarc['eigenmodes_fitted'] = {'x':[],'y':[],'z':[],'a':[]}

    for i in range (3):
        caarc['eigenmodes']['x'].append(np.zeros(60))
        caarc['eigenmodes']['y'].append(np.flip(np.loadtxt(src, usecols = (3+3*i,)))) 
        caarc['eigenmodes']['z'].append(np.flip(np.loadtxt(src, usecols = (2+3*i,)))) 
        caarc['eigenmodes']['a'].append(np.flip(np.loadtxt(src, usecols = (4+3*i,)))) 
  
    if not evaluate_at:
        x = caarc['storey_level']
    else:
        x = evaluate_at 

    for dof_label in ['y', 'z', 'a']:
        for mode_id in range(3):
            y = caarc['eigenmodes'][dof_label][mode_id]
            current_polynomial = np.poly1d(np.polyfit(caarc['storey_level'] ,y , interpolation_degree))
            values = []
            for x_i in x:# evaluate the fitted eigenmode at certain intervals
                values.append(current_polynomial(x_i))
            caarc['eigenmodes_fitted'][dof_label].append(np.asarray(values))

    return caarc

def get_CAARC_eigenform_polyfit (CAARC_eigenmodes, evaluate_at = None, degree = 5):
    '''
    retruns the values of a fitted caarc eigenmode.
    evaluate_at must be a list of x coordiantes at which the fitted curve should be evaluated.
    if it is not provided the fitted curve is evaluated at each storey level of caarc.
    '''
    eigenmodes_fitted = {} 
    #CAARC_eigenmodes = self.structure_model.CAARC_eigenmodes
    # returns the fitted polynomial and the discrete array of displacements
    if not evaluate_at:
        x = CAARC_eigenmodes['storey_level']
    else:
        x = evaluate_at 
    eigenmodes_fitted['storey_level'] = np.copy(x)
    eigenmodes_fitted['eigenmodes'] = []

    for mode_id in range(1,4):
        eigenmodes_fitted['eigenmodes'].append({})
        for dof_label in ['y', 'z', 'a']:
            y = CAARC_eigenmodes['eigenmodes'][mode_id][dof_label]
            current_polynomial = np.poly1d(np.polyfit(CAARC_eigenmodes['storey_level'],y,degree))
            values = []
            for x_i in x:# evaluate the fitted eigenmode at certain intervals
                values.append(current_polynomial(x_i))
            eigenmodes_fitted['eigenmodes'][mode_id][dof_label] = np.asarray(values)

    return eigenmodes_fitted



def get_m_eff(eigenmodes_dict, mode_id, main_direction_only, print_to_console):
    '''
    retruns the generalized mass and the participation factor of a mode 
    prints the effective mass that should be around 60% of the total mass (first modes)
    '''

    mass = eigenmodes_dict['mass'] # constant over height
    phi_y = eigenmodes_dict['eigenmodes']['y'][mode_id]
    phi_z = eigenmodes_dict['eigenmodes']['z'][mode_id]

    if main_direction_only:
        if mode_id == 1:
            participation_factor = (mass * sum(phi_y))**2 # mass not in the sum since it is constant
        elif mode_id == 2:
            participation_factor = (mass * sum(phi_z))**2
    else:
        participation_factor = (mass * sum(np.add(phi_y, phi_z)))**2


    if main_direction_only:
        if mode_id == 1:
            generalized_mass = mass * sum(np.square(phi_y))
        elif mode_id == 2:
            generalized_mass = mass * sum(np.square(phi_z))
    else:
        generalized_mass = mass * sum(np.add(np.square(phi_y), np.square(phi_z)))

    total_mass = 60*mass
    m_eff = participation_factor/generalized_mass 

    if print_to_console:
        print('m_eff of mode_id', mode_id,  ':', round(m_eff/total_mass, 4), 'of m_tot')  
        print ('generalized_mass:', round(generalized_mass, 2), 'should be 1 if mass normalized')     
        print ('participation_factor:', round(participation_factor,2)) 
        print () 

    return participation_factor, generalized_mass  

def plot_caarc_eigenmodes(caarc_dict,  number_of_modes = 3, dofs_to_plot = ['y','a','z'],
                            max_normed = False, do_rad_scale =False, use_caarc_fitted = False,
                            savefig = False, savefig_latex = False,
                            fig_title = '', filename_for_save = '0_no_name'):

    c_norm = 1
    rad_scale = np.sqrt(caarc_dict['dimensons']['y'] *caarc_dict['dimensons']['z'])
    if max_normed:
        do_rad_scale = False
    
    if use_caarc_fitted:
        c_modes = caarc_dict['eigenmodes_fitted']
    else:
        c_modes = caarc_dict['eigenmodes']
    

    if number_of_modes == 1:
        raise Exception('for 1 mode not implemented yet')
        fig, ax = plt.subplots(ncols = number_of_modes,  num='eigenmode results')#figsize=(2,3.5),
        if not self.savefig and not self.savefig_latex:
            fig.suptitle(fig_title)

        x = beam_model.nodal_coordinates['x0']
        ax.plot( beam_model.nodal_coordinates['y0'],
                    x,
                    #label = r'$structure$',
                    color = 'grey',
                    linestyle = '--')
        ax.set_title(r'$mode\,1$ '+'\n' +r'$f = $ ' + r'${}$'.format(str(round(beam_model.eigenfrequencies[0],3))) + r' $Hz$'  + weights)

        for d_i, dof in enumerate(dofs_to_plot):
            scale=1.0
            if do_rad_scale:
                if dof in ['a','b','g']:
                    scale = rad_scale
                else:
                    scale = 1.0
            y = utils.check_and_flip_sign_array(beam_model.eigenmodes[dof][0])
            if include_caarc:
                c_y = utils.check_and_flip_sign_array(c_modes[dof][0])

            if max_normed:
                norm = 1/max(y)
                c_norm = 1/max(c_y)
            
            if opt_targets:
                if dof in opt_targets.keys():
                    y2 = opt_targets[dof] *scale
                    ax.plot(y2*norm,
                                x,
                                label =  r'${}$'.format(GD.greek[dof]) + r'$_{target}$',
                                linestyle = '--',
                                color = COLORS[d_i])
            if initial:
                if dof in initial.keys():
                    y3 = initial[dof]
                    ax.plot(y3*norm*scale,
                                x,
                                label =  r'${}$'.format(GD.greek[dof]) + r'$_{inital}$',
                                linestyle = ':',
                                color = COLORS[d_i])
            lab = GD.greek[dof]
            ax.plot(y*norm*scale,
                        x,
                        label = r'${}$'.format(GD.greek[dof]) + r'$_{max}:$' + r'${0:.2e}$'.format(max(y)),#'max: ' +
                        linestyle = '-',
                        color = COLORS[d_i])
            if include_caarc:
                ax.plot(c_y*c_norm*scale,
                        c_x,
                        label = r'${}$'.format(GD.greek[dof]) + r'$benchmark$',#'max: ' +
                        linestyle = ':',
                        color = caarc_cols[d_i])
        ax.legend()
        ax.grid()
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r'$deflection$')
        ax.set_ylabel(r'$x \, [m]$') 

        ratio = max(utils.check_and_flip_sign_array(beam_model.eigenmodes['a'][0])) / max(utils.check_and_flip_sign_array(beam_model.eigenmodes['y'][0]))
        ax.plot(0,0, label = r'$\alpha_{max}/y_{max}: $' + str(round(ratio,3)), linestyle = 'None')    
        ax.legend()

    else:
        fig, ax = plt.subplots(ncols = number_of_modes, sharey=True,  num='eigenmode results')#figsize=(5,4),
        
        for i in range(number_of_modes):
            x = caarc_dict['storey_level']
            ax[i].plot( np.zeros(len(x)),
                        x,
                        #label = r'$structure$',
                        color = 'grey',
                        linestyle = '--')
            ax[i].set_title(r'$mode$ ' + r'${}$'.format(str(i+1)) + '\n' +r'$f=$ ' + r'${}$'.format(str(round(caarc_dict['frequencies'][i],3))) +r' $Hz$')

            
            for d_i, dof in enumerate(dofs_to_plot):
                scale=1.0
                if do_rad_scale:
                    if dof in ['a','b','g']:
                        scale = rad_scale
                    else:
                        scale = 1.0
                    
                y = utils.check_and_flip_sign_array(c_modes[dof][i])

                if max_normed:
                    c_norm = 1/max(y)
                
                ax[i].plot(y*c_norm*scale,
                            x,
                            label = r'${}$'.format(GD.greek[dof]),# + r'$_{max}:$ ' + r'${0:.2e}$'.format(max(y)),
                            linestyle = '-',
                            color = COLORS[d_i])

            ax[i].legend(loc = 'lower right')
            ax[i].grid()
            ax[i].set_ylim(bottom=0)
            ax[i].set_xlabel(r'$deflection$')
            ax[0].set_ylabel(r'$x \, [m]$') 

        ratio = max(utils.check_and_flip_sign_array(caarc_dict['eigenmodes']['a'][0])) / max(utils.check_and_flip_sign_array(caarc_dict['eigenmodes']['y'][0]))
        #ax[0].plot(0,0, label = r'$\alpha_{max}/y_{max} = $' + str(round(ratio,3)), linestyle = 'None')    
        ax[0].legend(loc= 'lower right')
        #plt.tight_layout()
    
    #if self.show_plots:
    #plt.grid()
    plt.show()

    plt.close()


if __name__ == '__main__':
    plt.rcParams.update({'axes.formatter.limits':(-3,3)}) 
    plt.rcParams.update(params)
    caarc_dict = get_CAARC_properties()
    plot_caarc_eigenmodes(caarc_dict, do_rad_scale=True,
                            savefig=1, savefig_latex=1,
                            filename_for_save='caarc_modes')
