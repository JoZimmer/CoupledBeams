import string
import numpy as np
import json
import os
import matplotlib.pyplot as plt 
from os.path import join as os_join
from os.path import sep as os_sep
from source.utilities import statistics_utilities as stats_utils
from source.utilities import global_definitions as GD


caarc_freqs = [0.231, 0.429, 0.536]
# VALUES COPIED WITH FULL MATRICIES CALCULATED
eigenmodes_target_2D_3_elems = {'y':[np.array([0.        , 0.0040305 , 0.013317  , 0.02434817]),
                                    np.array([ 0.        , -0.01444802, -0.01037197,  0.02449357]),
                                    np.array([ 0.        , -0.01819986,  0.01604842, -0.02446053])],
                                'g':[np.array([ 0.        , -0.0004894 , -0.00070823, -0.00074479]),
                                    np.array([ 0.        ,  0.00095991, -0.00161083, -0.00260447]),
                                    np.array([ 0.        , -0.0009034 , -0.00068957,  0.00432069])]}

eigenmodes_target_2D_10_elems ={'y':[np.array([0.0, -0.000223647046255273, -0.0008516138734941783, -0.0018197756020471587, -0.0030651302400258222, -0.00452698257707041, -0.006148471228742031, -0.007878364112695452, -0.009673052432583382, -0.011498680245678633, -0.01333335614230312]),
                                    np.array([ 0.        ,  0.00123514,  0.00401433,  0.00701557,  0.00911353, 0.00951618,  0.0078602 ,  0.00422764, -0.00093387, -0.00698381, -0.01333421]),
                                    np.array([ 0.        ,  0.00304246,  0.00806418,  0.01008837,  0.007016  , 0.00026276, -0.00632   , -0.00877015, -0.00526778,  0.00304816, 0.01334002])],
                                'g':[np.array([0.00000000e+00, 2.91028048e-05, 5.39127464e-05, 7.44741342e-05,9.08970510e-05, 1.03382795e-04, 1.12244221e-04, 1.17921246e-04,1.20991889e-04, 1.22179408e-04, 1.22356253e-04]),
                                    np.array([ 0.00000000e+00, -1.49126242e-04, -2.06595585e-04, -1.80905410e-04,-8.99100753e-05,  4.02818206e-05,  1.79513859e-04,  2.99658399e-04,3.81148479e-04,  4.18652321e-04,  4.24986410e-04]),
                                    np.array([ 0.00000000e+00, -3.34883542e-04, -2.77308592e-04,  3.15745412e-05,3.61060122e-04,  4.93762038e-04,  3.37172087e-04, -3.17237701e-05,-4.21134643e-04, -6.52640493e-04, -6.98021127e-04])]} 

def evaluate_residual(a_cur, a_tar):
    residual = np.linalg.norm(np.subtract(a_cur, a_tar)) / np.amax(np.absolute(a_tar))

    # print ('current: ', a_cur)
    # print ('target: ', a_tar)
    # print('residual:', residual)
    # print()
    return residual

def cm2inch(value):
    return value/2.54
    
def increasing_by(val_old, val_new):
    ''' 
    returns the increase in % from the origin = old
    ''' 
    increase = (val_new - val_old)/val_old * 100
    return round(increase,2)

def check_and_flip_sign_dict(eigenmodes_dict):
    '''
    flips the sign of y and a deformation of modes to be positive at the first node after ground 
    dependend/coupled dofs are flip accordingly
    '''

    for idx in range(len(eigenmodes_dict['y'])):
        if eigenmodes_dict['y'][idx][1] < 0:
            eigenmodes_dict['y'][idx] = np.negative(eigenmodes_dict['y'][idx])
            try:
                eigenmodes_dict['g'][idx] = np.negative(eigenmodes_dict['g'][idx])
            except KeyError:
                pass
        try:
            if eigenmodes_dict['a'][idx][1] < 0:
                eigenmodes_dict['a'][idx] = np.negative(eigenmodes_dict['a'][idx])
        except KeyError:
            pass

    return eigenmodes_dict

def check_and_flip_sign_array(mode_shape_array):
    '''
    check_and_change_sign
    change the sign of the mode shape such that the first entry is positive
    '''
    if mode_shape_array[1] < 0:
        return mode_shape_array * -1
    elif mode_shape_array [1] > 0:
        return mode_shape_array

def analytic_function_static_disp(parameters, x, load_type = 'single'):
    l = parameters['lx_total_beam']
    EI = parameters['E_Modul'] * parameters['Iy']
    magnitude = parameters['static_load_magnitude']
    if load_type == 'single':
        #print ('  w_max soll:', magnitude*l**3/(3*EI))
        return (magnitude/EI) * (0.5 * l * (x**2) - (x**3)/6) 
    elif load_type == 'continous':
        #print ('  w_max soll:', magnitude*l**4/(8*EI))
        return -(magnitude/EI) * (-x**4/24 + l * x**3 /6 - l**2 * x**2 /4)

def analytic_eigenfrequencies(beam_model):
    # von https://me-lrt.de/eigenfrequenzen-eigenformen-beim-balken 
    parameters = beam_model.parameters
    l = parameters['lx_total_beam']
    EI = parameters['E_Modul'] * parameters['Iy']
    A = parameters['cross_section_area']
    rho = parameters['material_density']
    lambdas = [1.875, 4.694, 7.855]
    f_j = np.zeros(3)
    for i, l_i in enumerate(lambdas):
        f_j[i] = np.sqrt((l_i**4 * EI)/(l**4 * rho *A)) / (2*np.pi)
    
    return f_j

def analytic_eigenmode_shapes(beam_model):
    # von https://me-lrt.de/eigenfrequenzen-eigenformen-beim-balken 
    #parameters = beam_model.parameters
    l = beam_model.parameters['lx_total_beam']
    x = beam_model.nodal_coordinates['x0']
    lambdas = [1.875, 4.694, 7.855] # could also be computed as seen in the link

    w = []
    for j in range(3):
        zeta = lambdas[j] * x / l
        a = (np.sin(lambdas[j]) - np.sinh(lambdas[j])) / (np.cos(lambdas[j]) + np.cosh(lambdas[j])) 
        
        w_j = np.cos(zeta) - np.cosh(zeta) + a*(np.sin(zeta) -np.sinh(zeta))
        w.append(w_j)
        
        reduced_m = beam_model.m[::2,::2] 
        gen_mass = np.dot(np.dot(w_j, reduced_m), w_j)
        norm_fac = np.sqrt(gen_mass)
        w_j /= norm_fac
        is_unity = np.dot(np.dot(w_j, reduced_m), w_j)
        if round(is_unity, 4) != 1.0:
            raise Exception ('analytic mode shape normalization failed')

    return w

def get_eigenform_polyfit(modeshape_i, z_coords,  evaluate_at, degree = 5, plot_compare = False):
    ''' 
    - modeshape_i: all modal deformations als 2D array, each column belongs to one dof
    - z_coords: the original floor levels from the 45 floor model
    - evaluate_at: nodal coordiantes at whcih the fitted curve should be evaluated 
                   -> this is returned
    ''' 
    eigenmodes_fitted = {} 
    #CAARC_eigenmodes = self.structure_model.CAARC_eigenmodes
    # returns the fitted polynomial and the discrete array of displacements
    if not evaluate_at.any():
        raise Exception('provied evaluation coordiantes of the eigenform')
    else:
        x = evaluate_at 

    eigenmodes_fitted['storey_level'] = np.copy(x)
    eigenmodes_fitted['eigenmodes'] = {}

    dof_direction_map = {'y':4, 'z':3,'a':2}

    for dof_label in ['y', 'z', 'a']:
        y = modeshape_i[:, dof_direction_map[dof_label]]
        current_polynomial = np.poly1d(np.polyfit(z_coords,y,degree))
        values = []
        for x_i in x:# evaluate the fitted eigenmode at certain intervals
            values.append(current_polynomial(x_i))
        if values[0] != 0.0:
            values[0] = 0.0
        eigenmodes_fitted['eigenmodes'][dof_label] = np.asarray(values)

    if plot_compare:
        fig, ax = plt.subplots(ncols=3, num='fitted compared')
        for d_i, dof_label in enumerate(['y', 'z', 'a']):
            ax[d_i].plot(modeshape_i[:, dof_direction_map[dof_label]], z_coords, label = 'origin ' + dof_label)
            ax[d_i].plot(eigenmodes_fitted['eigenmodes'][dof_label], x, label = 'fitted ' + dof_label)
            ax[d_i].legend()
            ax[d_i].grid()

        plt.show()

    return eigenmodes_fitted

def save_optimized_beam_parameters(opt_beam_model, fname):
    new = 'optimized_parameters'+os_sep+fname+'.json'
    if os.path.isfile(new):
        print('WARNING', new, 'already exists!')
        new = 'optimized_parameters'+os_sep+fname+'_1.json'

    f = open(new,'w')
    
    json.dump(opt_beam_model.parameters, f)
    f.close()
    print('\nsaved:', new)

def get_targets(beam_model, target='semi_realistic', opt_params =None):
    ''' 
    just used to have them especially for the initial comparison
    ''' 
    if target == 'realistic':
        modi = np.load(os_join(*['inputs', 'EigenvectorsGid.npy']))
        z_coords = np.load(os_join(*['inputs', 'z_coords_gid_45.npy']))
        
        modi_fitted = get_eigenform_polyfit(modi[0], z_coords, beam_model.nodal_coordinates['x0'], plot_compare=False)
        eigenmodes_target_y = modi_fitted['eigenmodes']['y']
        eigenmodes_target_a = -1*modi_fitted['eigenmodes']['a']

    elif target == 'semi_realistic':
        ratio_a_y = opt_params['ratio_a_y_tar']
        factor_y = opt_params['factor_y']
        eigenmodes_target_y = beam_model.eigenmodes['y'][0]*factor_y

        a_factor = ratio_a_y * max(eigenmodes_target_y)/max(beam_model.eigenmodes['a'][0])
        eigenmodes_target_a = beam_model.eigenmodes['a'][0] * a_factor

    return {'y':eigenmodes_target_y, 'a':eigenmodes_target_a}

def prepare_string_for_latex(string):
    greek = {'ya':'y'+r'\alpha','ga':r'\gamma' + r'\alpha'}
    if '_' in string:
        var, label = string.split('_')[0], string.split('_')[1]
        latex = r'${}$'.format(var) + r'$_{{{}}}$'.format(greek[label])
        #return string.replace('_','')
        return latex
    else:
        return string

    
def join_whitespaced_string(string):
    return string.replace(' ','_')

def parse_load_signal_backwards(signal):
    '''
    signal kommt als dictionary mit den Richtungen als keys (müssen nicht alle 6 Richtugne sein)
    output soll row vector sein mit dofs * n_nodes einträgen
    '''
    shape = GD.DOFS_PER_NODE['3D'] * len(list(signal.values())[0])
    signal_raw = np.zeros(shape)
    for label in signal:
        dof_label_id = GD.DOF_LABELS['3D'].index(label)
        for i, val in enumerate(signal[label]):
            sort_id = i * GD.DOFS_PER_NODE['3D'] + dof_label_id
            signal_raw[sort_id] = val

    #TODO: must this be transposed?!
    return signal_raw

def generate_static_load_vector_file(load_vector):
    '''
    return a npy file and saves it for a given load.

    If load is given as dictionary with directions as keys it parses it to a 1D array that is needed in the ParOptBeam
    '''
    if isinstance(load_vector, dict):
        load = parse_load_signal_backwards(load_vector)
    else:
        load = load_vector

    src_path = os_join('input','loads','static')
    if not os.path.isdir(src_path):
        os.makedirs(src_path)

    file_name = src_path + 'eswl_' + str(int(len(load)/GD.DOFS_PER_NODE['3D'])) + '_nodes.npy'

    np.save(file_name, load)

    return file_name

def generate_nodal_force_file(number_of_nodes, node_of_load_application, force_direction:str, magnitude, domain_size = '3D', file_base_name = 'static_load_'):
    '''
    creating a force .npy file with a nodal force at given node, direction and magnitude
    number_of_nodes: anzahl knoten des models 
    node_of_load_application: Knoten an dem die Last aufgebracht werden soll
    force_direction: string der richtung
    '''
    src_path = os_join(*['input','loads','static'])
    if not os.path.isdir(src_path):
        os.makedirs(src_path)

    # -1 ist notwendig da sozusagen vom knoten unter dem von interesse angefangen wird die dof nodes dazu zu addieren
    loaded_dof = (node_of_load_application-1)*GD.DOFS_PER_NODE[domain_size] + GD.DOF_LABELS[domain_size].index(force_direction)

    force_data = np.zeros(GD.DOFS_PER_NODE[domain_size]*number_of_nodes)
    force_data[loaded_dof] += magnitude
    force_data = force_data.reshape(GD.DOFS_PER_NODE[domain_size]*number_of_nodes,1)

    force_file_name = os_join(src_path, file_base_name + str(number_of_nodes) + '_nodes_at_' + str(node_of_load_application) + \
                        '_in_' + force_direction+'.npy')
    np.save(force_file_name, force_data)

    return force_file_name

# # DYNAMIC ANALYSIS
def get_fft(given_series, sampling_freq):

    '''
    The function get_fft estimates the Fast Fourier transform of the given signal 
    sampling_freq = 1/dt
    '''

    signal_length=len(given_series)

    freq_half =  np.arange(0, 
                           sampling_freq/2 - sampling_freq/signal_length + sampling_freq/signal_length, 
                           sampling_freq/signal_length)

    # single sided fourier
    series_fft = np.fft.fft(given_series)
    series_fft = np.abs(series_fft[0:int(np.floor(signal_length/2))])/np.floor(signal_length/2)  
    
    max_length = len(freq_half)
    if max_length < len(series_fft):
        max_length = len(series_fft)
    
    freq_half = freq_half[:max_length-1]
    series_fft = series_fft[:max_length-1]
    
    return freq_half, series_fft


def extreme_value_analysis_nist(given_series, dt, response_label = None, type_of_return = 'estimate', P1 = 0.98):
    ''' 
    dynamic_analysi_solved: dynamic_analysis.solver object
    response: label given as dof_label, if given as response label it is convertedd
    type_of_return: wheter the estimated or the quantile value of P1 is returned (both are computed)
    ''' 

    T_series = dt * len(given_series)
    dur_ratio = 600 / T_series
    # # MAXMINEST NIST
    #P1 = 0.98
    max_qnt, min_qnt, max_est, min_est, max_std, min_std, Nupcross = stats_utils.maxmin_qnt_est(given_series	, 
                                                                        cdf_p_max = P1 , cdf_p_min = 0.0001, cdf_qnt = P1, dur_ratio = dur_ratio)
    
    abs_max_est = max([abs(max_est[0][0]), abs(min_est[0][0])])
    abs_max_qnt = max([abs(max_qnt[0]), abs(min_qnt[0])])

    if type_of_return == 'estimate':
        extreme_response = abs_max_est
    elif type_of_return == 'quantile':
        extreme_response = abs_max_qnt
    
    glob_max = max(abs(given_series))

    return abs_max_qnt, abs_max_est

# SONSTIGES

def read_xl_column(xl_worksheet, start_cell:str = None, start_row:int = None, start_col:int = None, end_row:int=None, end_cell:str = None):
    '''
    reading excel colum startin from start_row, start_col ODER start_cell als Excel Zellen name
    start_row = Excel Zeilen Zahl (wird intern dann angepasst auf das 0 zählsystem von xlwings)
    start_col = A entspricht 0 
    end_row_optional: wenn nicht gegebn wird spalte bis zur nächsten leeren Zeile gelesen NOTE klappt nicht immer 
    '''
    def excel_cell_name_to_row_col(cell_name):

        row = int(cell_name[1:]) - 1 
        col = string.ascii_uppercase.index(cell_name[0])
        return row, col

    if start_cell:
        if start_col or start_row:
            raise Exception('entweder Start Zelle als String oder Zeile und Spalte als nummern geben')

        start_row, start_col = excel_cell_name_to_row_col(start_cell)

    else:
        start_row -=1
         
    if not end_row and not end_cell:
        end_row = xl_worksheet.range(start_row, start_col).end('down').row

    if end_cell:
        if end_row:
            raise Exception('End Zelle als String oder Zeilen nummer angeben')

        end_row, end_col = excel_cell_name_to_row_col(end_cell)
    
    column_values = xl_worksheet[start_row:end_row, start_col].value
    return column_values

