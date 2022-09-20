import pickle
import numpy as np
import os
from os.path import join as os_join
from os.path import sep as os_sep
import string 
import xlwings as xl
import pandas as pd
import copy

# from source.utilities import statistics_utilities as stats_utils
from source.utilities import global_definitions as GD

#_____________ ALTES VOM URSPRÜNGLICHEN BEAM_____________________

def evaluate_residual(a_cur, a_tar):
    residual = np.linalg.norm(np.subtract(a_cur, a_tar)) / np.amax(np.absolute(a_tar))

    # print ('current: ', a_cur)
    # print ('target: ', a_tar)
    # print('residual:', residual)
    # print()
    return residual

def cm2inch(value):
    return value/2.54

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
    EI = parameters['E_Modul'] * parameters['Iz']
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
    EI = parameters['E_Modul'] * parameters['z']
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

def save_optimized_beam_parameters(opt_beam_model, fname):
    import json
    new = 'optimized_parameters'+os_sep+fname+'.json'
    if os.path.isfile(new):
        print('WARNING', new, 'already exists!')
        new = 'optimized_parameters'+os_sep+fname+'_1.json'

    f = open(new,'w')
    
    json.dump(opt_beam_model.parameters, f)
    f.close()
    print('\nsaved:', new)

def remove_small_values_from_array(arr, tolerance=1e-05):

    for i, val in enumerate(arr):
        if abs(val) <= tolerance:
            arr[i] = 0
    return arr

#________________ LASTEN ____________________________________

def parse_load_signal_backwards(signal, domain_size):
    '''
    signal kommt als dictionary mit den Richtungen als keys (müssen nicht alle 6 Richtugne sein)
    output soll row vector sein mit dofs * n_nodes einträgen
    '''
    from source.utilities import global_definitions as GD

    shape = GD.DOFS_PER_NODE[domain_size] * len(list(signal.values())[0])
    signal_raw = np.zeros(shape)
    for label in signal:
        dof_label_id = GD.DOF_LABELS[domain_size].index(label)
        for i, val in enumerate(signal[label]):
            sort_id = i * GD.DOFS_PER_NODE[domain_size] + dof_label_id
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
    Einheit der Kraft = [N]
    '''
    from source.utilities import global_definitions as GD
    src_path = os_join(*['inputs','loads','static'])
    if not os.path.isdir(src_path):
        os.makedirs(src_path)

    force_file_name = os_join(src_path, file_base_name + str(number_of_nodes) + '_nodes_at_' + str(node_of_load_application) + \
                        '_in_' + force_direction+'.npy')

    if os.path.isfile(force_file_name):
        return force_file_name
    
    else:
        # -1 ist notwendig da sozusagen vom knoten unter dem von interesse angefangen wird die dof nodes dazu zu addieren
        loaded_dof = (node_of_load_application-1)*GD.DOFS_PER_NODE[domain_size] + GD.DOF_LABELS[domain_size].index(force_direction)

        force_data = np.zeros(GD.DOFS_PER_NODE[domain_size]*number_of_nodes)
        force_data[loaded_dof] += magnitude
        force_data = force_data.reshape(GD.DOFS_PER_NODE[domain_size]*number_of_nodes,1)

        np.save(force_file_name, force_data)

        return force_file_name

def generate_kopflasten_file(number_of_nodes, loads_dict, file_base_name = 'kopflasten_'):
    '''
    creating a force .npy file with a nodal force at given node, direction and magnitude
    number_of_nodes: anzahl knoten des models 
    node_of_load_application: Knoten an dem die Last aufgebracht werden soll
    force_direction: string der richtung
    Einheit der Kraft = [N]
    '''
    from source.utilities import global_definitions as GD
    src_path = os_join(*['inputs','loads','static'])
    if not os.path.isdir(src_path):
        os.makedirs(src_path)

    node_of_load_application = number_of_nodes
    force_file_name = os_join(src_path, file_base_name + '.npy')

    # if os.path.isfile(force_file_name):
    #     return force_file_name
    
    #else:
    # -1 ist notwendig da sozusagen vom knoten unter dem von interesse angefangen wird die dof nodes dazu zu addieren
    force_data = np.zeros(GD.DOFS_PER_NODE['2D']*number_of_nodes)
    for load_direction, magnitude in loads_dict.items():

        force_direction = GD.LOAD_DOF_MAP[load_direction]
        loaded_dof = (node_of_load_application-1)*GD.DOFS_PER_NODE['2D'] + GD.DOF_LABELS['2D'].index(force_direction)
        force_data[loaded_dof] += magnitude

    force_data = force_data.reshape(GD.DOFS_PER_NODE['2D']*number_of_nodes,1)

    np.save(force_file_name, force_data)

    return force_file_name

def update_lasten_dict(lasten_dict_base, wind_kraft, kopflasten):
    '''
    merged die Kopflasten und die windkraft nach DIN 
    lasten_dict_base: gibt größe und bezeichnugn der kräfte vor
    wenn windkraft = None wird nur die Kopflast in das richtige format gebracht
    '''
    lasten_dict_updated = copy.deepcopy(lasten_dict_base)
    if isinstance(wind_kraft, np.ndarray):
        lasten_dict_updated['Fy'] += wind_kraft
    for i, direction in enumerate(lasten_dict_base):
        lasten_dict_updated[direction][-1] += kopflasten[direction]
    return lasten_dict_updated

def generate_lasten_file(number_of_nodes, lasten_dict, file_base_name):

    from source.utilities import global_definitions as GD
    src_path = os_join(*['inputs','loads','static'])
    if not os.path.isdir(src_path):
        os.makedirs(src_path)

    force_file_name = os_join(src_path, file_base_name + '.npy')

    force_data = np.zeros(GD.DOFS_PER_NODE['2D']*number_of_nodes)
    for load_direction, magnitude in lasten_dict.items():
        force_direction = GD.LOAD_DOF_MAP[load_direction]
        start = GD.DOF_LABELS['2D'].index(force_direction)
        step = len(GD.DOF_LABELS['2D'])

        sub = force_data[start::step]

        force_data[start::step] += magnitude

    force_data = force_data.reshape(GD.DOFS_PER_NODE['2D']*number_of_nodes,1)

    np.save(force_file_name, force_data)

    return force_file_name

def linien_last_to_knoten_last(last_vektor, knoten_z, gamma_q=1.5):

    '''
    vekotr der eine Last pro meter darstellt in knotenlasten vektor überführen 
    '''
    if not isinstance(knoten_z, np.ndarray):
        knoten_z = np.array(knoten_z)

    d_i = np.diff(knoten_z)
    
    d_i[0] = d_i[1]/2
    d_i = np.append(d_i, d_i[1]/2)

    F_z = gamma_q * last_vektor * np.array(d_i)
    return F_z

def parse_schnittgrößen_labels(lasten_dict):
    '''
    lasten dict mit ersten keys gleich den einwirkungsdauern
    '''
    lasten_dict_parsed = copy.deepcopy(lasten_dict)
    for dauer, lasten in lasten_dict.items():
        for direction in lasten:
            lasten_dict_parsed[dauer][GD.DOF_RESPONSE_MAP[direction]] = lasten[direction]
            del lasten_dict_parsed[dauer][direction]

    return lasten_dict_parsed


#____________________DYNAMIC ANALYSIS________________________________
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

#______________________SONSTIGES_____________________________________

def load_object_from_pkl(pkl_file):
    '''full path mit datei name zum jobejt file'''
    with open(pkl_file, 'rb') as handle:
        data = pickle.load(handle)

    return data

def add_model_data_from_dict(section_dict, model_parameters,):
    '''
    Die Sektionsweiße definierten Querschnittswerte werden so an das Model gegebn.

    - Knoten Anzahl wird an die Anzahl an Ebenen angepasst
    - Koordinaten definition wird hier an den Beam angepasst
    - Jedem Interval wird der Mittelwert der jeweiligen Querschnittswerte zugeordnet 
        (dieser kommt schon aus dem QS dictionary )
    '''
    
    model_parameters['defined_on_intervals'] = []

    I = 'Iz'

    model_parameters['lx_total_beam'] = section_dict['section_absolute_heights'][-1]
    model_parameters['n_elements'] = section_dict['n_sections']

    for section in range(section_dict['n_sections']):
        model_parameters['defined_on_intervals'].append(
            {
            'interval_bounds':[section_dict['section_absolute_heights'][section], section_dict['section_absolute_heights'][section+1]],
            'area': [section_dict['A'][section]],
            'Iz':[section_dict[I][section]], # anpassen der Koordianten Richtung
            'D':[section_dict['d_achse'][section]]
            }
            )
        if section == section_dict['n_sections']:
            model_parameters['defined_on_intervals']['interval_bounds'][-1] = 'End'

    return model_parameters

#_____________________OPENFAST / RFEM_______________________________________

def convert_coordinate_system_and_consider_einwirkungsdauer(loads, direction = 'fore-aft'):
    '''
    Lasten aus IEA, OpenFAST in die Richtung vom Balken konvertieren
    Sortiert die Lasten nach einwirkungsdauer: Annahme -> Vertikallast = ständig die anderen kurz
    Dauer 'egal' gibt alle lasten zurück
    '''
    einwirkungsdauer = {'ständig':{'Fy':0, 'Fx':1, 'Mz':0}, 'kurz':{'Fy':1, 'Fx':0, 'Mz':1}, 'egal':{'Fy':1, 'Fx':1, 'Mz':1}}
    directions_beam = ['Fy', 'Fx', 'Mz']
    loads_beam = {'ständig':{}, 'kurz':{}, 'egal':{}}
    if isinstance(loads, dict):
        converter_beam_fast = {'Fy':'Fx','Fx':'Fz', 'Mz':'My'}
        # NOTE rein theoretisch ist dieser Vorzeichen wechsel hier falsch 
        sign_beam_fast = {'Fy':1,'Fx':1, 'Mz':-1}
        for dauer in einwirkungsdauer:
            for direction in directions_beam:
                factor = einwirkungsdauer[dauer][direction]
                loads_beam[dauer][direction] = loads[converter_beam_fast[direction]] * factor * sign_beam_fast[direction]

    return loads_beam

def teilsicherheitsbeiwert_gravity_IEC(DLC, F_gravity, Fk):
    '''
    nach IEC 61400-1 Kapitel 7.6.2.2
    '''
    if DLC == 1.1:
        phi = 0.15
    else:
        phi = 0.25

    if abs(F_gravity) > abs(Fk):
        zeta = 0
    else:
        zeta = 1 - abs(F_gravity/Fk)

    gamma_f = 1.1 + phi*zeta**2
    return gamma_f


#____________________  ALLGEMEINES __________________________________________________

def convert_for_latex(string):
    l = list(string.split('_'))
    label = r'${}$'.format(l[0])
    for i in range(1,len(l)):
        label += r'$_{{{}}}$'.format(l[i])

    return label

def unit_conversion(unit_in:str, unit_out:str) -> float:
    '''
    gibt faktor zurück um einen Wert der Einheit unit_in in einen der Einheit unit_out umzurechnen
    '''
    converter = {'N/mm²':{'N/cm²':1E+02, 'N/m²':1E+06, 'kN/mm²':1E-03,'kN/m²':1E+03,'MPa':1},
                'N/m²':{'N/mm²':1E-06,'N/cm²':1E-04,'kN/mm²':1E-09,'kN/m²':1E-03,'Pa':1, 'MPa':1E-06},
                'N/cm²':{'N/mm²':1E-02},
                'kN/m²':{'N/mm²':1E-03,'N/cm²':1E-01,'kN/mm²':1E-06,'kN/m²':1,'MPa':1E-03},
                'N':{'kN':1E-03, 'MN':1E-06},
                'kN':{'N':1E+03, 'MN':1E-03},
                'MN':{'N':1E+06, 'kN':1E+03},
                'Nm':{'kNm':1E-03, 'MNm':1E-06},
                'kNm':{'Nm':1E+03, 'MNm':1E-03},
                'MNm':{'Nm':1E+06, 'kNm':1E+03},
                'm':{'cm':1E+02, 'mm':1E+03},
                'cm':{'m':1E-02, 'mm':1E+01},
                'mm':{'m':1E-02, 'cm':1E-01}}

    return converter[unit_in][unit_out]


#________________________ EXCEL STUFF ________________________

def get_excel_column_names():
    xl_cols = list(string.ascii_uppercase)
    ABC = list(string.ascii_uppercase)

    for char in list(string.ascii_uppercase):
        next_level = [char + abc for abc in ABC]
        xl_cols.extend(next_level)
    return xl_cols

EXCEL_COLUMNS = get_excel_column_names()

def excel_cell_name_to_row_col(cell_name):

    row = int(cell_name[1:]) - 1 
    col = string.ascii_uppercase.index(cell_name[0])
    return row, col

def read_xl_column(xl_worksheet, start_cell:str = None, start_row:int = None, start_col:int = None, end_row:int=None, end_cell:str = None):
    '''
    reading excel colum startin from start_row, start_col ODER start_cell als Excel Zellen name
    start_row = Excel Zeilen Zahl (wird intern dann angepasst auf das 0 zählsystem von xlwings)
    start_col = A entspricht 0 
    end_row_optional: wenn nicht gegebn wird spalte bis zur nächsten leeren Zeile gelesen NOTE klappt nicht immer 
    '''
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

def read_xl_row(xl_worksheet, start_row, start_col):
    '''
    reading excel row startin from start_row, start_col
    '''
    #sheet["A1"].expand("right").last_cell.column
    end_col = xl_worksheet[start_row, start_col].expand('right').last_cell.column
    row_values = xl_worksheet[start_row, start_col:end_col].value
    return row_values

def write_to_excel(worksheet, start_cell, values, header = ['Hallo'], is_row = False, is_column = True):
    '''
    Werte in excel blatt schreibenbeginnend ab start cell
    worksheet: ein mit xlwings geöffnetes Excel Arbeitsblatt
    start_cell kann als Excel notation ('A1') oder rwo und column liste kommen 
    '''
    ws = worksheet

    for i in values:
        if isinstance(i, list):
            i = np.asarray(list)

    if isinstance(values, list):
        values = np.asarray(values)

    # convert values to a pandas dataframe -> das kann sehr gut in Excel verwendet werden
    if not isinstance(values, pd.DataFrame):
        values = values.transpose()
        df = pd.DataFrame(values, columns=header)

    if isinstance(start_cell, str):
        ws.range(start_cell).value = df
    else:
        ws.range(start_cell[0],start_cell[1]).value = df

    print ('saved', header, 'in', ws.name)

def get_simulation_history_id(case_name):

    wb = xl.Book(os_join(*['..','1_Simulations','simulation_history.xlsx']))
    ws = wb.sheets['parameter']
    header = read_xl_row(ws, 1, 2)
    col = header.index('outfile name')+2
    all_cases = ws[0:300, col].value
    row = all_cases.index(case_name)
    case_id = int(ws[row,1].value)

    return case_id

def write_stats_to_excel(excel_file, worksheet, start_cell, statistics_to_write, stats_results, case_name):
    '''
    TODO variablen die noch nicht drinnen sind im Excel händsich hinzufügen oder hier automatisch
    '''
    wb = xl.Book(excel_file)
    ws = wb.sheets[worksheet]
    header = read_xl_row(ws, 4, 2)

    if ws[start_cell[0], start_cell[1]].value != None:
        weiter = input('start_cell ist schon voll ' + worksheet + ' ' +str(start_cell[0]) + ',' +str(start_cell[1]) + ' -> wenn y wirds überschrieben? [y]/[n] ')
        if weiter == 'y':
            pass
        else:
            raise Exception('start_cell neu wählen -> row = letzte beschriebene Reihe + 1')

    ws[start_cell[0], 0].value = case_name 
    
    # id der simulation history hinzufügen
    case_id = get_simulation_history_id(case_name)
    ws[start_cell[0], 1].value = case_id

    for var in header:
        col = start_cell[1] + header.index(var)
        value = stats_results[var][statistics_to_write]
        ws[start_cell[0], col].value = value 
    wb.save()

def zellen_groeße_formatieren(excel_file, worksheet, cell_width, n_cols,  start_col = 'A', cell_height = None):
    ''' 
    excel_file: datei pfad
    Worksheet: sheet name
    cell_width: breite der zell
    n_cols: anzahl der spalten die formatiert werden sollen 
    start_col: erste Spalte ab der das formatieren los gehen soll (als Excel Buschtabe oder nummer)
    '''
    from openpyxl import load_workbook

    wb = load_workbook(excel_file)
    ws = wb[worksheet]

    if isinstance(start_col, str):
        start_col = EXCEL_COLUMNS.index(start_col)

    for i in range(start_col, n_cols):
        excel_col = EXCEL_COLUMNS[i]
        ws.column_dimensions[excel_col].width = cell_width
    wb.save(excel_file)
    wb.close()

def add_databar_color(excel_file, worksheet, columns):
    '''
    columns liste aus columns: ['B1:B12', ...]
    https://prosperocoder.com/posts/science-with-python/openpyxl-part-16-styling-with-databar/ 
    '''

    from openpyxl import load_workbook
    from openpyxl.formatting.rule import DataBarRule
    from openpyxl.styles import colors

    wb = load_workbook(excel_file)
    ws = wb[worksheet]

    ROT = '00FF0000'
    rule = DataBarRule(start_type='num', 
                   start_value=0, 
                   end_type='num', 
                   end_value=1, 
                   color=ROT) # Rot:'00FF0000'

    for col in columns:
        ws.conditional_formatting.add(col, rule)
    
    wb.save(excel_file)
    wb.close()

def zellen_farbe_formatieren(excel_file, worksheet, n_cols,  start_col = 'A', color = None):
    ''' 
    TODO funktioniert nicht
    '''
    from openpyxl import load_workbook

    wb = load_workbook(excel_file)
    ws = wb[worksheet]

    if isinstance(start_col, str):
        start_col = EXCEL_COLUMNS.index(start_col)

    for i in range(start_col, n_cols):
        excel_col = EXCEL_COLUMNS[i]
        c = ws[excel_col].value
        print()
    wb.save(excel_file)
    wb.close()

def zelle_beschriften(excel_file, worksheet, cell, wert, merge_from_to = None):
    '''
    merge_from_to: 'A2:D2'
    '''
    from openpyxl import load_workbook
    from openpyxl.styles import PatternFill

    wb = load_workbook(excel_file)
    ws = wb[worksheet]

    ws[cell].value = wert
    ws[cell].fill = PatternFill(fill_type='solid', fgColor='FFD3D3D3')
    if merge_from_to:
        ws.merge_cells(merge_from_to)
    wb.save(excel_file)
    wb.close()

def get_spalten_ausnutzung(df, df_header, start_row, start_col):

    n_rows = df.shape[0]
    end_row = start_row + n_rows + len(df_header) +1
    cols = []
    for i, val in enumerate(df_header[-1]):
        if 'Ausnutzung' in val:
            nth_col = i

    nth_col += 1 # da index mit dabei und 
    
    for i in range(df.columns.levshape[0]*df.columns.levshape[1]):
        cols.append(nth_col + i*df.columns.levshape[2])

    xl_col = [EXCEL_COLUMNS[i] for i in cols]

    cols_rows = [i + str(start_row) + ':' + i+ str(end_row) for i in xl_col]

    return cols_rows

#______________________ BERECHNUNGEN __________________________

def twr_poly_eval(x, coeffs):
    ''' beginnt erst am ^2 geht bis ^6 und c0 = 0'''
    power = 2
    y = np.zeros(x.shape)
    for c_i, c in enumerate(coeffs):
        y += c * x ** (power+c_i)

    return y

def collect_mode_shape_params(pyfast_file):

    mode_shapes_poly_coeffs = {'fore-aft':[[],[]],'side-side':[[],[]]}
    for label in pyfast_file.labels:
        if 'TwFAM' in label:
            if '1Sh' in label:
                mode_shapes_poly_coeffs['fore-aft'][0].append(pyfast_file[label])
            else:
                mode_shapes_poly_coeffs['fore-aft'][1].append(pyfast_file[label])
        if 'TwSSM' in label:
            if '1Sh' in label:
                mode_shapes_poly_coeffs['side-side'][0].append(pyfast_file[label])
            else:
                mode_shapes_poly_coeffs['side-side'][1].append(pyfast_file[label])
    return mode_shapes_poly_coeffs

def compute_sectional_mean(parameter):
        
    if not isinstance(parameter, list):
        parameter = list(parameter)
    sectional_sum = [j+i for i, j in zip(parameter[:-1], parameter[1:])]
    return np.array(sectional_sum)/2


# # DATEN MANGAMENT

def save_specific_columns_from_fast_output(cols_to_keep:list, fast_output:pd.DataFrame, fname:str):
    '''
    cols_to_keep: liste an output variablen in der FAST out datei die gespeichert werden soll
    fast_output: mit FAST_inputoutput Klasse gelesens und in dataframe konvertiertes out file von FAST
    fname: name/pfad der datei die als pkl gespeichert weren soll
    '''
    cols_to_keep.append('Time_[s]')
    fast_output.loc[:, cols_to_keep].to_pickle(fname)

    print ('\nsaved', cols_to_keep, 'in', fname)
    
# # ZÜBLIN 

def lgN_R(kfat, R, lgN_case = 1):
    '''
    lgN: erste guess dann wirds berechnet und ggf ne andere Formel genutzt
    '''
    def get_case(lgN):
        if lgN > 0 and lgN <= 5:
            return 1
        elif lgN > 5 and lgN <= 6:
            return 2
        elif lgN > 6:
            return 3

    if lgN_case == 1:
        result_lg_init = (kfat - 1) / (0.01395*R**2 + 0.004765*R - 0.06160)
    
    elif lgN_case == 2:
        result_lg_init = (kfat - (0.05494* R**2 - 0.06043*R + 1.00549)) / (0.0029*R**2 + 0.05974*R -0.0627)

    elif lgN_case == 3:
        result_lg_init = (kfat - (-0.40735* R**2 + 0.03202*R + 1.37532)) / (0.08333*R**2 + 0.04367*R -0.127)

    case = get_case(result_lg_init)
    if case == lgN_case:
        result_N = 10**result_lg_init
        return result_lg_init, result_N

    else:
        return lgN_R(kfat, R, case)


