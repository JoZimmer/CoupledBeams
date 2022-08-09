import string
import numpy as np
import xlwings as xl
import pandas as pd
from os.path import join  as os_join

# # ALLGEMEINES

EXCEL_COLUMNS = list(string.ascii_uppercase)

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
                'MN':{'N':1E+06, 'kN':1E+03}}

    return converter[unit_in][unit_out]


# # EXCEL STUFF

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


# BERECHNUNGEN

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


