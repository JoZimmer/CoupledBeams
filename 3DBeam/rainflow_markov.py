from statistics import mean
import matplotlib.pyplot as plt 
import numpy as np
from os.path import join as os_join
from pyFAST.input_output import FASTOutputFile 
import pandas as pd
import fatpack
from source.utilities.global_definitions import VariableNames

''' 
Skript zur Bestimmung von Ermüdungswirksamen Einwirkungen anhand von Last/SGR Zeitreihen
  Begriffe Ermüdung:
    Sm = Mittelspannung
    SD = Dauerfestigkeit
    So, Su = Ober- Unterspannung
    dS/range = Schwingbreite = So - Su -> oft range gennant
    Sa = Spannungsamplitude

Daten sind bisher aus FAST Simulationen.

dataFame verwenden um mit Excel zu kommunizieren: https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_excel.html 
'''

def get_mean_range_N_dataframe(mean, range, rainflow, signal = None):
  '''
  jedem mean - range paar wird die zugehörige Schwingspielzahl zugewiesen (wenn ungleich 0) 
  result = [[id, mean, range, N], ...]
  '''
  #from utilities.Querschnitte import nEck

  mean_range_N = []
  idx = 0
  max_N = [0,0]
  for i, N_i_mean in enumerate(rainflow): # zeile = ein fester mean wert
    for j, N_ij in enumerate(N_i_mean): # spalte = ein fester range wert 
      if N_ij != 0:
        mean_range_N.append([mean[i], range[j], N_ij])
        # TODO
        # aus Funktion hier gesamte Schädigungsberechnungstabelle erstellen
        # Spannungen mean und Range berechnen -> append
        # Ober und Unterspannung, R Verhältnis berechnen -> append
        # kfat = Ober/fik berechnen, daraus Max N mit lgN_R
        # Schädigung i mit N/max_N
        if N_ij > max_N[1]:
          max_N = [idx, N_ij]
        idx += 1


  df = pd.DataFrame(mean_range_N, columns=['mean', 'range', 'N'])
  return df, max_N


# FLAGS
show_plots = True
save_plots = False

# _____________________________________________________________________________________________________________________________________________
# INPUT LESEN AUS OPENFAST OUTPUT
fast_version = '_fastv3.1.0'#'_fastv2.4.0'

case_dir = 'DLCs'#'elastoDyn_tower_params'#'base_case10_servo_off' #'windspeed10'#'base_case' #'TwFADOF1_2_True_rest_False'# 'turbsim_dlc_testing'# 

output_directory = os_join(*['..','..','OpenFAST','1_Simulations','Referenz_modelle','IEAWindTask37','IEA-3.4-130-RWT'+fast_version,'openfast', 'output', case_dir]) # ordner in dem die fast output dateien liegen
data_directory = 'fast_output_files'

discard_ramp_up =True
dt = 0.00625
ramp_up_time = 10 # [s]
discard_n_steps = int(ramp_up_time/dt)

case_file = 'dlc13_12'
file_to_analyse = os_join(*[output_directory, 'IEA-3.4-130-RWT_' + case_file + '.outb'])
#fast_output = FASTOutputFile(file_to_analyse)
#df = FASTOutputFile(file_to_analyse).toDataFrame()[discard_n_steps:]
#utils.save_specific_columns_from_fast_output(forces, df, data_directory + 'moments.pkl')

# damit nicht immer das gesamte out file gelesen werden muss
df = pd.read_pickle(data_directory + '_moments.pkl')

# _____________________________________________________________________________________________________________________________________________
# SPEZIFISCHE DATEN ZUR ANALYSE WÄHLEN ODER LOOP STARTEN
# RAINFLOW ZÄHLUNG 
#df.plot(x='Time_[s]', y=VariableNames.forces, subplots=True)
my = df[VariableNames.forces[1]].values
mx = df[VariableNames.forces[0]].values
print ('size series:', my.size)

# mit FATPACK Bibliothek
# signal
signal = np.concatenate((my,my,my,my,my,my,my,my,my,my,my,my,my,my))
print ('mean, max, min von signal', mean(signal), max(signal), min(signal))
# plt.plot(signal, label=VariableNames.forces[1])
# plt.legend()
# plt.show()
intervals = 64 #TODO was ist das genau und was beeinflusst es?
# -->  ist auf jedenfall was vom rainflow zähl algorithumus was mit dem Residuuen Problem zusammenhängt
# --> sollte auf die Gesamte Berechnung der SChädigung dann nicht soo den einfluss haben?!

# extract the ranges 
range, Sm = fatpack.find_rainflow_ranges(signal, k=64, return_means=True)
# fig, ax= plt.subplots(2)
# ax[0].plot(dS, label='range')
# ax[1].plot(Sm, label='mean')
# plt.legend()
#plt.show()

# create the bins to sort the range and mean values into.
step_size = 300
mean_bin = np.arange(Sm.min()-step_size, Sm.max()+step_size, step_size, dtype=int) # column 
range_bin = np.arange(0, range.max() + step_size, step_size, dtype=int) # row

# create array with range in the first column and mean in the second column and extract rainflow matrix
data_arr = np.array([range, Sm]).transpose()

# get mean-range matrix
rfcmat = fatpack.find_rainflow_matrix(data_arr, range_bin, mean_bin)
df, max_N = get_mean_range_N_dataframe(Sm, range, rfcmat)

# mean - Range - N als DataFrame in Excel schreibbar
# NOTE: im mode 'a' können vorhandene Tabellen ergänzt werden -> so kann ein Kopf im Excel z.B. bleiben
# with pd.ExcelWriter('Ermüdungslasten.xlsx', mode= 'a', engine="openpyxl", if_sheet_exists='overlay') as writer:
#   df.to_excel(writer, sheet_name='my', startrow=3, startcol=0)
# with pd.ExcelWriter('output.xlsx') as writer:  
#     df1.to_excel(writer, sheet_name='my')
#     df2.to_excel(writer, sheet_name='mx')


X, Y = np.meshgrid(range_bin, mean_bin, indexing='ij')
plt.figure(dpi=96)
C = plt.pcolormesh(X, Y, rfcmat, cmap='jet')
plt.colorbar(C)
plt.title("Rainflow matrix")
plt.xlim((0,6000))
plt.ylim((-8000,0))
plt.xlabel("Range")
plt.ylabel("Mean")

if show_plots:
  plt.show()

