import numpy as np
from os.path import join as os_join
import pandas as pd
import pickle
import matplotlib.pyplot as plt

import source.utilities.global_definitions as GD
import source.utilities.utilities as utils


lokal = os_join(*['C:\\','Users','jz','Documents','WEA_lokal'])
server = 'output'
base_dir = server

# Allgemeine Informationen
ramp_up_time = 60 # in sekunden
dt = 0.01
ramp_up = int(ramp_up_time / dt)

# in den Ordner mit den Simulations ergebnissen gehen
model_directory = os_join(*['..','..','OpenFAST','1_Simulations','Referenz_modelle','IEAWindTask37','IEA-3.4-130-RWT_fastv3.1.0','openfast'])
zeitreihen_file = os_join(*['..','..','OpenFAST','1_Simulations','Referenz_modelle','IEAWindTask37','IEA-3.4-130-RWT_fastv3.1.0',
                        'openfast', 'output', 'DLCs', 'pickles', 'Zeitreihen','dlc12.pkl'])
case_dir = os_join(*['230118_holz_dlc13','output','DLCs', 'pickles' ,'dlc13.pkl'])
zeitreihen_file = os_join(*[model_directory, case_dir])
wind_case = 9
skalierung_IEA_kleinwind = 707 / 13273 *2# Rotorflächen x2


with open(zeitreihen_file, 'rb') as handle:
    kopflasten_zeitreihen = pickle.load(handle)

# jetzt muss hier mit den Kopflasten eine dynamisceh windbelastung über die gesamte Höhe gemacht werden
n_nodes = 11
h = np.array([0.0, 3.5, 7.0, 10.5, 14.0, 17.5, 21.0, 24.5, 28.0, 31.5, 35.0])
d_außen = np.array([3.0, 2.88, 2.7, 2.6, 2.52, 2.40, 2.2, 2.1, 2.03, 1.91, 1.8])

force_data, units = utils.parse_fast_labels(kopflasten_zeitreihen)


force_data = {k:force_data[k][wind_case]['time_series'][ramp_up:]*skalierung_IEA_kleinwind for k in force_data}
time_steps = force_data[list(force_data.keys())[0]].shape[0]



#fast = np.load(os_join(*['inputs','loads','dynamic','dlc13_9ms.npy']))
#beam = np.load(os_join(*['inputs','loads','dynamic','dynamic_force_11_nodes.npy']))

#time_steps = 1000

import inputs.DIN_Windlasten as wind_DIN

terrain_kategorie = 'II'
windzone = 2
# 25 m/s = cutout windspeed
basis_windgeschwindigkeit = wind_DIN.vb_von_v_nabenhöhe(wind_case, terrain_kategorie, 35) #17
vm, Iv, qp, z = wind_DIN.DIN_potenz_profil(basis_windgeschwindigkeit, terrain_kategorie, h)
#Fw_i, z = wind_DIN.wind_kraft(basis_windgeschwindigkeit)
#t_w = np.random.normal(loc=vm, scale=Iv*vm, size=(n_nodes,time_steps))
#fig, ax = plt.subplots(nrows=n_nodes)
fw_t_z = {'Fx':np.zeros((n_nodes, time_steps))}
for node in range(n_nodes):
    turb_wind = np.random.normal(loc=vm[node], scale=Iv[node]*vm[node], size=time_steps)
    fw_t_z['Fx'][node] = wind_DIN.winddruck(turb_wind, cd=1, Aref=d_außen[node])

fw_t_z['Fx'][-1] += force_data['Fx']

utils.generate_lasten_file(n_nodes, fw_t_z, 'dlc13_9ms', dimension='3D', lasten_typ='dynamic')
m = [np.mean(fw_t_z['Fx'][i]) for i in range(n_nodes)]
plt.plot(m, h)
plt.show()
# for node in range(n_nodes):
#     ax[node].plot(fw_t_z[node], label= 'mean' + str(round(np.mean(fw_t_z[node]),1)))
#     ax[node].legend()
#     ax[node].set_xlabel('time')
#     ax[node].set_ylabel('Force')
#     ax[node].hlines(np.mean(fw_t_z[node]),0,time_steps)

#plt.plot()
print('fertig')


