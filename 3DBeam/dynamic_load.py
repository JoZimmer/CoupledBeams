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
jahre = 20
ramp_up_time = 60 # in sekunden
dt = 0.01
T_simulation = 720 # s 
n_seeds = 3 # je geschwindigkeit
T_vorhanden = (T_simulation - ramp_up_time)/60 * n_seeds # minuten

# in den Ordner mit den Simulations ergebnissen gehen
zeitreihen_file = os_join(*['..','..','OpenFAST','1_Simulations','Referenz_modelle','IEAWindTask37','IEA-3.4-130-RWT_fastv3.1.0',
                        'openfast', 'output', 'DLCs', 'pickles', 'Zeitreihen','dlc12.pkl'])

with open(zeitreihen_file, 'rb') as handle:
    kopflasten_zeitreihen = pickle.load(handle)

# jetzt muss hier mit den Kopflasten eine dynamisceh windbelastung über die gesamte Höhe gemacht werden
n_nodes = 10
print(kopflasten_zeitreihen.keys())
utils.convert_coordinate_system_and_consider_einwirkungsdauer


