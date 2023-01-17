# from math import log10
import source.utilities.utilities as utils
import source.utilities.fatigue_utilities as fat
# from inputs import holz
import numpy as np
# import inputs.DIN_Windlasten as windDIN
#import pandas as pd
from os.path import join as os_join
import matplotlib.pyplot as plt

bin_size = 0.1

fatigue_dict = {'R' : np.array(   [0.39, 0.14, -0.29, -0.304, 1.0, 1.0, 0.33, 1.0, 0.6, 0.33, 0.0, 1.0, 0.71, 0.5]),
                'N_ist': np.array([1,     1,     1,    1,   1,   1,    1,   1,   1,   1,    1,   1,   1,    1])}

b_e =       [-0.3, -0.2, -0.1, 0.0,  0.1,  0.2,  0,3,  0.4,  0.5,  0.6,   0.7,   0.8,  0.9,  1.0]
bin_numer = [     0,    1,    2,   3,    4,    5,    6,    7,     8,    9,    10,    11,   12]
anzahlen =  [     1,    0,    1,   0,    1,    0,    3,    0,     1,    0,    1,     0,   5]

min_R = min(fatigue_dict['R'])
max_R = max(fatigue_dict['R'])
n_bins = round((max_R-min_R)/bin_size + 0.5) # aufrunden
bin_edges = np.around(np.linspace(min_R, max_R, n_bins+1),1)
bin_centers = np.around((bin_edges[:-1] + bin_edges[1:]) / 2,2)

hist, bin_edges = np.histogram(fatigue_dict['R'], bins=n_bins)
bin_edges[0] -= 0.0001


N_R = np.zeros(len(bin_edges)-1)

for i, R_i in enumerate(fatigue_dict['R']):
    b = np.where(bin_edges < R_i)
    bin_ist = int(b[0][-1])
    N_R[bin_ist] += fatigue_dict['N_ist'][i]
