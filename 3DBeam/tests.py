# from math import log10
import source.utilities.utilities as utils
#from scipy.stats import rayleigh
#import source.utilities.fatigue_utilities as fat
# from inputs import holz
import numpy as np
import inputs.DIN_Windlasten as windDIN
import inputs.wind_modelle_IEC as windIEC
# import pandas as pd
# from os.path import join as os_join
import matplotlib.pyplot as plt

import numpy as np

def nonlinear_stiffness_matrix_3d(E, Ix, Iy, Iz, L, theta_x = 0, theta_y = 0, theta_z=0):
    # Berechne die Jacobi-Matrix
    dN_dX = np.array([
        [-1.0,  1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -1.0,  1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -1.0,  1.0],
        [-L/2, L/2, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, -L/2, L/2, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -L/2, L/2]
    ])
    dN_dX_inv = np.linalg.inv(dN_dX)
    dN_dU = np.dot(dN_dX_inv, np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, Ix, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, Iy, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, Iz]
    ]))
    dN_dU *= E
    
    # Berechne die Steifigkeitsmatrix
    K = np.dot(dN_dU.T, dN_dU)
    return K

#print (nonlinear_stiffness_matrix_3d(2000, 10, 200, 200,5))