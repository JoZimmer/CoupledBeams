# from math import log10
import source.utilities.utilities as utils
from scipy.stats import rayleigh
#import source.utilities.fatigue_utilities as fat
# from inputs import holz
import numpy as np
import inputs.DIN_Windlasten as windDIN
# import pandas as pd
# from os.path import join as os_join
import matplotlib.pyplot as plt
# from source.writer import from_nested_dict

fig, ax = plt.subplots(1, 1)
# x = np.linspace(rayleigh.ppf(0.01),
#                 rayleigh.ppf(0.99), 100)
x = np.linspace(0,20, 100)

v_ave = 5.5
sigma = np.pi / (2*v_ave**2)
print (sigma)

# TODO das stimmt wohl och nicht so richtig 
ax.plot(x, rayleigh.pdf(x, loc = sigma, scale = 3), 'r-', lw=5, alpha=0.6, label='rayleigh pdf')
plt.legend()
plt.grid()
plt.show()