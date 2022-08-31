from os.path import join as os_join
#from source.utilities import utilities as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import string
from scipy.optimize import minimize_scalar, minimize
from functools import partial

def func_hi(h, hi):
     return h/hi - int(h/hi)

opt_func = partial(func_hi, 110)
b = tuple([12,15])
minimization_result = minimize_scalar(opt_func,
                                    method='Bounded',
                                    bounds=b)

h_sec = minimization_result.x
n = 110/h_sec


print(n)




