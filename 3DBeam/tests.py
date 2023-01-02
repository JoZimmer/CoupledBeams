# from math import log10
# from source.utilities import utilities as utils
# from inputs import holz
# import matplotlib.pyplot as plt
import numpy as np
# import inputs.DIN_Windlasten as windDIN
import pandas as pd
from os.path import join as os_join

a = [{'ID': [0, 1, 2, 3]}, {'DLC': [5.1, 5.1, 5.1, 5.1]}]
b = {'AeroDyn':{}}

for i in a:
    k, v = list(i.items())[0]
    b['AeroDyn'][k] = v
print()