from os.path import join as os_join
from source.utilities import utilities as utils
import matplotlib.pyplot as plt
import numpy as np
h = 160
h1 = 40



k = [[2,0] , [2,3]]
u = [[0.5], [0.5]]
f = np.dot(k, u)
print (f)
