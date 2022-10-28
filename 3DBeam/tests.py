# from math import log10
# from source.utilities import utilities as utils
# from inputs import holz
# import matplotlib.pyplot as plt
# import numpy as np
# import inputs.DIN_Windlasten as windDIN



a= [0,0,0,0,8,6,6,5,5,4]
n_int_summe = []
n_sum_ist = 0
for i, n in enumerate(a):
    n_sum_ist += n
    n_int_summe.append(n_sum_ist)
print (n_int_summe)

