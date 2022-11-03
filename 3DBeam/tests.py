# from math import log10
# from source.utilities import utilities as utils
# from inputs import holz
# import matplotlib.pyplot as plt
import numpy as np
# import inputs.DIN_Windlasten as windDIN

a =  [87.1, 82.7, 77.2, 71.6, 65.0, 58.3, 49.5, 49.5, 49.5, 49.5,]
x_eb = np.linspace(0,130, len(a))
x_fe = np.linspace(0,130, 20)
P_ist_fuge = np.zeros(20)

for i, x in enumerate(x_fe):
    for j, x_ebene in enumerate(x_eb[1:]):
        j += 1
        grenze_unten = x_eb[j-1]
        if x <= x_ebene and x >= grenze_unten:
            P_dif_ist = a[j]
            P_ist_fuge[i] += a[j]

    print ('an Stelle x_FE', round(x,1), 'ist P', P_dif_ist)

print(P_ist_fuge)
