import numpy as np
''' 
LF1: max Biegemoment + Querkraft und Torsionsmoment 
LF2: max Querkraft + zeitgleiches Biegemoment und Torsionsmoment
LF3: max Torsionsmoment + zeitgleiches Biegemoment und querkraft 
LF4: min Torsionsmoment + zeitgleiches Biegemoment und Querkraft 
TODO ist noch nicht vollstädnig. Oder mit panda machen
'''
IEA = {'@maxFx': {'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':6.81E+06, 'My':3.24E+06, 'Mz':2.31E+06},
       '@minFx': {'Fx':-1.02E+06,'Fy':-3.38E+04, 'Fz':-2861512.038, 'Mx':5.83E+06, 'My':5.83E+06, 'Mz':2.31E+06},
       '@maxFxy': {'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':6.81E+06, 'My':3.24E+06, 'Mz':2.31E+06},
       '@maxFxy': {'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':6.81E+06, 'My':3.24E+06, 'Mz':2.31E+06},
       '@maxFxy': {'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':6.81E+06, 'My':3.24E+06, 'Mz':2.31E+06},
       '@maxFxy': {'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':6.81E+06, 'My':3.24E+06, 'Mz':2.31E+06}
        }# NOTE Mx = Mxy #Masse charakt. 2.62E+06

import pandas as pd
from os.path import join as os_join

iea_xl = os_join(*['inputs','loads','IEAonshore.xlsx'])

iea_df = pd.read_excel(iea_xl, 'WTEnvelopes', header= [30,31], index_col=0, nrows=14)
iea_transpose = iea_df.transpose()

#print(list(iea_df.columns))
print (list(iea_transpose.columns))
print (iea_transpose['Loads @ max Fx'])


