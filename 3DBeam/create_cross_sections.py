import numpy as np
from os.path import join as os_join
import source.utilities.holz as holz 
from source.Querschnitte import KreisRing, nEck

ecken = [8,10,12]
höhen = [110,130,140,130,160]
d_oben = 3.4 # durchmesser am Kopf immer so
d_unten = 12 
n_ebenen = 14
höhen_parameter = {}
d_achse = np.linspace(12, 3.4, 14)
t_wand = 0.4
einheiten_input = {'Kraft':'N', 'Moment':'Nm', 'Festigkeit':'N/mm²', 'Länge':'m'}

lagen_aufbau = [{'ortho':'X','ti':0.08, 'di':d_achse + 0.03 + 0.04 + 0.04},
                {'ortho':'Y','ti':0.04, 'di':d_achse + 0.03 +0.02},
                {'ortho':'X','ti':0.06, 'di':d_achse},
                {'ortho':'Y','ti':0.04, 'di':d_achse -  0.03 - 0.02},
                {'ortho':'X','ti':0.08, 'di':d_achse - 0.03 -  0.04 -  0.04}]

destination = os_join(*['inputs','geometry'])

for höhe in höhen:
    höhen_parameter['absolute_höhen'] = np.linspace(0, höhe, n_ebenen)
    höhen_parameter['hfract'] = höhen_parameter['absolute_höhen']/höhen_parameter['absolute_höhen'][-1]

    kreis_ring = KreisRing(d_achse, 0.4, lagen_aufbau=lagen_aufbau,
                            holz_parameter = holz.charakteristische_werte['BSP_RFEM'], hoehen_parameter= höhen_parameter, einheiten=einheiten_input)
    ring_name = 'Ring_' + str(höhe) + '_du12.pkl'
    kreis_ring.export_to_dict_pkl(dest_file = os_join(*[destination, ring_name]))

    for n_ecken in ecken:                      
        current_nEck = nEck(n_ecken, d_achse, t_wand, lagen_aufbau=lagen_aufbau, 
                            holz_parameter = holz.charakteristische_werte['BSP_RFEM'], hoehen_parameter= höhen_parameter, einheiten=einheiten_input)
        nEck_name = str(n_ecken) + 'Eck_' + 'h'+ str(höhe) + '_du12.pkl'
        #current_nEck.plot_properties_along_height(['d_achse', 'Iy'], title=str(n_ecken) + 'Eck_' + 'h'+ str(höhe) + '_du12')
        current_nEck.export_to_dict_pkl(dest_file = os_join(*[destination, nEck_name]))