import numpy as np
from os.path import join as os_join
import inputs.holz as holz 
from source.Querschnitte import KreisRing, nEck

'''
TODO hier muss die Eingabe noch auf die neue Version der Höhenparameter angepasst werden
'''

ecken = [8,10,12]
höhen = [110,130,140,150,160]
cd_ecken = {8:1.7,10:1.5,12:1.3}
d_oben = 3.4 # durchmesser am Kopf immer so
d_unten = 12 
n_ebenen = 14
section_höhe = 12 
höhen_parameter = {}
d_achse = np.linspace(d_unten, d_oben, n_ebenen)
t_wand = 0.4
einheiten_input = {'Kraft':'N', 'Moment':'Nm', 'Festigkeit':'N/mm²', 'Länge':'m'}
holzgüte = 'C30'

t_längslagen = [[0.16, 0.08],[0.16,0.16],[0.2,0.16]]
t_querlagen = 0.04

destination = os_join(*['inputs','geometry'])

for höhe in höhen:
    for tX in t_längslagen:
        lagen_aufbau = [{'ortho':'X','ti':tX[0]}, # 0.16, 0.2
                        {'ortho':'Y','ti':t_querlagen},
                        {'ortho':'X','ti':tX[1]}, # 0.16
                        {'ortho':'Y','ti':t_querlagen},
                        {'ortho':'X','ti':tX[0]}] # 0.16, 0.2
        höhen_parameter['absolute_höhen'] = np.linspace(0, höhe, n_ebenen)

        kreis_ring = KreisRing(d_achse, cd=1.1, lagen_aufbau=lagen_aufbau,
                                holz_parameter = holz.charakteristische_werte[holzgüte], 
                                nachweis_parameter = holz.HOLZBAU,
                                hoehen_parameter= höhen_parameter, einheiten=einheiten_input)
        ring_name = 'Ring_' + str(höhe) + '_du' +str(int(d_achse[0])) + '_' + holzgüte + '_tX' + str(int(kreis_ring.t_laengslagen*100)) + '_tY' + str(int(kreis_ring.t_querlagen*100)) +'.pkl'
        kreis_ring.export_object_to_pkl(dest_file = os_join(*[destination, 'objekte', ring_name]))

        for n_ecken in ecken:                      
            current_nEck = nEck(n_ecken, d_achse, cd=cd_ecken[n_ecken], lagen_aufbau=lagen_aufbau, 
                                holz_parameter = holz.charakteristische_werte[holzgüte], 
                                nachweis_parameter = holz.HOLZBAU,
                                hoehen_parameter= höhen_parameter, einheiten=einheiten_input)

            nEck_name = str(n_ecken) + 'Eck_' + 'h'+ str(höhe) + '_du' +str(int(d_achse[0])) + '_' + holzgüte + '_tX' + str(int(current_nEck.t_laengslagen*100)) + '_tY' + str(int(current_nEck.t_querlagen*100)) +'.pkl'
            #current_nEck.plot_properties_along_height(['d_achse', 'Iy'], title=str(n_ecken) + 'Eck_' + 'h'+ str(höhe) + '_du12')
            current_nEck.export_object_to_pkl(dest_file = os_join(*[destination, 'objekte', nEck_name]))