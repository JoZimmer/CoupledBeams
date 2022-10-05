import unittest
from source.model import BeamModel
from source.Querschnitte import KreisRing
import inputs.holz as holz
import numpy as np
from source.utilities import utilities as utils
import matplotlib.pyplot as plt

'''
Testen des Codes. 
Am besten so wie das mit python gedacht ist. 
'''

holzgüte = 'C24'

werkstoff_parameter = holz.charakteristische_werte[holzgüte]

einheiten_input = {'Kraft':'N', 'Moment':'Nm', 'Festigkeit':'N/mm²', 'Länge':'m'}

höhen_parameter = {
                    'nabenhöhe' :110,
                    'h_knick_von_oben' :70,
                    'höhe_sektionen' :[12,15], # Range angeben
                    'd_knick' :5,#None, für keinen Knick
                    'd_unten_oben' :[12, 3.4],
                    'd_unten_angepasst':12, # damit ein Übergangsknick entseht muss dieser hier kleiner als 'd_unten' sein
                    'd_knick_übergang':'automatisch',#7.5,# # bei automatisch nimmt es den schon vorhandenen an dieser Stelle
                    'n_sektionen_übergang':1,
                    }

lagen_aufbau = [{'ortho':'X','ti':0.16},
                        {'ortho':'Y','ti':0.04},
                        {'ortho':'X','ti':0.08},
                        {'ortho':'Y','ti':0.04},
                        {'ortho':'X','ti':0.16}]


kreis_ring = KreisRing(cd = 1.1, lagen_aufbau=lagen_aufbau,
                        holz_parameter = holz.charakteristische_werte['C24'], 
                        nachweis_parameter = holz.HOLZBAU,
                        hoehen_parameter= höhen_parameter, 
                        einheiten=einheiten_input)


parameters = {
               'dimension': '2D',
                'n_elements': 20, #NOTE wird überschrieben wenn die Anzahl der Ebenen verschieden ist -> in add_model_data_from_dict
                'lx_total_beam': 110,
                'material_density': 904,# 500,#42, #für dynamische Berechnung äquivalenter Wert
                'total_mass_tower': 668390,# aus RFEM
                'nacelle_mass': 267910 ,# RFEM: 267910,# optimized: 287920.5, # Gondel Masse in kg
                'vorspannkraft':12*2.5E+06, # N
                'imperfektion':0.005, # m/m 
                'E_Modul': 12000E+06,# N/m²
                'nu': 0.1, # querdehnung
                'damping_coeff': 0.025,
                'nodes_per_elem': 2,
                'Iz': 51.0,# Not used if defined on intervals
                'dofs_of_bc':[0,1,2], # Einspannung
                'type_of_bc':'clamped',#'clamped',# or 'spring'
                'spring_stiffness':[1E+13,2E+13], # Federsteifigkeit am Boden in u und gamma richtung Bögl 40 GNm/rad
                'dynamic_load_file': '',
                'eigen_freqs_target':[0.133,0.79,3.0], 
                'defined_on_intervals':[] # kann gefüllt werden mit einer pickle datei 
            }

holzgüte = 'C24'

werkstoff_parameter = holz.charakteristische_werte[holzgüte]
parameters['material_density'] = werkstoff_parameter['rhok']
parameters['E_Modul'] = werkstoff_parameter['E0mean'] 

# NOTE hier das Vorzeichen ansich falsch rum aber Fehler wird nicht gefunden
lasten_dict_base = {'Fx':np.zeros(parameters['n_elements']+1), 'Fy':np.zeros(parameters['n_elements']+1), 'Mz':np.zeros(parameters['n_elements']+1)}

kopflast_test = {'Fy':1000, 'Mz':-1000, 'Fx':-1000}
lasten_dict_test = utils.update_lasten_dict(lasten_dict_base, wind_kraft = None, kopflasten = kopflast_test)
lasten_file_test = utils.generate_lasten_file(parameters['n_elements']+1, lasten_dict_test, 'test_kopf_1000')

section_properties = kreis_ring.section_parameters
parameters = utils.add_model_data_from_dict(section_properties, parameters)

beam = BeamModel(parameters, adjust_mass_density_for_total = False, 
                optimize_frequencies_init=False , apply_k_geo=False)

def test_element_section_number(beam_model, kreis):
    element_values = {'Iz':[], 'A':[]}
    querschnitt_values = {'Iz':kreis.Iz, 'A':kreis.A}

    beam_x = []
    kreis_ring_x = kreis.section_absolute_heights

    for elem in beam_model.elements:
        element_values['Iz'].append(elem.Iz)
        element_values['A'].append(elem.A)
        beam_x.append(elem.x_mid)

    fig, ax = plt.subplots(ncols=2)

    for i, key in enumerate(element_values):
        ax[i].plot(beam_model.nodal_coordinates['y0'],
            beam_model.nodal_coordinates['x0'],
            label = 'structure',
            marker = 'o',
            color = 'grey',
            linestyle = '--')
        ax[i].plot(element_values[key], beam_x, label='beam', marker = 'o')#, linestyle = '--')
        ax[i].plot(querschnitt_values[key], kreis_ring_x, label='QS', marker='v')
        ax[i].legend()
        ax[i].grid()
        ax[i].set_xlabel(key)
        ax[i].set_ylabel(r'height [$m$]')

    plt.show()

class TestModel(unittest.TestCase):

    def test_querschnitt(self):

        print ('Iz und A mit händischem Wert Verglichen')

        self.assertAlmostEqual(kreis_ring.A[0], 18.10,2)
        self.assertAlmostEqual(kreis_ring.Iz[0], 326.24,2)

        kopflast_IEA = {'Fy':1.17E+06, 'Mz':3.24E+06, 'Nx':-3.64E+06}
        SGR_IEA_und_wind = {'egal':{'Fy':1.17E+06, 'Mz':187E+06, 'Nx':-9.4E+06}}

        print ('Druckspannung händisch verglichen')
        kreis_ring.calculate_ausnutzung_normalspannung(lasten_design=SGR_IEA_und_wind)

        self.assertAlmostEqual(kreis_ring.sigma_druck_design[0], -4.096, 2)

    def test_windlasten(self):
        pass

    def test_static_analysis(self):

        beam.static_analysis_solve(lasten_file_test, add_eigengewicht=False, add_imperfektion=False)

        mz_soll_base = -kreis_ring.nabenhöhe * 1000 - 1000
        qy_soll_base = 1000
        nx_soll_base = -1000

        self.assertAlmostEqual(beam.internal_forces['g'][0], mz_soll_base,2)
        self.assertAlmostEqual(beam.internal_forces['y'][0], qy_soll_base,2)
        self.assertAlmostEqual(beam.internal_forces['x'][0], nx_soll_base,2)

#Summe Knotenwindkraft bei n_nodes: Ring gerade 14 nodes 648130.1479304307


if __name__ == '__main__':
    unittest.main()
    #test_element_section_number(beam, kreis_ring)
