import unittest
from source.model import BeamModel
from source.Querschnitte import KreisRing
import inputs.holz as holz
import numpy as np
from source.utilities import utilities as utils
import matplotlib.pyplot as plt
from copy import copy

'''
Vergeleich der Berechnung mit einfachen händischen ergebnissen
'''
parameters_init = {
               'dimension': '2D',
                'n_elements': 3, #NOTE wird überschrieben wenn die Anzahl der Ebenen verschieden ist -> in add_model_data_from_dict
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

einheiten_input = {'Kraft':'N', 'Moment':'Nm', 'Festigkeit':'N/mm²', 'Länge':'m'}

höhen_parameter = {
                    'nabenhöhe' :130,
                    'd_unten_oben' :[12, 12],
                    'höhe_sektionen' :[12,15], # Range angeben
                    'd_knick' :None,#None, für keinen Knick -> alle folgenden Parameter sind dann unrelevant
                    'h_knick_von_oben' :70,
                    'd_unten_angepasst':12, # damit ein Übergangsknick entseht muss dieser hier kleiner als 'd_unten' sein
                    'd_knick_übergang':'automatisch',#7.5,# # bei automatisch nimmt es den schon vorhandenen an dieser Stelle
                    'n_sektionen_übergang':1,
                    'transportbreite_max':3, # past nicht so ganz zur höhe aber mus irgendwo hin
                    } 

lagen_aufbau = [{'ortho':'X','ti':0.16},
                {'ortho':'Y','ti':0.04},
                {'ortho':'X','ti':0.08},
                {'ortho':'Y','ti':0.04},
                {'ortho':'X','ti':0.16}]

holzgüte = 'C30'
n_nodes = parameters_init['n_elements']+1

werkstoff_parameter = holz.charakteristische_werte[holzgüte]
parameters_init['material_density'] = werkstoff_parameter['rhok']
parameters_init['E_Modul'] = werkstoff_parameter['E0mean']


kreis_ring = KreisRing(cd = 1.0, lagen_aufbau=lagen_aufbau,
                        holz_parameter = holz.charakteristische_werte['C30'], 
                        nachweis_parameter = holz.HOLZBAU,
                        hoehen_parameter = höhen_parameter, 
                        einheiten=einheiten_input,
                        FE_elements=parameters_init['n_elements'])

# NOTE hier das Vorzeichen ansich falsch rum aber Fehler wird nicht gefunden
n_nodes = kreis_ring.n_sections +1 #parameters_init['n_elements']+1
lasten_dict = {'Fx':np.zeros(n_nodes), 'Fy':np.zeros(n_nodes), 'Mz':np.zeros(n_nodes)}

kopflast_test = {'Fy':1E+06, 'Mz':-1E+06, 'Fx':-1E+06}
for direction in kopflast_test: # utils.update_lasten_dict in kurz
    lasten_dict[direction][-1] += kopflast_test[direction]
lasten_file_test = utils.generate_lasten_file(n_nodes, lasten_dict, 'test_kopf_1000')

section_properties = kreis_ring.section_parameters
parameters = {'FE': utils.add_model_data_from_dict(section_properties['FE'], parameters_init)}
parameters['Ebenen'] = utils.add_model_data_from_dict(section_properties['Ebenen'], parameters_init)
parameters['Ebenen']['n_elements'] = kreis_ring.n_sections

beam = BeamModel(parameters['Ebenen'], adjust_mass_density_for_total = False, 
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

    def test_querschnittswerte(self):

        print ('Iz und A mit händischem Wert Verglichen')

        lagen_aufbau = [{'ortho':'X','ti':0.16},
                        {'ortho':'Y','ti':0.04},
                        {'ortho':'X','ti':0.08},
                        {'ortho':'Y','ti':0.04},
                        {'ortho':'X','ti':0.16}]

        kreis_ring = KreisRing(cd = 1.0, lagen_aufbau=lagen_aufbau,
                        holz_parameter = holz.charakteristische_werte['C30'], 
                        nachweis_parameter = holz.HOLZBAU,
                        hoehen_parameter = höhen_parameter, 
                        einheiten=einheiten_input,
                        FE_elements=parameters_init['n_elements'])

        self.assertAlmostEqual(kreis_ring.A['Ebenen'][0], 18.10,2)
        self.assertAlmostEqual(kreis_ring.Iz['Ebenen'][0], 326.24,2)

    def test_querschnitt_nachweise(self):

        pass

        # kopflast_IEA = {'Fy':1.17E+06, 'Mz':3.24E+06, 'Nx':-3.64E+06}
        # SGR_IEA_und_wind = {'egal':{'Fy':1.17E+06, 'Mz':187E+06, 'Nx':-9.4E+06}}

        # print ('Druckspannung händisch verglichen')
        # kreis_ring.calculate_ausnutzung_normalspannung(lasten_design=SGR_IEA_und_wind)

        # self.assertAlmostEqual(kreis_ring.sigma_druck_design[0], -4.096, 2)

    def test_model(self):
        '''
        Kragarm mit kontinuierlichem Querschnitt wird die Statische Berechnung mit dem Balken FE Modell mit einer Analytischen Lösung verglichen
        Hier: maximale Durchbiegung
        '''
        t = 0.3  # t: beam thickness (y) [m]
        h = 0.3  # h: beam height (z) [m]
        E = 2.861e8  # E: Young's modulus of steel [N/mˆ2]
        nu = 3/10  # nu: Poisson's ratio

        l = 5  # l: beam element length
        A = t*h  # beam area [mˆ2]
        Iz = 1/12*t**3*h  # second moments of area [mˆ4]
        
        parameters_simple = copy(parameters_init)
        parameters_el = {'E_Modul':E,'A':A, 'Iz':Iz, 'D':0, 'nu':nu, 'lx_total_beam':l, 'n_elements': 3}
        parameters_simple.update(parameters_el)
        n_nodes = parameters_el['n_elements'] + 1

        lasten_dict_base = {'Fx':np.zeros(n_nodes), 'Fy':np.zeros(n_nodes), 'Mz':np.zeros(n_nodes)}
        Fy = 100 # N
        einzellast = {'Fy':Fy, 'Mz':0, 'Fx':0}
        for direction in einzellast: # utils.update_lasten_dict in kurz
            lasten_dict_base[direction][-1] += einzellast[direction]

        lasten_file_einzellast = utils.generate_lasten_file(n_nodes, lasten_dict_base, 'einzellast_1000')

        w_max = Fy * l**3 / 3 / (E*Iz)
        M_max = -Fy * l

        simple_beam = BeamModel(parameters_simple)
        simple_beam.static_analysis_solve(lasten_file_einzellast)
        y_max = simple_beam.static_deformation['y'][-1][0]
        g_max = simple_beam.internal_forces['g'][0]

        self.assertAlmostEqual(w_max,y_max)
        self.assertAlmostEqual(M_max,g_max)

    def test_windlasten(self):
        pass

    def test_static_analysis(self):

        beam.static_analysis_solve(lasten_file_test)

        mz_soll_base = -kreis_ring.nabenhöhe * 1E+06 - 1E+06
        qy_soll_base = 1E+06
        nx_soll_base = -1E+06

        self.assertAlmostEqual(beam.internal_forces['g'][0], mz_soll_base,2)
        # NOTE Fehler kommt von der nicht kontinuierlichen definition von EI? --> müsste mit mehr elementen besser werden?
        self.assertAlmostEqual(beam.internal_forces['y'][0], qy_soll_base,2)
        self.assertAlmostEqual(beam.internal_forces['x'][0], nx_soll_base,2)


if __name__ == '__main__':
    unittest.main()

