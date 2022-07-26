import matplotlib.pyplot as plt
from os.path import join as os_join
from source.model import BeamModel
from source.optimizations import Optimizations
#from source.postprocess import Postprocess
from source.utilities import utilities as utils
#from plot_settings import plot_settings
from source.dynamic_analysis import DynamicAnalysis
from inputs import model_parameters
import source.postprocess as postprocess
import xlwings as xl
import numpy as np

section_properties_pkl = os_join(*['input', 'nEck_geometry', 'nEck_Lagen.pkl'])

parameters = model_parameters.parameters['holzturm']

parameters = utils.load_model_data_from_pkl(section_properties_pkl, parameters)

beam = BeamModel(parameters, adjust_mass_density = True, optimize_frequencies_init=False , apply_k_geo=False)

beam_tuner = Optimizations(beam)

target_freqs = parameters['eigen_freqs_target']

# beam_tuner.adjust_sway_z_stiffness_for_target_eigenfreq(target_freqs[0], 
#                                                         target_mode = 0,
#                                                         print_to_console=True)

# beam_tuner.adjust_mass_density_for_target_eigenfreq(target_freqs[0], 
#                                                     target_mode = 0,
#                                                     print_to_console=True)

beam_tuner.adjust_nacelle_mass_for_target_eigenfreq(target_freqs[0], 
                                                    target_mode = 0,
                                                    print_to_console=True)
print ('Frenquncies: ')
for i in range (3):
    print (round(beam.eigenfrequencies[i],4))

static_load_vector = utils.generate_nodal_force_file(beam.n_nodes, node_of_load_application=11, force_direction='y', magnitude=100000, domain_size='2D')

beam.static_analysis_solve(load_vector_file=static_load_vector)
#postprocess.plot_static_result(beam, ['y'], unit='mm')
#postprocess.plot_eigenmodes_3D(beam, beam.eigenfrequencies, beam.eigenmodes, dofs_to_plot=['y','g'])
print ('finished')

# VERGLEICHE
excel_file = os_join (*['..','..','RFEM','Ergebnisse', 'Einheitslast_ohne Vorspannung.xlsx'])
wb = xl.Book(excel_file)
ws = wb.sheets['LF6 - 4.2 Knoten - Verformungen']

ux_rfem = utils.read_xl_column(ws, start_cell = 'N8', end_row= 23)
result_rfem = {'y': np.flip((np.array(ux_rfem)))}
result_rfem['knoten'] = np.flip((np.array(utils.read_xl_column(ws, start_cell = 'K8', end_row= 23))))

postprocess.plot_static_result(beam, ['y'], unit='mm', rfem_result=result_rfem)




