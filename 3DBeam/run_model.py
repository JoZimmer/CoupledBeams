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

ecken = [12] #8,10,
höhen = [110,130,140,150,160]
cross_sections = []#, 'Kreisring_d12.pkl']#'8Eck.pkl','10Eck.pkl', 
for höhe in höhen:
    #cross_sections.append('Ring_'+str(höhe)+'_du12.pkl')
    for n_ecken in ecken:
        cross_sections.append(str(n_ecken) + 'Eck_' + 'h'+ str(höhe) + '_du12.pkl')

kopf_lasten_IEA = {'Fx':3.47E+05,'Fy':5.87E+04, 'Fz':-3.51E+06, 'Mx':4.22E+06, 'My':1.28E+07, 'Mz':1.06E+06}
kopf_lasten_beam = utils.convert_coordinate_system(kopf_lasten_IEA)

for cross_section in cross_sections:
    print ('\n', cross_section)
    section_properties_pkl = os_join(*['inputs', 'geometry', cross_section])

    parameters = model_parameters.parameters['holzturm']

    parameters = utils.add_model_data_from_pkl(section_properties_pkl, parameters, set_I_eff=True)

    beam = BeamModel(parameters, adjust_mass_density_for_total = False, optimize_frequencies_init=False , apply_k_geo=False)

    # beam_tuner = Optimizations(beam)

    # target_freqs = parameters['eigen_freqs_target']

    # beam_tuner.adjust_sway_z_stiffness_for_target_eigenfreq(target_freqs[0], 
    #                                                         target_mode = 0,
    #                                                         print_to_console=True)

    # beam_tuner.adjust_mass_density_for_target_eigenfreq(target_freqs[0], 
    #                                                     target_mode = 0,
    #                                                     print_to_console=True)

    # beam_tuner.adjust_nacelle_mass_for_target_eigenfreq(target_freqs[0], 
    #                                                     target_mode = 0,
    #                                                     print_to_console=True)

    print ('  Frequencies: ')
    for i in range (3):
        print ('   ', round(beam.eigenfrequencies[i],3))

    static_load_vector = utils.generate_kopflasten_file(beam.n_nodes, kopf_lasten_beam, file_base_name='kopflasten_IEA_my_max')
    load_vector = np.load(static_load_vector)

    beam.static_analysis_solve(load_vector_file=static_load_vector)
    #postprocess.plot_static_result(beam, 'deformation', ['y'], unit='mm')
    #postprocess.plot_static_result(beam, 'internal_forces', ['y','g'], unit='kN')
    postprocess.plot_static_result_forces(beam, 'internal_forces', ['y','g'], unit='kN')
    #postprocess.plot_eigenmodes_3D(beam, beam.eigenfrequencies, beam.eigenmodes, dofs_to_plot=['y','g'])
    print ('finished')

    # VERGLEICHE
    # excel_file = os_join (*['..','..','RFEM','Ergebnisse', 'Einheitslast_ohne Vorspannung.xlsx'])
    # wb = xl.Book(excel_file)
    # ws = wb.sheets['LF6 - 4.2 Knoten - Verformungen']

    # ux_rfem = utils.read_xl_column(ws, start_cell = 'N8', end_row= 23)
    # result_rfem = {'y': np.flip((np.array(ux_rfem)))}
    # result_rfem['knoten'] = np.flip((np.array(utils.read_xl_column(ws, start_cell = 'K8', end_row= 23))))

    # postprocess.plot_static_result(beam, ['y'], unit='mm')




