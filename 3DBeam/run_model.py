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

parameters = model_parameters.parameters['holzturm']

beam = BeamModel(parameters, coupled=False, optimize_frequencies_init=True)


print ('Frenquncies: ')
for i in range (3):
    print (round(beam.eigenfrequencies[i],4))

static_load_vector = utils.generate_nodal_force_file(beam.n_nodes, node_of_load_application=11, force_direction='y', magnitude=100000)

beam.static_analysis_solve(load_vector_file=static_load_vector)
postprocess.plot_static_result(beam, ['y'])
print ('finished')






