from source.utilities import utilities as utils
import matplotlib.pyplot as plt
import numpy as np
import inputs.DIN_Windlasten as windDIN

v_z, Iv_z, qp_z, z = windDIN.DIN_potenz_profil(30, 'II', 160)
plt.plot(v_z,z)
plt.show()


