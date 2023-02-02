from os.path import join as os_join
import matplotlib.pyplot as plt
import numpy as np


source = os_join(*['inputs', 'airfoil_data'])
source = os_join(*['..', '..', 'OpenFAST','1_Simulations','Referenz_modelle','IEAWindTask37','IEA-3.4-130-RWT_fastv3.1.0','openfast','Airfoils'])

datei = 'AH 93-W-145.dat' #['MH102.dat', 'MH104.dat','MH106.dat', 'MH108.dat', 'MH110.dat']

for num in range(30):

    if num < 10:
        datei = 'AF0' +str(num) + '_Coords.txt'
    else:
        datei = 'AF' +str(num) + '_Coords.txt'

    airfoildata = np.asarray([i.strip().split() for i in open(os_join(source,datei)).readlines()[9:]], dtype=float)

    x = airfoildata[:,0]

    plt.plot(airfoildata[:,0], airfoildata[:,1], label = datei)

plt.legend()  
plt.ylim((-0.5,0.5))
plt.show()



