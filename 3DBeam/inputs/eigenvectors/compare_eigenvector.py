import sys
from matplotlib.pylab import *
import numpy as np

# pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
# matplotlib.rcParams.update(pgf_with_rc_fonts)

# measurements from padova
z_measured = [0, 0.115, 0.23, 0.345, 0.460]
z_adapted = [0 * 180/.46, 0.115 * 180/.46, 0.23 * 180/.46, 0.345 * 180/.46, 0.460 * 180/.46]
phi_1_measured = [0, 0.228, 0.498, 0.782, 1.0]      # weak axis
phi_2_measured = [0, 0.22, 0.535, 0.732, 1.0]       # strong axis

#results from eigenvalue analysis
eigenvector_matrix = 'inputs\\EigenvectorMatrix_mod_jz.dat'#'EigenvectorMatrix_conditioned.dat'
eigenvalues = 'inputs\\Eigenvalues.dat'#'Eigenvalues_conditioned.dat'
z = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = '),(', usecols = 0)

mode1_raw = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = '),(', usecols = 1, dtype=str)
mode2_raw = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = '),(', usecols = 2, dtype=str)
mode3_raw = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = '),(', usecols = 3, dtype=str)

modi_raw = [mode1_raw, mode2_raw, mode3_raw]


#mode1, mode2, mode3 = np.zeros(len(z),6), np.zeros(len(z),6), np.zeros(len(z),6)
z = np.insert(z, 0, 0.0)
modi = [np.zeros((len(z),6)),np.zeros((len(z),6)),np.zeros((len(z),6))]

for i in range(3):
    for z_i in range(len(z)):
        if z_i == 0:
            continue
        cur = modi_raw[i][z_i-1].split(',')
        modi[i][z_i] = np.asarray([float(val) for val in cur])



np.save('inputs\\z_coords_gid_45.npy',z)        

np.save('inputs\\EigenvectorsGid.npy', np.asarray(modi))
dof_direction_map = ['rotX', 'rotY','rotZ', 'x', 'y','z']
for i in range(3):
    fig, ax = plt.subplots(ncols=6, num='modes')
    plt.title('mode '+str(i+1))
    for dof in range(6):
        dof_z = modi[i][:,dof]
        ax[dof].plot(dof_z, z, label = 'dof ' + dof_direction_map[dof])
        ax[dof].grid()
        ax[dof].legend()

    plt.show()


# phi_1 = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = ',', usecols = 5)
# phi_2 = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = ',', usecols = 10)
# phi_3 = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = ',', usecols = 15)
# phi_4 = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = ',', usecols = 22)
# phi_5 = np.loadtxt(eigenvector_matrix, skiprows = 1, delimiter = ',', usecols = 29)

# freq_1 = round_((np.sqrt(np.loadtxt(eigenvalues, delimiter = ',', usecols = 0)) / (2*math.pi)), decimals = 3)        # doesn't work with default Eigenvalues.dat file
# freq_2 = round_((np.sqrt(np.loadtxt(eigenvalues, delimiter = ',', usecols = 1)) / (2*math.pi)), decimals = 3)
# freq_3 = round_((np.sqrt(np.loadtxt(eigenvalues, delimiter = ',', usecols = 2)) / (2*math.pi)), decimals = 3)
# freq_4 = round_((np.sqrt(np.loadtxt(eigenvalues, delimiter = ',', usecols = 3)) / (2*math.pi)), decimals = 3)
# freq_5 = round_((np.sqrt(np.loadtxt(eigenvalues, delimiter = ',', usecols = 4)) / (2*math.pi)), decimals = 3)

# for i in range(0,46):
#     print(phi_1[i] / phi_1[-1])
# for i in range(0,46):
#     print(phi_2[i] / phi_2[-1])
# for i in range(0,46):
#     print(phi_3[i] / phi_3[-1])
# for i in range(0,46):
#     print(phi_4[i] / phi_4[-1])    
# for i in range(0,46):
#     print(phi_5[i] / phi_5[-1])

#plot
# fig = plt.figure('eigenvector generic highrise', figsize=(5.85,3.5), frameon=True)
# plt.subplots_adjust(wspace=0.5)
# plt.subplot(1,8,(1,3))
# plt.plot(phi_1 / phi_1[-1], z, 'k', linewidth=1, label=r'$\Phi_1$')
# plt.plot(phi_4 / phi_4[-1], z, '--k', linewidth=1, label=r'$\Phi_4$')
# plt.plot(phi_1_measured, z_adapted, ':k', linewidth=1, label=r'$\Phi_{ref}$')
# plt.title('y-sway')
# plt.xlim(-1.1, 1.1)
# plt.grid(True)
# plt.yticks(np.arange(0, 200, 20))
# plt.xticks(ticks = [-1, -0.5, 0, 0.5, 1])
# plt.legend(loc="upper left")

# plt.subplot(1,8,(4,6))
# plt.plot(phi_2 / phi_2[-1], z, 'k', linewidth=1, label=r'$\Phi_2$')
# plt.plot(phi_5 / phi_5[-1], z, '--k', linewidth=1, label=r'$\Phi_5$')
# plt.plot(phi_2_measured, z_adapted, ':k', linewidth=1, label=r'$\Phi_{ref}$')
# plt.title('x-sway')
# plt.grid(True)
# plt.xticks(ticks = [-1, -0.5, 0, 0.5, 1])
# plt.yticks(np.arange(0, 200, 20))
# plt.xlim(-1.1, 1.1)
# plt.gca().axes.get_yaxis().set_ticklabels([])
# plt.legend(loc="upper left")

# plt.subplot(1,8,(7,8))
# plt.plot(phi_3 / phi_3[-1], z, 'k', linewidth=1, label=r'$\Phi_3$')
# plt.plot(np.array([0,0.5,1]), np.array([0,90,180]), '-.k', linewidth=1, label=r'$\Phi_{lin}$')
# plt.title('torsion')
# plt.grid(True)
# plt.xlim(-0.1, 1.1)
# plt.xticks(ticks = [0, 0.5, 1])
# plt.yticks(np.arange(0, 200, 20))
# plt.gca().axes.get_yaxis().set_ticklabels([])
# plt.legend(loc="upper left")
# fig.text(0.5, 0.01, r'normalized eigenform $\Phi_{normalized}$', ha='center')
# fig.text(0.04, 0.5, r'height $z$ [m]', va='center', rotation='vertical')
# plt.show()

# fig.savefig('/home/koenig/Desktop/Graphs/eigenvector_generic_highrise.pdf')
# fig.savefig('/home/koenig/Desktop/Graphs/eigenvector_generic_highrise.pgf')
# fig.savefig('/home/koenig/Desktop/Graphs/eigenvector_generic_highrise.svg')