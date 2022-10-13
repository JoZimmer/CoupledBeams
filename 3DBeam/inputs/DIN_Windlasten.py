import numpy as np 
import matplotlib.pyplot as plt

def wind_kraft(vb, category, height, cd, Aref):
    '''
    
    '''
    v_z, Iv_z, qp_z, z = DIN_potenz_profil(vb, category, height)
    we = qp_z * cd
    Fw_i = {'Fy': we * Aref}

    return Fw_i, z

def flaechen_windkraft(vb, category, height, cp_max):
    v_z, Iv_z, qp_z, z = DIN_potenz_profil(vb, category, height)
    
    F_flaeche= qp_z*cp_max
    return F_flaeche, z


def DIN_potenz_profil (vb, category, height,   at_height = 100):
    '''
    height: entweder absoluten wert geben oder diskrete Höhen Punkte als vektor oder liste
    '''

    # profile functions
    #categories = ['I','II','III','IV']
    if isinstance(height, int) or isinstance(height, float):
        z = np.arange(0,height,0.5)
    else:
        z = height

    qb = 0.5*1.25*vb**2
    
    if category == 'I':
        zmin, vmin, Imin, qbmin = 2, 0.97, 0.17, 1.9
        a, b = 1.18, 0.12
        aI, bI = 0.14, -0.12
        aq, bq = 2.6, 0.19
    elif category == 'II':
        zmin, vmin, Imin, qbmin = 4, 0.86, 0.22, 1.7
        a, b = 1.0, 0.16
        aI, bI = 0.19, -0.16
        aq, bq = 2.1, 0.24
    elif category == 'III':
        zmin, vmin, Imin, qbmin = 8, 0.73, 0.29, 1.5
        a, b = 0.77, 0.22
        aI, bI = 0.28, -0.22
        aq, bq = 1.6, 0.31
    elif category == 'IV':
        zmin, vmin, Imin, qbmin = 16, 0.64, 0.37, 1.3
        a, b = 0.56, 0.30
        aI, bI = 0.43, -0.30
        aq, bq = 1.1, 0.4
    
    #if z[1] < zmin:
    
    v_z0 = np.array([vmin*vb for i in z if i < zmin])
    Iv_z0 = np.array([Imin for i in z if i < zmin])
    qb_z0 = np.array([qbmin*qb for i in z if i < zmin])
    z1 = z[len(v_z0):]
    v_z = np.concatenate((v_z0, a*vb*(z1/10)**b), axis = 0)
    Iv_z = np.concatenate((Iv_z0, aI*(z1/10)**bI), axis = 0)
    qp_z = np.concatenate((qb_z0, qb*aq*(z1/10)**bq), axis = 0)

    return v_z, Iv_z, qp_z, z

def plot_DIN_all (vb, v_ref=1, z_ref=1, categories= ['I','II','III','IV']):

    # profile functions
    fig = plt.figure('DIN profiles', figsize=(3,3.8))
    ax = fig.add_subplot(111)
    colors = ['tab:blue','tab:orange','tab:green', 'tab:red']
    z = np.arange(0,160,0.5)
    for i, category in enumerate(categories):
        #TODO: split z array an stelle von z min und erstelle array 0 bis zmin mit vz0 wert--> concatenate them
        if category == 'I':
            zmin, vmin, Imin = 2, 0.97, 0.17
            a, b = 1.18, 0.12
            aI, bI = 0.14, -0.12
        elif category == 'II':
            zmin, vmin, Imin = 4, 0.86, 0.22
            a, b = 1.0, 0.16
            aI, bI = 0.19, -0.16
        elif category == 'III':
            zmin, vmin, Imin = 8, 0.73, 0.29
            a, b = 0.77, 0.22
            aI, bI = 0.28, -0.22
        elif category == 'IV':
            zmin, vmin, Imin = 16, 0.64, 0.37
            a, b = 0.56, 0.30
            aI, bI = 0.43, -0.30
        
        z1 = z[2*zmin+1:]
        v_z0 = np.full(2*zmin+1, vmin*vb)
        Iv_z0 = np.full(2*zmin+1, Imin)
        v_z = np.concatenate((v_z0, a*vb*(z1/10)**b), axis = 0)
        Iv_z = np.concatenate((Iv_z0, aI*(z1/10)**bI), axis = 0)


        #ax.plot(Iv_z, z/z_ref, color = colors[i])
        #ax.plot(v_z/v_ref, z/z_ref, label = 'DIN - category '+category, color = colors[i])
        ax.plot(v_z, z, label = 'DIN - category '+category, color = colors[i])

    ax.set_xlabel('v(z) [m/s]')
    ax.set_ylabel('z [m]')
    ax.legend()
    ax.grid()
    plt.show()



#plot_DIN_all(16,20,10)

