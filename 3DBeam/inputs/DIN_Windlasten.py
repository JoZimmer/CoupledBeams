import numpy as np 
import matplotlib.pyplot as plt


VB_WINDZONEN = {1: 22.5, 2:25, 3:27.5, 4:30} # NA.A.1
RHO = 1.225 # nach DIBt 7.3.1

def wind_kraft(vb, category, height, cd, Aref):
    '''
    
    '''
    v_z, Iv_z, qp_z, z = DIN_potenz_profil(vb, category, height)
    we = qp_z * cd
    Fw_i = {'Fx': we * Aref}

    return Fw_i, z

def flaechen_windkraft(vb, category, height, cp_max):
    v_z, Iv_z, qp_z, z = DIN_potenz_profil(vb, category, height)
    
    F_flaeche= qp_z*cp_max
    return F_flaeche, z

def DIN_potenz_profil (vb, category, height, return_at = None):
    '''
    height: entweder absoluten wert geben oder diskrete Höhen Punkte als vektor oder liste
    dann ist z in 0.5 er schritten
    '''

    # profile functions
    #categories = ['I','II','III','IV']
    if isinstance(height, int) or isinstance(height, float):
        z = np.arange(0,height,0.5)
    else:
        z = height

    qb = 0.5*RHO*vb**2
    
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
    elif category == 'dibt':
        #zmin, vmin, Imin, qbmin = 4, 0.86, 0.22, 1.7 aus category II
        zmin, vmin, Imin, qbmin = 0,0,0,0
        a, b = 1.15, 0.121
        aI, bI = 0.128, -0.05
    
    #if z[1] < zmin:
    
    vm_z0 = np.array([vmin*vb for i in z if i < zmin])
    Iv_z0 = np.array([Imin for i in z if i < zmin])
    qb_z0 = np.array([qbmin*qb for i in z if i < zmin])
    z1 = z[len(vm_z0):]
    vm_z = np.concatenate((vm_z0, a*vb*(z1/10)**b), axis = 0)
    Iv_z = np.concatenate((Iv_z0, aI*(z1/10)**bI), axis = 0)
    if category == 'dibt':
        Iv_z[Iv_z == np.inf] = 0
        qp_z = (1+7*Iv_z) * 0.5 * vm_z**2
    else:
        qp_z = np.concatenate((qb_z0, qb*aq*(z1/10)**bq), axis = 0)

    if return_at:
        v_at = vm_z[list(z).index(return_at)]
        I_at = Iv_z[list(z).index(return_at)]
        print ('Windgeschwindigkeit auf Höhe', return_at, 'm beträgt', round(v_at,2), 'm/s')
        return v_at, I_at, vm_z, Iv_z, qp_z, z
    else:
        return vm_z, Iv_z, qp_z, z

def vb_von_v_nabenhöhe (vh, category, nabenhöhe):
    '''
    Rückrechnung der basis Windgeschwidigkeit in 10 m Höhe bie gegebender Windgeschwindigkeit in Nabenhöhe
    '''
    
    if category == 'I':
        a, b = 1.18, 0.12
    elif category == 'II':
        a, b = 1.0, 0.16
    elif category == 'III':
        a, b = 0.77, 0.22
    elif category == 'IV':
        a, b = 0.56, 0.30
    elif category == 'dibt':
        a, b = 1.15, 0.121
    
    vb = vh / (a * (nabenhöhe/10)**b)

    return vb

def plot_DIN_profiles (vb, h_max ,v_ref=1, z_ref=1, categories= ['I','II','III','IV','dibt'], values = ['vm(z)','Iv(z)','qp(z)'], hline = 0):

    fig, ax = plt.subplots(ncols=len(values),num='DIN profiles')#, figsize=(3,3.8)
    
    colors = ['tab:blue','tab:orange','tab:green', 'tab:red']
    untis = {'vm(z)':' [m/s]','Iv(z)':' [-]','qp(z)':' [N/m²]'}

    for i, category in enumerate(categories):
        vm_z, Iv_z, qp_z, z = DIN_potenz_profil (vb, category, h_max)
        res = {'vm(z)':vm_z,'Iv(z)':Iv_z,'qp(z)':qp_z}


        if v_ref != 1:
            norm_by = v_ref
        else:
            norm_by = 1
        for j, value in enumerate(values):
            label = 'DIN - category '+category
            linestyle = '-'
            if category == 'dibt':
                label = 'DIBt'
                linestyle = '--'
            ax[j].plot(res[value], z, label = label, color = colors[i], linestyle=linestyle)

            ax[j].set_xlabel(value + untis[value])
            
            ax[j].grid(True)
    if hline:
        val_at = DIN_potenz_profil(vb, category, h_max, hline)
        ax[0].hlines(hline, min(res['vm(z)']), max(res['vm(z)']), label='Nabenhöhe v = ' + str(round(val_at[0],2)), color = colors[i+1])
        ax[1].hlines(hline, min(res['Iv(z)']), max(res['Iv(z)']), label='Nabenhöhe Iv = ' + str(round(val_at[1],2)), color = colors[i+1])

    ax[0].legend()
    ax[1].legend()
    ax[0].set_ylabel('z [m]')
    plt.tight_layout()
    plt.show()


#plot_DIN_profiles(vb=16, h_max=130,categories=['I','II', 'dibt'])

class Querschwingung(object):
    '''
    Wirbelerreget Querschwingung nach DIN EN 1991-1-4 Anhang E
    '''
    #from source import Querschnitte 
    # def __init__(self, querschnitt:Querschnitte.Querschnitt, n_1, vm_z, eigenform):

       
    #     self.St = 0.18 # Zylinder
    #     self.nu = 15E-06 # kinematische Viskosität Luft m²/s
        
    #     self.v_crit1 = querschnitt.d_außen[-1] * n_1 * self.St # An stelle des Turms der Maximalen modalen auslenkung

    #     ist_anfällig = False
    #     if self.v_crit1 > 1.25* vm_z[-1]:
    #         ist_anfällig = True

    #     self.Re_vcrit1 = querschnitt.d_außen[-1] * self.v_crit1 / self.nu


    def get_c_lat(self):
        '''
        Nach Tabelle E.3
        clat0 abschnittsweise konstant angenommen
        (je größer clat desto größer die Kraft - also auf der sicheren Seite die Linearen Bereich in den größeren clat0 verschieben)
        '''
        if self.Re_vcrit1 < 5E+06:
            c_lat0 = 0.7
        elif self.Re_vcrit1 > 5E+06 and self.Re_vcrit1 < 5E+07:
            c_lat0 = 0.2
        elif self.Re_vcrit1 > 5E+07:
            c_lat0 = 0.3
        

if __name__ == '__main__':
    vb =  VB_WINDZONEN[2] #3.7
    cat = 'II'
    nabenhöhe = 50
    
    DIN_potenz_profil(vb, cat, 100, nabenhöhe)
    plot_DIN_profiles(vb, 80, categories=[cat], hline=nabenhöhe, values=['vm(z)','Iv(z)'])




