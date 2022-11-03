#from inputs import holz
import numpy as np
import matplotlib.pyplot as plt

class holz():
    HOLZBAU = {'k_mod':{'kurz':0.9, 'mittel':0.8, 'ständig':0.6}}

def lgN_R_zü(kfat, R, lgN_case = 1):
    '''
    Formeln von Züblin
    lgN: erst guess dann wirds berechnet und ggf ne andere Formel genutzt
    '''
    def get_case(lgN):
        if lgN > 0 and lgN <= 5:
            return 1
        elif lgN > 5 and lgN <= 6:
            return 2
        elif lgN > 6:
            return 3

    if lgN_case == 1:
        result_lg_init = (kfat - 1) / (0.01395*R**2 + 0.004765*R - 0.06160)
    
    elif lgN_case == 2:
        result_lg_init = (kfat - (0.05494* R**2 - 0.06043*R + 1.00549)) / (0.0029*R**2 + 0.05974*R -0.0627)

    elif lgN_case == 3:
        result_lg_init = (kfat - (-0.40735* R**2 + 0.03202*R + 1.37532)) / (0.08333*R**2 + 0.04367*R -0.127)

    case = get_case(result_lg_init)
    if case == lgN_case:
        result_N = 10**result_lg_init
        return result_lg_init, result_N

    else:
        return lgN_R_zü(kfat, R, case)

def kfat_din(last=None,R= 0.01,N= 1E+07):

    a = {'druck':2.0, 'zug':9.5, 'schub':6.7}
    b = {'druck':9.0, 'zug':1.1, 'schub':1.3}

    #print ('mit R=',R ,'und einer Lastspielzahl von', N , 'ergeben sich folgende kfat beiwerte')
    kfat = {}
    for last in a:
        kfat[last] = 1- ((1-R)/(a[last]*(b[last]-R))) * np.log10(N)
        #print ('    ',last, ' kfat:', round(kfat[last], 2))

    return kfat

def lgN_R_kfat_din(kfat, R, last):
    '''
    lgN als Funktion von kfat und R (Din Gleichung mit a,b einfach umgestellt)
    '''

    a = {'druck':2.0, 'zug':9.5, 'schub':6.7}
    b = {'druck':9.0, 'zug':1.1, 'schub':1.3}

    lgN = ((1-kfat) * a[last] * (b[last]-R))/(1-R)

    N = 10**lgN
    '''
    und das jetzt alles mal anschauen für
      kfat = x * kmod/gamma_m (also ansatz Schwingung um X geringer als maximale statische last)
      R = 0.1 - 0.9 (untere Grenze Validieren bzw. argumentieren)
    '''

    return lgN, N

def X_N_R(R,N,dauer):
    gamma_M = 1.25
    k_mod = holz.HOLZBAU['k_mod'][dauer]

    a = {'druck':2.0, 'zug':9.5, 'schub':6.7}
    b = {'druck':9.0, 'zug':1.1, 'schub':1.3}

    kfat = {}
    X ={}
    for last in a:
        kfat[last] = 1- ((1-R)/(a[last]*(b[last]-R))) * np.log10(N) 
        X[last] = kfat[last] * (gamma_M/k_mod)

    return X
    

def beispiel_züblin_bemessungskonzept():
    lgN = 5.5
    kfat = 0.018 # --> berechnet aus oberspannung / fik
    R = -0.25 # --> Verhältnis Ober-/Unterspannung

    res_lg, res_N = lgN_R_zü(kfat, R)

    N_auf = 9194071900
    N_Markov = 135145
    print()
    print ('log Züblin:', np.log10(N_auf))
    print ('log elf:', res_lg)
    print ('Züblin - res_N', N_auf - res_N)

    d_z = 1.47E-05
    d_i_z = N_Markov/N_auf
    d_i_s = N_Markov/res_N
    print ('\nN_Markov/N_auf', N_Markov/N_auf) 
    print ('N_Markov/result_N', N_Markov/res_N) 

    d = d_i_z - d_i_s
    rate = d/d_i_z
    print ('differenz Rate von Züblin Wert:', rate)
    print ('Faktor * selber = Züblin', d_i_s/d_i_z) # -> das eigene ist nur 74% dessen was bei Züblin drinnen steht

def plot_kfat_R(R_erwartbar = np.linspace(0,1,20), 
                N_erwartbar = np.asarray([1E+06, 1E+07, 1E+08, 1E+09])
):

    fig, ax = plt.subplots(ncols=len(N_erwartbar), num = 'kfat (N,R)')

    for i, N_i in enumerate(N_erwartbar):
        kfat_last = {'druck':[], 'zug':[], 'schub':[]}

        ax[i].set_title('N = '+ str(N_i)[0] + 'E+0' + str(int(np.log10(N_i))))
        
        for Ri in R_erwartbar:
            kfat_Ri = kfat_din(R=Ri, N=N_i)
            for last in kfat_Ri:
                kfat_last[last].append(kfat_Ri[last])

        for last in kfat_last:
            ax[i].plot(R_erwartbar, kfat_last[last], label=last)
        ax[i].grid()
        ax[i].set_xlabel('R')
    ax[0].set_ylabel('kfat')

    plt.tight_layout()
    plt.show()

def plot_X_N_R():

    N_erwartbar = np.asarray([1E+06, 1E+07, 1E+08, 1E+09])
    R_erwartbar = [-0.8]# 

    X_last = {'druck':[], 'zug':[], 'schub':[]}
    for i, N_i in enumerate(N_erwartbar):    
        for Ri in R_erwartbar:
            X = X_N_R(R=Ri, N=N_i, dauer='mittel')
            for last in X:
                X_last[last].append(X[last])
                print (N_i, last, X[last])

    fig, ax_x = plt.subplots(ncols=1, num = 'X (N,R'+str(R_erwartbar[0])+')', sharey=True)

    for last in X_last:
        ax_x.set_title('R: '+str(R_erwartbar[0]))
        ax_x.scatter(np.log10(N_erwartbar), X_last[last], label=last)
        ax_x.grid()
        ax_x.set_xlabel('lgN [-]')
        ax_x.set_ylabel('Faktor X')

    ax_x.legend()

    plt.tight_layout()
    plt.show()

#plot_X_N_R()
#plot_kfat_R()#R_erwartbar=np.linspace(-1, 0, 10))

def plot_N_ertragbar_X():
    lasten = ['druck', 'zug', 'schub']

    kmod = holz.HOLZBAU['k_mod']['mittel']
    gamma_m = 1.25 # 'neuer EC BSP

    X = np.linspace(0.5,0.9,10)

    R_erwartbar = [-0.8, -0.1, 0.1, 0.8]
    fig, ax_x = plt.subplots(ncols=len(R_erwartbar), num = 'X (N,R'+str(R_erwartbar[0])+')', sharey=False)

    for i, R_i in enumerate(R_erwartbar):
        N_ertragbar	= {'druck':[], 'zug':[], 'schub':[]}
        
        for faktor in X:
            kfat = faktor * kmod / gamma_m
            for last in lasten:
                lgN, N = lgN_R_kfat_din(kfat, R_i, last)
                N_ertragbar[last].append(N)
        for last in N_ertragbar:
            ax_x[i].set_title('R: '+str(R_erwartbar[i]))
            ax_x[i].plot(X, np.log10(np.asarray(N_ertragbar[last])), label=last)
            ax_x[i].grid()
            ax_x[i].set_xlabel('Faktor X')

    ax_x[-1].legend()
    ax_x[0].set_ylabel('lgN [-]')

    #plt.tight_layout()
    plt.show()

def plot_kfat_N():
    N_erwartbar = np.linspace(1E+05, 1E+09, 20)
    R_erwartbar = [-0.8, -0.1, 0.1, 0.8]
    fig, ax = plt.subplots(ncols=len(R_erwartbar), num = 'kfat (N,R)')

    for i, R_i in enumerate(R_erwartbar):
        kfat_last = {'druck':[], 'zug':[], 'schub':[]}

        ax[i].set_title('R = ' + str(R_i))
        
        for N_i in N_erwartbar:
            kfat_Ri = kfat_din(R=R_i, N=N_i)
            for last in kfat_Ri:
                kfat_last[last].append(kfat_Ri[last])

        for last in kfat_last:
            ax[i].plot(np.log10(N_erwartbar), kfat_last[last], label=last)
        ax[i].grid()
        ax[i].set_xlabel('lgN')

    ax[-1].legend()    
    ax[0].set_ylabel('kfat')

    plt.tight_layout()
    plt.show()

plot_kfat_N()
#plot_N_ertragbar_X()

