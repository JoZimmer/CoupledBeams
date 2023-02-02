import matplotlib.pyplot as plt 
import numpy as np

#colors = plt.cm.rainbow(np.linspace(0, 1, 7))
colors = ['tab:green','tab:blue',]

Analgen = {
    #'E-70':{'Leistung':2.3,'drehzahl':[6,21]},
    #'E-82 E2':{'Leistung':2.3,'drehzahl':[6,18]}, 
    #'E-82 E3':{'Leistung':3.0,'drehzahl':[6.1,18]}, 
    #'E-82 E4':{'Leistung':3.0,'drehzahl':[6,18]}, 
    #'E-92':{'Leistung':2.35,'drehzahl':[6,16]}, 
    #'E-101':{'Leistung':3.05,'drehzahl':[6,14,5]}, 
    #'E-115 EP3':{'Leistung':4.2,'drehzahl':[4.4,12.9]}, 
    #'E-126 EP3':{'Leistung':4.0,'drehzahl':[4.4,12.4]}, 
    #'E-138 EP3':{'Leistung':4.2,'drehzahl':[4.4,11.1]}, 
    #'IEA_37':{'Leistung':3.4,'drehzahl':[3.8,12.9]}
    'WTN250':{'drehzahl':[26 , 40]}
    } 

drehzahl = np.linspace(0,20,10)
drehzahl = np.linspace(0,40,10)

# Turmeigenfrequenzen
eigenfrequenzen = {130:0.28, 140:0.25, 150:0.22, 160:0.20}
eigenfrequenzen = {130:0.28, 160:0.20}
eigenfrequenzen = {35:4.5}#, 160:0.20}
n1 = 0.26 
n1_IEA_FA = 0.305
n1_IEA_SS = 0.310


def campell_diagram(drehzahl):
    P1 = drehzahl/60
    P3 = drehzahl*3/60

    return P1, P3

plt.plot(drehzahl, campell_diagram(drehzahl)[0], label = '1P', color='k')
plt.plot(drehzahl, campell_diagram(drehzahl)[1], label = '3P', color='darkgrey')

i = 0
f_max_show = 0.8
f_max_show = 5.5
for typ, daten in Analgen.items():
    plt.vlines(daten['drehzahl'][0], ymin=0, ymax=f_max_show, label=typ, colors=colors[i])
    plt.vlines(daten['drehzahl'][1], ymin=0, ymax=f_max_show, colors=colors[i])
    i += 1

color_n = ['tab:red','tab:orange','orange','gold',]
i= 0
for h, ni in eigenfrequenzen.items():
    plt.hlines(ni, xmin=0, xmax=max(drehzahl), colors=color_n[i], label= 'h' + str(h) + ' n1='+str(ni), linestyles='--')
    i+=1

# plt.hlines(n1, xmin=0, xmax=max(drehzahl), colors='grey', label='n1='+str(n1), linestyles='--')
# plt.hlines(n1_IEA_FA, xmin=0, xmax=max(drehzahl), colors='tab:green', label='n1_IEA_FA='+str(n1_IEA_FA), linestyles='--')
# plt.hlines(n1_IEA_SS, xmin=0, xmax=max(drehzahl), colors='tab:green', label='n1_IEA_SS='+str(n1_IEA_SS), linestyles='--')
# fsize = 16
# plt.rcParams['font.size'] = str(fsize)
fsize = 14
#plt.xlim(right=16)
plt.title('Campell Diagramm', fontsize=fsize)
plt.xlabel('Drehzahl [U/min]', fontsize = fsize)
plt.ylabel('Frequenz [Hz]', fontsize = fsize)
plt.xlim(0)
plt.ylim(0)
plt.grid()
plt.legend()
plt.show()


