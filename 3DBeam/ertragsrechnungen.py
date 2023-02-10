import source.utilities.global_definitions as GD 
import inputs.DIN_Windlasten as windDIN
import source.utilities.utilities as utils
import numpy as np
import pandas as pd
from os.path import join as os_join
import matplotlib.pyplot as plt
from inputs.leistungskurven import anlagen, anlagen_dict
from source import writer 

def rayleigh_pdf(v_hub,v_ave):
    '''
    https://www.volker-quaschning.de/software/windertrag/index.php 
    https://www.youtube.com/watch?v=P1A-XgyP730
    '''
    scale = 2*v_ave / np.sqrt(np.pi)
    p_v = 2*v_hub / scale**2 * np.exp(-(v_hub/scale)**2)
    return p_v

def rayleigh_cdf(v_hub, v_ave):
    '''
    CDF
    mittelwert der Windgeschwindigkeitsverteilung über 10 min muss als so eine Rayleigh Verteilung angenommen werden 
    Gl. (8)
    v_ave ist das Jahres mittel der Windgeschwindigkeit auf Nabenhöhe
    v_hub sind sämtliche Windgeschwindigkeiten auf Nabenhöhe
    -> somit ist P_R die Wahrscheinlichkeit, dass v_hub gemessen wird
    '''

    P_R = 1 - np.exp(-np.pi*(v_hub/(2*v_ave))**2)

    return P_R

def leistungsbeiwert(P_el, A_rotor, v):
    '''
    v=v_hub 10 min mittel
    '''
    cp = 2*P_el /(GD.RHO_AIR * A_rotor, v)
    return cp

colors = {'WTN250':'tab:cyan','E30':'tab:green'}

show_plots = True
to_excel = False

excel = os_join(*['..','..','Kleinwindanlagen','Kleinwindanlagen.xlsx'])
ertrags_dict = {}
# NOTE WTN250 Kurve für Anlage mit 50 m Nabenhöhe

fig, (ax,ax1, axPv) = plt.subplots(3, 1, sharex=False)#, axv)
axPel = ax1.twinx()
for ai, anlage in enumerate(anlagen):
    ertrags_dict[anlage] = {'vm,hub [m/s]':[],'p(vm)':[],'T(vm) [h/a]':[], 'Pel [kW]':[]}
    axes = (ax,axPel, axPv, ax1)#, axv)
    v_hub = anlagen[anlage][0] #np.linspace(0,20, 21, True)

    v_aves = [4.5]#, 4.9]

    P_R_v0 = rayleigh_cdf(np.asarray(v_hub), v_ave=v_aves[0])

    for v_ave in v_aves:
        Pv2_v1 = []
        edges,vals = [],[]
        for v in v_hub:
            Pv1 = rayleigh_cdf(v-0.5, v_ave=v_ave)
            Pv2 = rayleigh_cdf(v+0.5, v_ave=v_ave)
            Pv2_v1.append(Pv2-Pv1)
            edges.extend([v-0.5,v+0.5])
            vals.extend([Pv2-Pv1]*2)

        if anlage == 'WTN250':
            leistung = anlagen[anlage][3] # gefittet und ausgewertetes polynom -> gibt nicht für jeden v wert in 1er schritten einen wert
        else:
            leistung = anlagen[anlage][1]  
        
        cut_in = anlagen[anlage][1].index(next(x for x in anlagen[anlage][1] if x != 0))

        leistungs_dichte = np.array(Pv2_v1[cut_in:]) * np.array(leistung[cut_in:]) # Leisutng bei bestimmter windgeschwindigkeit
        Ea = leistungs_dichte.sum() * 8760 * utils.unit_conversion('kWh', 'MWh') # kWh
        print ('Der Jahres Ertrag Ea von ', anlage, 'beträgt', round(Ea,2), 'MWh')

        axPel.plot(v_hub[cut_in:], leistung[cut_in:], label = anlage, color= colors[anlage])
        axPel.vlines(v_hub[cut_in], 0, leistung[-1], linestyles='--', color= colors[anlage], label= 'cut in')
        axPv.plot(v_hub[cut_in:], leistungs_dichte, label = 'Ea ' +str(round(Ea, 2)) + ' MWh',color= colors[anlage])

        histc = 'tab:gray'
        if ai == 0:
            ax.plot(v_hub, Pv2_v1, color = 'k', marker = 'o',markerfacecolor=histc,markeredgecolor=histc)
            ax.plot(edges, vals, label = 'pdf rayleigh für v_m ' + str(v_ave) + 'm/s', color = histc)
            #axv.plot(v_hub, np.array(Pv2_v1) * 8760, label = 'v/Jahr für v_m ' + str(v_ave) + 'm/s', color = histc)
            #axPel.bar(v_hub, np.array(Pv2_v1)/max(Pv2_v1), width=1, color = histc, label = 'Rel. Häufigkeit v_hub (max norm)')
            ax1.bar(v_hub, np.array(Pv2_v1)* 8760, width=1, color = histc, alpha = 0.5, label = 'T (v_hub [m/s])')
            #ax1.plot(v_hub, P_R_v0, color = histc, label = 'CDF(v_hub)')#, alpha = 0.5,width=1, 
        
        y = rayleigh_pdf(v_hub, v_ave)

        ertrags_dict[anlage]['vm,hub [m/s]'] = v_hub
        ertrags_dict[anlage]['p(vm)'] = Pv2_v1
        ertrags_dict[anlage]['T(vm) [h/a]'] = np.array(Pv2_v1)* 8760
        ertrags_dict[anlage]['Pel [kW]'] = leistung
        
    ax.set(ylabel='P(v)',xlim= (0,26))
    #axv.set(ylabel='T(v)/Jahr')
    axPel.set( ylabel='Pel [kW]',xlim= (0,26), ylim= (0,350))
    ax1.set( ylabel='T [h/a]', xlim= (0,26))
    axPv.set(xlabel = 'v_hub [m/s]', ylabel='P(v) x Pel(v)',xlim= (0,26))

if to_excel:
    ertrags_df = writer.df_from_nested_dict(ertrags_dict)
    with pd.ExcelWriter(excel, mode= 'w', engine="openpyxl") as xl_writer:
        ertrags_df.to_excel(xl_writer, 'Ertrag',index=False)
    utils.zellen_groeße_formatieren(excel, worksheet= 'Ertrag', cell_width=13, n_cols=len(ertrags_df))

if show_plots:
    for a in axes:
        a.set_xticks(np.arange(v_hub[-1],dtype=int),np.arange(v_hub[-1],dtype=int))
        a.grid()
    for a in (ax, axPv):
        a.legend()
    
    ax1.legend(loc=2)
    axPel.legend(loc=1)

    plt.tight_layout()
    #plt.legend()
    #plt.grid()
    plt.show()