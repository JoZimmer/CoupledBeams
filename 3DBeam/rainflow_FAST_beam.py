import numpy as np
from os.path import join as os_join
import pandas as pd
import pickle
import fatpack
import matplotlib.pyplot as plt

import source.postprocess as postprocess
import source.utilities.utilities as utils
import source.utilities.fatigue_utilities as fat
from source.Querschnitte import KreisRing
import source.utilities.global_definitions as GD
from inputs import holz

'''
Input: 
    - Aus x (meist 3) seeds je Windgeschwindigkeit kombinierte Zeitreihen von DLC 1.2 Simuliert mit OpenFAST als pickle / dictionary
    - Definierter Kreisquerschnitt eines Turms

Berechnung:
    - Mit Einflusslinien werden die Schnittgrößen im Turm anhand der Kopflasten Zeitreihen durch eine zustäzliche, konstante Windbelastung
      des Turms and n (hier den defierten Ebenen) bestimmt
      TODO -> die konstante Belastung wirkt sich nicht unbeding ungünstig auf die Ermüdungs berechnung aus (Höherer Mittelwert bei gleich bleibenden ranges)
    - Die Windlast wird auf Basis der Input Windgeschwindigkeit der entsprechenden Simulation nach DIN berechnet  

    - Querkraft und  Torsion werden schon auf Zeitreihen ebene in eine Schubspannungszeitreihe umgerechnet 
      NOTE Kraft kommt von OpenFAST in [kN] Fläche im beam in [m²] -> Spannungen werden in [kN/m²] gespeichert
        - ggf. getrennt nach Lagen wenn mit furnier ebenen

    - Rainflow Zählung jeder Zeitreihe: Jeder Schnittgrößen Komponete an jedem Knoten
    - Extrapolation der 'N'-Werte der resultierenden Markov Matrix anhand der Input Windgeschwindigkeit und der Rayleigh-Verteilung

Output: 
    - Markov Matritzen der relevanten SGR (Mz, Mx, Qy) für Biegung und Schub, im Fall von Schub schon Spannungen 
    - Gespeichert als dicts in pickles im Ordner output//Markov

Plots:
    - Markov Matrix als colormesh in 2D
    - Markov Matrix als 3D colorbarmesh
    - Reversal Darstellungen in den Zeitreihen

'''
lokal = os_join(*['C:\\','Users','jz','Documents','WEA_lokal'])
server = 'output'
base_dir = server

# Allgemeine Informationen
jahre = 20
ramp_up_time = 60 # in sekunden
dt = 0.01
T_simulation = 720 # s 
n_seeds = 3 # je geschwindigkeit
T_vorhanden = (T_simulation - ramp_up_time)/60 * n_seeds # minuten

# Kopflasten Zeitreihen werden für vieles gebraucht
# in den Ordner mit den Simulations ergebnissen gehen
zeitreihen_file = os_join(*['..','..','OpenFAST','1_Simulations','Referenz_modelle','IEAWindTask37','IEA-3.4-130-RWT_fastv3.1.0',
                        'openfast', 'output', 'DLCs', 'pickles', 'Zeitreihen','dlc12.pkl'])

#zeitreihen_file = os_join(*[base_dir, 'Zeitreihen','dlc12.pkl'])
with open(zeitreihen_file, 'rb') as handle:
    kopflasten_zeitreihen = pickle.load(handle)

fatigue_excel = os_join(*['output','Berechnungs_Ergebnisse_Rainflow.xlsx'])

SGR_zeitreihen_vorhanden = True
add_windkraft_konst = False
zeitreihen_file_name = 'schnittgrößen_all_furnier_ohne_wind.pkl'
save_outfile = False
mit_furnier = True
# Markov Zählung
markov = True
v_list_custom = [9,14.5,20,25] # oder False 
plot2d = False
plot3d = True
show_plots = True
save_plots = False
to_excel = False # nur die Querschnitts und Skalierungsinformationen werden geschrieben

if not SGR_zeitreihen_vorhanden:

    parameters_init = {
                'dimension': '2D',
                    'n_elements': 20, #
                    'lx_total_beam': 130,
                    'material_density': 904,# 500,#42, #für dynamische Berechnung äquivalenter Wert
                    'total_mass_tower': 668390,# aus RFEM
                    'nacelle_mass': 267910,# IEA37: Gondel Masse in kg # aus optimierung 287920.5
                    'vorspannkraft':12*2.5E+06, # N
                    'imperfektion':0.008, # m/m Schiefstellung + Fundament schief
                    'E_Modul': 12000E+06,# N/m²
                    'nu': 0.1, # querdehnung
                    'damping_coeff': 0.025,
                    'nodes_per_elem': 2,
                    'Iz': 51.0,# Not used if defined on intervals
                    'dofs_of_bc':[0,1,2], # Einspannung
                    'type_of_bc':'Feder', #'clamped',#'Eingespannt',# or 
                    'spring_stiffness':[1E+13,2E+13], # Federsteifigkeit am Boden in u und gamma richtung Bögl 40 GNm/rad
                    'dynamic_load_file': os_join(*["inputs","forces","dynamic_force_11_nodes.npy"]),
                    'eigen_freqs_target':[0.133,0.79,3.0], 
                    'defined_on_intervals':[] # kann gefüllt werden mit einer pickle datei 
                }
    holzgüte = 'C30'
    n_nodes = parameters_init['n_elements']+1

    werkstoff_parameter = holz.charakteristische_werte[holzgüte]
    parameters_init['material_density'] = werkstoff_parameter['rhok']
    parameters_init['E_Modul'] = werkstoff_parameter['E0mean']

    nachweis_parameter = holz.HOLZBAU

    '''
    QUERSCHNITTS DEFINITION
    '''

    querschnitte,qs_werte, qs_werte_FE, qs_labels  = [],[],[],[]
    höhe = 130
    dicken = [44]#[44]#[40]# ,, 44]#,48,56,64] # key für lagenaufbau dictonaire in holz 
    
    #mit_furnier = True # dicke 36 sonst 44

    werkstoff_parameter['Furnierebene'] = mit_furnier
    furnier_dict = {True:'BSP_furnier',False:'BSP_normal'}
    #knicke = {'gerade':[None,10], '1Knick':[4.5,10], 'Knick mit Übergang':[4.5,8]}
    knicke = {'Gerade':[None, 11]}
    for label, knick in knicke.items():
        for t in dicken:
                höhen_parameter = {
                                'nabenhöhe' :höhe,
                                'd_unten_oben' :[11, 3.4],
                                'höhe_sektionen' :[12,15], # Range angeben
                                'd_knick' :knick[0],#None, für keinen Knick -> alle folgenden Parameter sind dann unrelevant
                                'h_knick_von_oben' :70,
                                'd_unten_angepasst':knick[1], # damit ein Übergangsknick entseht muss dieser hier kleiner als 'd_unten' sein
                                'd_knick_übergang':'automatisch',#7.5,# # bei automatisch nimmt es den schon vorhandenen an dieser Stelle
                                'n_sektionen_übergang':1,
                                'transportbreite_max':3, # past nicht so ganz zur höhe aber mus irgendwo hin
                                } 

                lagen_aufbau = holz.lagenaufbauten[furnier_dict[mit_furnier]][t]
                t_ges = sum([lage['ti'] for lage in lagen_aufbau])

                einheiten_input = {'Kraft':'N', 'Moment':'Nm', 'Festigkeit':'N/mm²', 'Länge':'m'}
                
                cd_zylinder = 1.0
                kreis_ring = KreisRing(cd = cd_zylinder, cp_max=1.8, lagen_aufbau=lagen_aufbau,
                                        holz_parameter = werkstoff_parameter, 
                                        nachweis_parameter = nachweis_parameter,
                                        hoehen_parameter= höhen_parameter, einheiten=einheiten_input,
                                        FE_elements=parameters_init['n_elements'])

                kreis_ring.name += ' ' + label + ' t' + str(t)
                qs_label_openfast = kreis_ring.name + '_h' + str(höhen_parameter['nabenhöhe']) + '_du' + str(höhen_parameter['d_unten_oben'][0]) +\
                                '_do' + str(höhen_parameter['d_unten_oben'][1])
                qs_labels.append(kreis_ring.name)
                qs_values_digits = 2
                for values in kreis_ring.querschnitts_werte['Ebenen'].values():
                    qs_werte.append(np.around(values, qs_values_digits))
                for values in kreis_ring.querschnitts_werte['FE'].values():
                    qs_werte_FE.append(values)

                querschnitte.append(kreis_ring)
    
        qs_header = pd.MultiIndex.from_product([qs_labels, list(kreis_ring.querschnitts_werte['Ebenen'].keys())])
        qs_header_FE = pd.MultiIndex.from_product([qs_labels, list(kreis_ring.querschnitts_werte['FE'].keys())])
        #lagen_header = pd.MultiIndex.from_product([['Lagenaufbau [m]'], list(lagen_aufbau[0].keys())])

        qs_werte=np.array(qs_werte).transpose()
        qs_werte_FE=np.array(qs_werte_FE).transpose()
        qs_df = pd.DataFrame(qs_werte, columns=qs_header)
        qs_FE_df = pd.DataFrame(qs_werte_FE, columns=qs_header_FE)
        holz_df = pd.DataFrame(werkstoff_parameter, index=[0])
        lagen_df = pd.DataFrame(lagen_aufbau,  index=list(range(1,len(lagen_aufbau)+1)))
        lagen_df.rename({'ti': 'ti [' + kreis_ring.einheiten['Länge'] +']', 'di': 'di [' + kreis_ring.einheiten['Länge']+']'}, axis='columns',inplace=True)

    '''
    ÄUßERE BELASTUNG FATIGUE + mittlere Turmbelastung TODO: Alternative SGR Dynmaisch berechnen
    NOTE: hier wird jetzt mal einfach nur von einem Querschnitt ausgegangen und nicht mehreren
    '''
    # sehr händisch die länge der zeitreihen holen
    zeitreihen_len = kopflasten_zeitreihen[list(kopflasten_zeitreihen.keys())[0]][4]['time_series'].shape[0]
    v_in_simuliert = list(kopflasten_zeitreihen[list(kopflasten_zeitreihen.keys())[0]])

    import inputs.DIN_Windlasten as wind_DIN

    terrain_kategorie = 'II'
    knoten_typ = 'Ebenen'
    qs = kreis_ring
    gamma_q_fatigue = 1.0 # TODO

    # NOTE Achtung: OpenFASt Last kommt in kN hier ist alles in N
    if add_windkraft_konst:
        knoten_wind_kraft_z = {}
        for v_in in v_in_simuliert:
            vb_10m = wind_DIN.vb_von_v_nabenhöhe(v_in, terrain_kategorie, kreis_ring.nabenhöhe)
            wind_kraft_z, z_coords = wind_DIN.wind_kraft(vb=vb_10m, category=terrain_kategorie, height=qs.section_absolute_heights[knoten_typ], cd = qs.cd, Aref=qs.d_achse[knoten_typ]) # ergibt N/m
            knoten_wind_kraft_z[v_in] = utils.linien_last_to_knoten_last(wind_kraft_z, qs.section_absolute_heights[knoten_typ], gamma_q = gamma_q_fatigue)
            
    einfluss_funktion =  {}
    # Einfluss Funktion der Last in "last_richtung" auf die Schnittgröße "reaktion" am knoten "at_node"
    for reaktion in ['Mz', 'Mx', 'Qy']:
        einfluss_funktion[reaktion] = {}
        for last_richtung in ['y','g','a']:
            einfluss_funktion[reaktion][last_richtung] = []
            for at_node in range(kreis_ring.n_ebenen):
                einfluss_funktion[reaktion][last_richtung].append(np.zeros(kreis_ring.n_ebenen))
                for node in range(kreis_ring.n_ebenen):
                    einfluss_funktion[reaktion][last_richtung][at_node][node] = utils.kragarm_einflusslinien(kreis_ring.section_absolute_heights['Ebenen'],
                                                load_direction=last_richtung, node_id=node, response=reaktion, response_node_id=at_node)

    # schnittgrößen_veränderlich = {knoten_typ : {kreis_ring.name: {kreis_ring.nabenhöhe: {}}}} NOTE die unterscheidung der QS und nabenhöhen macht das ganze unnötg tief -> verwende gerade eh nur einen QS
    schnittgrößen_zeitreihen = {}
    # direkt resultierende verwenden # die LAsten von OpenFASt sind ja schon schnittgrößen
    last_openFAST = {'y':'YawBrFxyp_[kN]', 'a':'YawBrMzp_[kN-m]', 'g':'YawBrMxyp_[kN-m]'} 

    # SGR - reihen: Höhe, Spalten: Zeit
    '''
    Berechnung der Schnittgrößen simpel mit Einfluss linien
    TODO Alternativ Dynamisch
    '''
    for reaktion in ['Mz', 'Mx', 'Qy']:
        schnittgrößen_zeitreihen[reaktion] = {}

        for i, v_in in enumerate(v_in_simuliert):
            schnittgrößen_zeitreihen[reaktion][v_in] = np.zeros((kreis_ring.n_ebenen, zeitreihen_len)) #reihen: Höhe, Spalten: Zeit

            for at_node in range(kreis_ring.n_ebenen):        
                for last_richtung in ['y','g','a']:

                    einfluss= einfluss_funktion[reaktion][last_richtung][at_node] 
                    kopflast_t = kopflasten_zeitreihen[last_openFAST[last_richtung]][v_in]['time_series']
                    sgr_kopf_t = einfluss[-1] * kopflast_t# kopflast ist der oberste knoten
                    if add_windkraft_konst and last_richtung == 'y': #knoten windkraft nur in y-Richtung
                        sgr_konst = sum(einfluss * knoten_wind_kraft_z[v_in][GD.DOF_LOAD_MAP[last_richtung]] * utils.unit_conversion('N', 'kN'))
                    else:
                        sgr_konst = 0.
                    schnittgrößen_zeitreihen[reaktion][v_in][at_node] += sgr_kopf_t + sgr_konst

    # stimmt jetzt nicht so ganz aber ja
    if kreis_ring.holz_parameter['Furnierebene']:
        schnittgrößen_zeitreihen['tau_Fe'], schnittgrößen_zeitreihen['tau_längs'] = {},{}
    else:
        schnittgrößen_zeitreihen['tau'] = {}
    for i, v_in in enumerate(v_in_simuliert):
        if kreis_ring.holz_parameter['Furnierebene']:
            schnittgrößen_zeitreihen['tau_Fe'][v_in] = np.zeros((kreis_ring.n_ebenen, zeitreihen_len))
            schnittgrößen_zeitreihen['tau_längs'][v_in] = np.zeros((kreis_ring.n_ebenen, zeitreihen_len))
        else:
            schnittgrößen_zeitreihen['tau'][v_in] = np.zeros((kreis_ring.n_ebenen, zeitreihen_len))  #reihen: Höhe, Spalten: Zeit
        for at_node in range(kreis_ring.n_ebenen):  
            sgr = {'Qy':schnittgrößen_zeitreihen['Qy'][v_in][at_node], 'Mx':schnittgrößen_zeitreihen['Mx'][v_in][at_node]}

            kreis_ring.calculate_schubspannung_scheibe_at_knoten(sgr, knoten=at_node, knoten_typ='Ebenen') #NOTE das wird knoten weise gemacht da an jedem knoten zeitreihen vorliegen
            
            if kreis_ring.holz_parameter['Furnierebene']:
                schnittgrößen_zeitreihen['tau_Fe'][v_in][at_node] += kreis_ring.tau_Fe_design['Ebenen']
                schnittgrößen_zeitreihen['tau_längs'][v_in][at_node] += kreis_ring.tau_längs_design['Ebenen']
            else:
                schnittgrößen_zeitreihen['tau'][v_in][at_node] += kreis_ring.tau_xy_design['Ebenen'] #* utils.unit_conversion('kN/m²', 'N/mm²') # SGR kommen in kN von openFAST, QS-Werte im QS-Objekt sind alle in [m]

    output_folder = os_join(*['output', 'Zeitreihen'])
    if save_outfile:
        outfile = os_join(*[output_folder, zeitreihen_file_name])
        with open(outfile, 'wb') as handle:
            pickle.dump(schnittgrößen_zeitreihen, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print ('\nSaved', outfile)

if SGR_zeitreihen_vorhanden:
    zeitreihen_sgr_file = os_join(*[base_dir,'Zeitreihen',zeitreihen_file_name]) # server: 'output'
    with open(zeitreihen_sgr_file, 'rb') as handle:
        schnittgrößen_zeitreihen = pickle.load(handle)
    v_in_simuliert = list(kopflasten_zeitreihen[list(kopflasten_zeitreihen.keys())[0]])

'''
Rainflow Zählung
'''
output_folder = os_join(*['output', 'Markov'])

einheiten = {'Mz':'_[MNm]', 'Mx':'_[kNm]', 'Qy':'_[kN]', 'tau':'_[kN/m²]','tau_Fe':'_[kN/m²]','tau_längs':'_[kN/m²]', 'YawBrFxyp_[kN]':''}

# ----- Die eigentliche rainflow Zählung und generierung der Markov Matritzen für die weiter verwendung
# -----------------------------------------------------------------------------------------------------
if markov:
    if v_list_custom:
        v_in_simuliert = v_list_custom
    markov_gesamt, markov_sammlung = {}, {}

    step_size_reaktion = {'Mz':600, 'Mx':100, 'Qy':10, 'tau':1, 'tau_Fe':10, 'tau_längs':1, 'YawBrFxyp_[kN]':10}#* utils.unit_conversion('kN/m²', 'N/mm²')
    
    for reaktion in ['Mz']:#['tau_längs']:#['tau_Fe']:#,['YawBrFxyp_[kN]']:#['Mx', 'Qy']:# , 
        markov_sammlung[reaktion], markov_gesamt[reaktion] = {}, {}
        for at_node in [3]:# range(10):#range(kreis_ring.n_ebenen):#
            markov_sammlung[reaktion][at_node], markov_gesamt[reaktion][at_node] = {}, {}

            skalierung = {  'v_in-[m/s]':v_in_simuliert, 'ramp_up-[s]':[ramp_up_time]*len(v_in_simuliert), 'T_seed-[s]':[T_simulation]*len(v_in_simuliert),
                            'T_vorh.-[min]':[T_vorhanden]*len(v_in_simuliert), 'T_soll.-[J]':[jahre]*len(v_in_simuliert), 'P_R(v1 < v < v2)':[], 'Skalierung-[-]':[]}

            print ('At node', at_node)
            for i, v_in in enumerate(v_in_simuliert):
                
                # aktuell ausgewertete Zeitreihe
                if reaktion != 'YawBrFxyp_[kN]':
                    current_signal = schnittgrößen_zeitreihen[reaktion][v_in][at_node]
                else:
                    current_signal = kopflasten_zeitreihen[reaktion][v_in]['time_series']

                # nach dem für einen knoten die gesamte schnittgröße(t) bestimmt ist kann auch gleich die Zählung stattfinden
                intervals = 100 #TODO was ist das genau und was beeinflusst es?
                # -->  ist auf jedenfall was vom rainflow zähl algorithumus was mit dem Residuuen Problem zusammenhängt
                # --> sollte auf die Gesamte Berechnung der Schädigung dann nicht so den einfluss haben?!

                #print ('   für v =', v_in, 'maximales tau:', round(schnittgrößen_zeitreihen[reaktion][v_in][at_node].max() * utils.unit_conversion('kN/m²', 'N/mm²'),2))
                # extract the ranges 
                rainflow_range, Sm = fatpack.find_rainflow_ranges(current_signal, k=intervals, return_means=True)

                # create the bins to sort the rainflow_range and mean values into.
                step_size = step_size_reaktion[reaktion] # TODO sollte gewählt werden mit sinn oder einfach mal vershciedene testen und gesamt outcome der schädigung vergleichen

                # Für einheitlich bin und ranges der verschiedenen Zeitreihen wird der start der means am nächst kleineren vielfachen der step size ausgehend von Sm.min gewählt
                Sm_start = fat.find_mean_start(int(Sm.min()), step_size) 
                mean_bin = np.arange(int(Sm_start), int(Sm.max()+step_size), step_size, dtype=int) # column NOTE arange funktion hat ein paar eigenheiten siehe definition im internet - so passt es jetzt aber
                range_bin = np.arange(0, rainflow_range.max() + step_size, step_size, dtype=int) # row -> durch transpose wird das umgedreht für die rfcamt bestimmung

                # create array with rainflow_range in the first column and mean in the second column and extract rainflow matrix
                data_arr = np.array([rainflow_range, Sm]).transpose() # -> mean: row, range: col
                T_soll_v = kopflasten_zeitreihen['YawBrFxyp_[kN]'][v_in]['rayleigh'] * jahre * 365 * 24 * 60 # minuten #bisschen schlecht die rayleigh werte für jede kraftrichtung zu speicher
                faktor_X_Jahre = T_soll_v/T_vorhanden
                #print ('    Für v_in', v_in, '[m/s] Skalierung auf 20 Jahre mit Faktor', round(faktor_X_Jahre,2))
                skalierung['Skalierung-[-]'].append(round(faktor_X_Jahre))
                skalierung['P_R(v1 < v < v2)'].append(round(kopflasten_zeitreihen['YawBrFxyp_[kN]'][v_in]['rayleigh'],8))

                # get mean-range matrix
                rfcmat = fatpack.find_rainflow_matrix(data_arr, range_bin, mean_bin)
                rfcmat *= int(faktor_X_Jahre)

                #postprocess.plot_markov(range_bin, mean_bin, rfcmat, show_plots = True, save_plots = False, fig_name_save='', title= '\n'+ reaktion + ' at node' + str(at_node) + ' v' + str(v_in)+ 'm/s')

                # sammeln der markov daten für jede geschwindigkeit
                if i == 0:
                    markov_sammlung[reaktion][at_node] = {'mean':[mean_bin],'range':[range_bin],'N':[rfcmat]} # N hier in matrix form rfcmat
                else:
                    markov_sammlung[reaktion][at_node]['mean'].append(mean_bin)
                    markov_sammlung[reaktion][at_node]['range'].append(range_bin)
                    markov_sammlung[reaktion][at_node]['N'].append(rfcmat)

            # kombinieren der einzelenen Markovs je v_in
            markov_gesamt[reaktion][at_node] = fat.combine_markovs(markov_sammlung[reaktion][at_node], step_size)

            if reaktion == 'Mz':
                unit_factor = utils.unit_conversion('kNm','MNm')
            else:
                unit_factor = 1
            if plot2d:
                postprocess.plot_markov(markov_gesamt[reaktion][at_node]['range'], 
                                        markov_gesamt[reaktion][at_node]['mean'], 
                                        markov_gesamt[reaktion][at_node]['N'], 
                                        show_plots = show_plots, save_plots = save_plots, 
                                        fig_name_save=os_join(*['output', 'plots','Markov_Matrizen', reaktion + '_at_node' + str(at_node) + '_v_all.png']), 
                                        title= 'bin size ' + str(step_size)+ '\n'+ reaktion + einheiten[reaktion] + ' at node' + str(at_node) + ' v all')
            if plot3d:
                postprocess.plot_markov_3d( markov_gesamt[reaktion][at_node]['range'], 
                                            markov_gesamt[reaktion][at_node]['mean'], 
                                            markov_gesamt[reaktion][at_node]['N'], 
                                            bin_size= step_size,
                                            unit_factor = unit_factor,
                                            show_plots = show_plots, save_plots = save_plots, 
                                            fsize = 16,
                                            fig_name_save=os_join(*['output', 'plots','Markov_Matrizen', '3d_' + reaktion + '_at_node' + str(at_node) + '_v_all.png']), 
                                            title= 'bin size ' + str(step_size)+ '\n'+ reaktion + einheiten[reaktion] + ' at node' + str(at_node) + ' v all')

    if save_outfile:
        spannungs_art = {'Mz':'druck', 'tau':'schub', 'tau_Fe':'schub_Fe', 'tau_längs':'schub_längs'}
        outfile = os_join(*[output_folder, 'markov_matrix_'+ spannungs_art[reaktion] +'_DLC12.pkl'])
        with open(outfile, 'wb') as handle:
            pickle.dump(markov_gesamt, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print ('\nSaved', outfile)

    if to_excel:
        skalierung_df = pd.DataFrame(skalierung)

        with pd.ExcelWriter(fatigue_excel, mode= 'w', engine="openpyxl") as writer:#, if_sheet_exists='replace'
            holz_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=0, startcol=0, index=False)#
            qs_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=4, startcol=0)
            lagen_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=qs_df.shape[0] + 10, startcol=0)
            skalierung_df.to_excel(writer, sheet_name= 'Rayleigh_Skalierung', startrow=0, startcol=0, index=False)#

        utils.zellen_groeße_formatieren(fatigue_excel, worksheet= 'QS_Werte', cell_width=13, n_cols=max(len(qs_df.columns), len(holz_df.columns))+1)
        utils.zellen_groeße_formatieren(fatigue_excel, worksheet= 'Rayleigh_Skalierung', cell_width=16, n_cols=len(skalierung_df.columns)+1)

# ----- Plots und tests zur Rainlfow Zählung
# ---------------------------------------------------------------------------
reversal = False

if reversal:
    for reaktion in ['Mz']:#['tau']:#['Mx', 'Qy']:# , ['YawBrFxyp_[kN]']
        for at_node in range(1):#range(kreis_ring.n_ebenen):
            for i, v_in in enumerate(v_in_simuliert):
                #current_signal = kopflasten_zeitreihen[reaktion][v_in]['time_series'][:1000]
                current_signal = schnittgrößen_zeitreihen[reaktion][v_in][at_node][:1000]
                np.random.seed(10)
                y_label = reaktion + '_node' + str(at_node)  + einheiten[reaktion]

                # Generate a signal
                dt = 0.01
                y = np.random.normal(size=1000) * 25.
                #current_signal=y
                T_ges = len(current_signal) * 0.01
                time = np.arange(0,T_ges, 0.01) # [s]

                # extract the ranges 
                # Find reversals (peaks and valleys), extract cycles and residue (open cycle
                # sequence), process and extract closed cycles from residue.
                reversals, reversals_ix = fatpack.find_reversals(current_signal, k=100)
                cycles, residue = fatpack.find_rainflow_cycles(reversals)
                processed_residue = fatpack.concatenate_reversals(residue, residue)
                cycles_residue, _ = fatpack.find_rainflow_cycles(processed_residue)
                cycles_total = np.concatenate((cycles, cycles_residue))

                #figsize = np.array([140., 140.]) / 25.
                #fig = plt.figure(dpi=96, figsize=figsize)

                # Plotting signal with reversals.
                ax_signal = plt.subplot()
                ax_signal.plot(time, current_signal)
                ax_signal.plot(reversals_ix*dt, current_signal[reversals_ix], 'ro', fillstyle='none', label='reversal')
                ax_signal.scatter(time, current_signal, color='tab:grey',label='Messpunkte',s=3)
                ax_signal.legend()
                ax_signal.set(title="Zeitreihe 3 Seeds", ylabel=y_label , xlabel="Time [s]", xlim=[0, T_ges])# 
                ax_signal.grid()
                plt.show()
