'''
Kopiert von run_model.py am 02.01.2023

Input:
    - Querschnitt
    - Markov matrix als pickle in dictionary Form.
        key1: reaktion (Bezeichnung der Schnittgrößenkomponente oder Spannung) -> Schub kommt direkt als spannung Moment wird erst in kombination mit Normalkraft in spannung umgerechnet
        key2: at_node: Knoten/Ebenen Nummer
        key3: mean, range, N -> arrays in 1 dimensionaler Form jeder eintrag i ist ein Eintrag in der Markov Matrix ij

wichtigeste Optionen: 
    - case: 'druck' oder 'schub' markov daten sind seperat gespeichert und hier erst berechnung der druckspannung aus Schnittgrößen

Berechnung und outputs:
    - Teilschädigung für jeden Eintrag aus der Markov Matrix mit kfat nach DIN oder Züblin (Eingabe)
    - Für jede Ebene werden die Zwischen ergebnisse für jeden Eintrag in ein dictionary geschrieben und dieses mittels dataframe in eine Excel

TODO Bezeichnungen der Kräfte anpassen siehe run model und bem ppt
'''
import numpy as np
from os.path import join as os_join
import pandas as pd
import pickle
import time 

import source.utilities.utilities as utils
import source.utilities.fatigue_utilities as fat
import source.utilities.global_definitions as GD
from source.Querschnitte import KreisRing
from inputs import holz
from source import postprocess


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

PC = 'server' # 'lokal

case = 'druck'#  'schub' #'schub_längs'# 'schub_Fe'# 

fatigue_excel = os_join(*['output','Berechnungs_Ergebnisse_Fatigue_DIN_'+ case + '.xlsx'])

to_excel = True
to_pickle = False
plot_häufigkeit = False
plot_N_R = False
mit_furnier = False

'''
QUERSCHNITTS DEFINITION
NOTE ACHTUNG: muss der selbe sein von dem auch die Schnittgrößen aus dem GZT aufgerufen werden
'''
start_t = time.time()
querschnitts_dateien_laden = False
if not querschnitts_dateien_laden:
    querschnitte,qs_werte, qs_werte_FE, qs_labels  = [],[],[],[]
    nabenhöhe = 130
    dicken = [44]#[36]#[40]# ,, 44]#,48,56,64] # key für lagenaufbau dictonaire in holz 
    
    #mit_furnier = True # dicke 36 sonst 44

    werkstoff_parameter['Furnierebene'] = mit_furnier
    furnier_dict = {True:'BSP_furnier',False:'BSP_normal'}
    #knicke = {'gerade':[None,10], '1Knick':[4.5,10], 'Knick mit Übergang':[4.5,8]}
    knicke = {'Gerade':[None, 11]}
    for label, knick in knicke.items():
        for t in dicken:
                höhen_parameter = {
                                'nabenhöhe' :nabenhöhe,
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

    kreis_ring.compute_effektive_festigkeiten_charakteristisch()
    f_i_k_eff = {'druck':kreis_ring.fc0k_eff_laengs, 'schub':kreis_ring.fvxyk_netto}
    werkstoff_parameter['fc0k_eff_längs'] = round(kreis_ring.fc0k_eff_laengs,2)
    werkstoff_parameter['fvxyk_eff'] = round(kreis_ring.fvxyk_netto,2)
   
    qs_header = pd.MultiIndex.from_product([qs_labels, list(kreis_ring.querschnitts_werte['Ebenen'].keys())])
    #lagen_header = pd.MultiIndex.from_product([['Lagenaufbau [m]'], list(lagen_aufbau[0].keys())])

    qs_werte = np.array(qs_werte).transpose()
    qs_df = pd.DataFrame(qs_werte, columns=qs_header)
    holz_df = pd.DataFrame(werkstoff_parameter, index=[0])
    lagen_df = pd.DataFrame(lagen_aufbau,  index=list(range(1,len(lagen_aufbau)+1)))
    lagen_df.rename({'ti': 'ti [' + kreis_ring.einheiten['Länge'] +']', 'di': 'di [' + kreis_ring.einheiten['Länge']+']'}, axis='columns',inplace=True)

end_t = time.time()
print ('Zeit zum Querschnitt erstellen:', round(end_t-start_t,3))   

if to_excel:
    with pd.ExcelWriter(fatigue_excel, mode= 'w', engine="openpyxl") as writer:#, if_sheet_exists='replace'
        holz_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=0, startcol=0, index=False)#
        qs_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=4, startcol=0)
        lagen_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=qs_df.shape[0] + 10, startcol=0)

    utils.zellen_groeße_formatieren(fatigue_excel, worksheet= 'QS_Werte', cell_width=13, n_cols=max(len(qs_df.columns), len(holz_df.columns))+1)

'''
ÄUßERE BELASTUNG FATIGUE + mittlere Turmbelastung TODO: Alternative SGR Dynmaisch berechnen
NOTE: hier wird jetzt mal einfach nur von einem Querschnitt ausgegangen und nicht mehreren
Einheiten handling hier nicht so eindeutig 
'''
# Datei die der output aus rainflow_FAST_beam.py ist
source_base = {'server':os_join(*['output','Markov']), 'lokal':os_join(*['C:\\','Users','jz','Documents','WEA_lokal','Markov'])}

#markov_file = 
markov_file = os_join(*[source_base[PC],'markov_matrix_'+case+'_DLC12.pkl'])
SGR_statisch_file = os_join(*['output','Berechnungs_Ergebnisse_dataframes','230105_SGR_Spannungen_Ebenen.pkl']) # für vorspannng etc.
output_folder = os_join(*['output','Palmgren_Miner'])

start_t = time.time()
with open(markov_file, 'rb') as handle:
    markov_matrix_z = pickle.load(handle) # Dict mit 1. key die Komponenten 2.key knoten nummer, values sind dann entsprechende Matritzen
if case == 'druck':
    with open(SGR_statisch_file, 'rb') as handle:
        schnittgrößen_GZT = pickle.load(handle)
end_t = time.time()
print ('Zeit zum Öffnen der markov pickles:', round(end_t-start_t,3))   

'''
Die Markov daten in Dataframes sammeln und anhand dessen Spannungen und Schädigung berechnen
für jede Schnittgröße und jeden Knoten einen df
- erst ein dictionary dann daraus direkt ein df
'''
spannung_greek = {  'Mz':GD.GREEK_UNICODE['sigma'], 'Mx':GD.GREEK_UNICODE['tau'], 'Qy':GD.GREEK_UNICODE['tau'], 
                    'tau':GD.GREEK_UNICODE['tau'], 'tau_Fe':GD.GREEK_UNICODE['tau'], 'tau_längs':GD.GREEK_UNICODE['tau']}

#Normalkraft 'ständig' Eigengewicht und Vorspannung aus GZT Berechnung verwenden
# results_df.loc[:,(nabenhöhe, QS_label, 'N [MN]')]
if case == 'druck':
    Nx_ständig = schnittgrößen_GZT[nabenhöhe][kreis_ring.name]['N [MN]'].values * utils.unit_conversion('MN','kN')

fatigue_dicts = {}
aktuelle_einheiten = {'Mz':'kNm', 'Mx':'kNm', 'Qy':'kN','kraft':'kN','fläche':'m²'} #TODO etwas händisch und nicht variable hier
einheit_spannung = utils.unit_conversion('kN/m²', 'N/mm²')
einheit_festigkeit = kreis_ring.einheiten['Festigkeit']

kreis_ring.compute_effektive_festigkeiten_charakteristisch()
f_i_k_eff = {'druck':kreis_ring.fc0k_eff_laengs, 'schub':kreis_ring.fvxyk_netto, 
             'schub_Fe':kreis_ring.fvFek, 'schub_längs':kreis_ring.fvk_brutto}#brutto, da spannung schon nur auf längs lagen bezogen ist

r_2 = 4 # rundungs ziffern
gesamt_ergebniss = {'Höhe [m]':[], 'D_ges':[]}
case_dict = {'druck':'Mz', 'schub':'tau', 'schub_Fe':'tau_Fe','schub_längs':'tau_längs'}
zuordnung = dict((v,k) for k,v in case_dict.items())

for reaktion in [case_dict[case]]: #['Mz']:#['tau']:#, 'Mx', 'Qy']:
    s = spannung_greek[reaktion]
    fatigue_dicts[reaktion] = {}
    for at_node in range(kreis_ring.n_ebenen):#range(1)[9]:# 
        if case == 'druck':
            fatigue_dicts[reaktion][at_node] = {
            s + '_mean_Mz-[N/mm²]':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            s + '_mean_Nx-[N/mm²]':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size))}

        fatigue_dicts[reaktion][at_node].update({
            s + '_mean-[N/mm²]':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            s + '_range-[N/mm²]':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            s + '_max-[N/mm²]':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            s + '_min-[N/mm²]':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            'R':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            'kfat':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            'N_ist':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            'N_ertragbar':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
            'Di':np.zeros((markov_matrix_z[reaktion][at_node]['N'].size)),
        })

        n_mean = len(markov_matrix_z[reaktion][at_node]['mean'])
        n_range = len(markov_matrix_z[reaktion][at_node]['range'])

        D_ges = 0
        zeile = 0

        for i, mean_val in enumerate(markov_matrix_z[reaktion][at_node]['mean']):
            if reaktion == 'Mz':
                # Umrechnung der SGR in Sannunge - gibt nur die druckspannung zurück
                mean_val, sigma_Mz, sigma_Nx = kreis_ring.calculate_normalspannung_at_knoten({'Mz':mean_val, 'Nx':Nx_ständig[at_node]}, at_node, knoten_typ= 'Ebenen') 
            for j, range_val in enumerate(markov_matrix_z[reaktion][at_node]['range']):
                Nij = markov_matrix_z[reaktion][at_node]['N'][j,i]
                if Nij == 0:
                    continue
                if reaktion == 'Mz':
                    range_val, range_Mz, range_Nx = kreis_ring.calculate_normalspannung_at_knoten({'Mz':range_val, 'Nx':0}, at_node, knoten_typ= 'Ebenen')

                sigma_1 = (mean_val - abs(range_val)/2) * einheit_spannung # beide negative daher der wahre maximale wert 
                sigma_2 = (mean_val + abs(range_val)/2) * einheit_spannung # bei druck der absolute maximale werte
                if sigma_1 == 0 and sigma_2 == 0:
                    continue
                # Maximal ist der größere wert unter Berücksichtigung des Vorzeichens
                R, sigma_max_abs, sigma_max, sigma_min = fat.get_R_and_abs_max(sigma_1, sigma_2) # je nach vorzeichen von mean ist s1 oder s2 maximal 
                kfat = abs(sigma_max_abs) / f_i_k_eff[case]
                # DIN
                lgN_ertragbar, N_ertragbar = fat.lgN_R_kfat_din(kfat, R, last=case)
                # ZÜ
                #lgN_ertragbar, N_ertragbar = fat.lgN_R_kfat_zü(kfat, R, last=case)
                Di = Nij/N_ertragbar
                D_ges += Di

                if case == 'druck':
                    fatigue_dicts[reaktion][at_node][s + '_mean_Mz-[N/mm²]'][zeile] = round(sigma_Mz * einheit_spannung, r_2)
                    fatigue_dicts[reaktion][at_node][s + '_mean_Nx-[N/mm²]'][zeile] = round(sigma_Nx * einheit_spannung, r_2)

                fatigue_dicts[reaktion][at_node][s + '_mean-[N/mm²]'][zeile] = round(mean_val * einheit_spannung, r_2)
                fatigue_dicts[reaktion][at_node][s + '_range-[N/mm²]'][zeile] = round(abs(range_val) * einheit_spannung, r_2)
                fatigue_dicts[reaktion][at_node]['N_ist'][zeile] = Nij

                fatigue_dicts[reaktion][at_node][s + '_max-[N/mm²]'][zeile] = round(sigma_max, r_2)
                fatigue_dicts[reaktion][at_node][s + '_min-[N/mm²]'][zeile] = round(sigma_min, r_2)
                fatigue_dicts[reaktion][at_node]['R'][zeile] = round(R,r_2)
                fatigue_dicts[reaktion][at_node]['kfat'][zeile] = round(kfat,r_2)
                fatigue_dicts[reaktion][at_node]['N_ertragbar'][zeile] = N_ertragbar
                fatigue_dicts[reaktion][at_node]['Di'][zeile] = Di

                zeile += 1
        
        # die restlichen einträge der initalisierten np.zeros sind 0 und werden abgeschnitten
        for key, arr in fatigue_dicts[reaktion][at_node].items():
            fatigue_dicts[reaktion][at_node][key] =arr[:zeile]


        if plot_häufigkeit:
            var = 'N_ist' #'R' # 
            postprocess.plot_häufigkeit(fatigue_dicts[reaktion][at_node][var], var=var, show_plots=True, save_plots=True,
                                        fig_name_save = os_join(*['output', 'plots','Rainflow','Häufigkeiten', var + '_'+ zuordnung[reaktion] + '_at_node' + str(at_node) + '.png']), 
                                        title= zuordnung[reaktion]  + ' at node' + str(at_node))

        if plot_N_R:
            bin_size_R = 0.1
            N_R, bin_centers = fat.N_von_R(fatigue_dicts[reaktion][at_node], bin_size=bin_size_R)
            postprocess.plot_N_von_R(N_R, bin_centers, show_plots=False, save_plots=True,
                                fig_name_save = os_join(*['output', 'plots','Rainflow','N_R', zuordnung[reaktion]  + '_at_node' + str(at_node) + '.png']), 
                                title= zuordnung[reaktion]  + ' at node' + str(at_node))

        df = pd.DataFrame(fatigue_dicts[reaktion][at_node])

        if to_pickle:
            outfile = os_join(*[output_folder, reaktion + '_Ebene' + str(at_node) + '_DLC12.pkl'])
            df.to_pickle(outfile)

        if to_excel:
            with pd.ExcelWriter(fatigue_excel, mode= 'a', engine="openpyxl", if_sheet_exists='replace') as writer:
                df.to_excel(writer, sheet_name= case + '_Ebene' + str(at_node)  , startrow=0, startcol=0, index=False)#

            utils.zellen_groeße_formatieren(fatigue_excel, worksheet= case + '_Ebene' + str(at_node), df=df, n_cols=len(df.columns)+1)

        print('Auf Höhe:', kreis_ring.section_absolute_heights['Ebenen'][at_node], 'm beträgt die gesamt Schädigung für', case, D_ges)
        gesamt_ergebniss['Höhe [m]'].append(kreis_ring.section_absolute_heights['Ebenen'][at_node])
        gesamt_ergebniss['D_ges'].append(D_ges)

    if to_excel:
        gesamt_df = pd.DataFrame(gesamt_ergebniss)
        with pd.ExcelWriter(fatigue_excel, mode= 'a', engine="openpyxl", if_sheet_exists='replace') as writer:
            gesamt_df.to_excel(writer, sheet_name= case + '_D_gesamt', startrow=0, startcol=0, index=False)#

        utils.zellen_groeße_formatieren(fatigue_excel, worksheet= case + '_D_gesamt', df=df, n_cols=len(df.columns)+1)

        from datetime import datetime
        now = datetime.now()
        dt_string = now.strftime("%d/%m/%Y %H:%M")
        print ('\nErgebnisse in ', fatigue_excel, 'geschrieben -', dt_string)

    

