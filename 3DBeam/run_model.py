'''
Das kommt von der Jupyter Version vom 26.10.2022
TODO:
    - run ausdünnen und übersichtlichere, am Anfang stehende Parameter Eingabe
    - dataframe und in excel ausgabe sachen in einen writer
    - Lastenfiles erzeugung auslagern -> Lasten Klasse
    - Querschnitt und Höhen ggf. besser trennen
    - ...
'''

import numpy as np
import os
from os.path import join as os_join
import pandas as pd

from source.model import BeamModel
import source.postprocess as postprocess
import source.utilities.utilities as utils
from source.Querschnitte import KreisRing
import source.utilities.global_definitions as GD
from inputs import holz
from inputs.spannglieder import Spannglieder


parameters_init = {
                'dimension': '3D',
                'n_elements': 10, #
                'nacelle_mass': 14400, # kg Kleinwind,# IEA37:267910,#  Gondel Masse in kg # aus optimierung 287920.5 
                'imperfektion':0.008, # m/m Schiefstellung + Fundament schief
                'E_Modul': 12000E+06,# N/m²
                'type_of_bc':'Eingespannt',#'Feder', #'clamped',# TODO Feder evlt falsch gemacht
                'spring_stiffness':[1E+13,2E+13], # Federsteifigkeit am Boden in u und gamma richtung Bögl 40 GNm/rad
            }

parameters_init = utils.add_defaults(parameters_init)

holzgüte = 'C30'

# # VORSPANNUNG ANGABEN
n_draht_suspa_ex, stahl_supsa = 20, 'St_1570/1770' # 20 # Nabenhöhe > 100: 84
n_mono_litzen, stahl_mono = 2, 'St_1860' # Nabenhöhe > 100: 5
# n_segmente_ext: Anzahl der Segmente die durch Externe Vorspannung Überdrückt sein sollen: Inklusive der Fuge Übergnag stahl holz
vorpsannungs_params = {'n_segmente_ext':3,  'Pd_ext':Spannglieder.suspa_draht_ex[stahl_supsa]['Pm0_n'][n_draht_suspa_ex], 
                                            'Pd_int':Spannglieder.monolitzen[stahl_mono]['Pmax_n'][n_mono_litzen]}
spannglieder_werte_df = pd.DataFrame(Spannglieder.parameter_dict(n_draht_suspa_ex, n_mono_litzen, stahl_supsa, stahl_mono))

n_nodes = parameters_init['n_elements']+1
werkstoff_parameter = holz.charakteristische_werte[holzgüte]
parameters_init['material_density'] = werkstoff_parameter['rhok']
parameters_init['E_Modul'] = werkstoff_parameter['E0mean']

nachweis_parameter = holz.HOLZBAU

spannkraft_verlust_pauschal = 20 # % 

results_excel = os_join(*['output','Berechnungs_Ergebnisse.xlsx'])

'''
QUERSCHNITTS DEFINITION
'''
knoten_typen = ['Ebenen','FE']
querschnitts_dateien_laden = False
if not querschnitts_dateien_laden:
    querschnitte,qs_werte, qs_werte_FE, qs_labels  = [],[],[],[]
    höhe = 35
    mit_furnier = True # t_furnier = 38,40,44; BSP mit dem Druck Nachweis aufgeht: t 44

    dicken = [20]#[44]# [44]#[40]#,, 44]#,48,56,64] # key für lagenaufbau dictonaire in holz 

    werkstoff_parameter['Furnierebene'] = mit_furnier
    furnier_dict = {True:'BSP_furnier',False:'BSP_normal'}
    #knicke = {'gerade':[None,10], '1Knick':[4.5,10], 'Knick mit Übergang':[4.5,8]}

    knicke = {'Gerade':[None, 6]}#{'Gerade':[None, 11]}
    for label, knick in knicke.items():
            for t in dicken:
                    höhen_parameter = {
                                    'nabenhöhe' :höhe,
                                    'd_unten_oben' :[3.0, 1.8],#[11, 3.4],
                                    'höhe_sektionen' :[11,12], # Range angeben
                                    'd_knick' :knick[0],#None, für keinen Knick -> alle folgenden Parameter sind dann unrelevant
                                    'h_knick_von_oben' :70,
                                    'd_unten_angepasst':knick[1], # damit ein Übergangsknick entseht muss dieser hier kleiner als 'd_unten' sein
                                    'd_knick_übergang':'automatisch',#7.5,# # bei automatisch nimmt es den schon vorhandenen an dieser Stelle
                                    'n_sektionen_übergang':1,
                                    #'transportbreite_max':3, # past nicht so ganz zur höhe aber mus irgendwo hin
                                    } 

                    #lagen_aufbau = holz.lagenaufbauten[furnier_dict[mit_furnier]][t]
                    lagenaufbauten = holz.Lagenaufbauten(furnier_dict[mit_furnier], t)
                    lagen_aufbau = lagenaufbauten.get()
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
   
    qs_header = pd.MultiIndex.from_product([qs_labels, list(kreis_ring.querschnitts_werte['Ebenen'].keys())], names=['', 'Ebene'])
    qs_header_FE = pd.MultiIndex.from_product([qs_labels, list(kreis_ring.querschnitts_werte['FE'].keys())], names=['', 'Knoten'])
    #lagen_header = pd.MultiIndex.from_product([['Lagenaufbau [m]'], list(lagen_aufbau[0].keys())])

    qs_werte=np.array(qs_werte).transpose()
    qs_werte_FE=np.array(qs_werte_FE).transpose()
    qs_df = pd.DataFrame(qs_werte, columns=qs_header)
    qs_FE_df = pd.DataFrame(qs_werte_FE, columns=qs_header_FE)
    holz_df = pd.DataFrame(werkstoff_parameter, index=[0])
    holz_units_df = pd.DataFrame(holz.charakteristische_werte['units'])
    lagen_df = pd.DataFrame(lagen_aufbau,  index=list(range(1,len(lagen_aufbau)+1)))
    lagen_df.rename({'ti': 'ti [' + kreis_ring.einheiten['Länge'] +']', 'di': 'di [' + kreis_ring.einheiten['Länge']+']'}, axis='columns',inplace=True)
    lagen_df.index.name = 'Lage'
    
    # Querschnittswerte speichern für openfast einlesen -> NOTE gehe hier davon aus das gerade nur ein QS
    qs_FE_df.to_pickle(os_join(*['output','Turmwerte_OpenFAST', qs_label_openfast + '.pkl']))
    print ('  Querschnittswerte in dictionary als pickle gespeichert in:', os_join(*['output','Turmwerte_OpenFAST', qs_label_openfast + '.pkl']))

'''
DATAFRAME
'''
# https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy

einheiten_out = {'Kraft':'kN','Moment':'kNm','Spannung':'N/mm²'} # TODO das als einheiten rechner von berechnung und ausgabe nutzen
df_results_header_list = [[],[], []] # 1. Liste = 1. Level usw...
kraft_komponenten =  []#['Fx [kN]', 'Fy [kN]', 'Mz [kNm]','Mx [kNm]'] #TODO hier auch die einheiten zentraler steuern
df_einwirkung_header_list = [[],[],['ständig','kurz','egal'], kraft_komponenten]
df_einwirkung_typ_header_list = [[],[],['Kopflast [kN]', 'windkraft [kN]', 'eigengewicht [kN]', 'F_h,ers,imp [kN]' ],kraft_komponenten]
df_vorspannung_header_list = [[],[],['Höhe [m]', 'd_achse [m]', 'P_erf,end [MN]', 'n EX-'+str(n_draht_suspa_ex), 'n Mono-'+str(n_mono_litzen), 'n int. Summe', 'P_ist Fuge [MN]']]
df_bauzustand_header_list = [[],[],[],[]]
include_tau = True
include_sigma_M_N = False
include_sgr = True
include_reibung = True
include_bauzustand_sgr = True
include_bauzustand_spannung = False
include_bauzustand_sgr_m = True # pro meter

for qs in querschnitte:
    if qs.nabenhöhe not in df_results_header_list[0]:
        df_results_header_list[0].append(qs.nabenhöhe)
        df_einwirkung_header_list[0].append(qs.nabenhöhe)
        df_bauzustand_header_list[0].append(qs.nabenhöhe)
    if qs.name not in df_results_header_list[1]:
        df_results_header_list[1].append(qs.name)
        df_einwirkung_header_list[1].append(qs.name)
        df_bauzustand_header_list[1].append(qs.name)

# Berechnungsergebnisse - dritte Ebene:
df_results_header_list[2].append('Höhe [' + qs.einheiten['Länge'] + ']') 
df_results_header_list[2].append('d_achse [' + qs.einheiten['Länge'] + ']') 
df_results_header_list[2].append('P_erf,end [MN]') 
df_results_header_list[2].append('P_ist [MN]') 

if include_sgr:
    df_results_header_list[2].append('Mz [MNm]') 
    df_results_header_list[2].append('N [MN]') 
    df_results_header_list[2].append('G [MN]') 
    df_results_header_list[2].append('Q [MN]') 
    df_results_header_list[2].append('Mx [MNm]') 
if include_sigma_M_N:
    df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_N [' + qs.einheiten['Normalspannung'] + ']') 
    df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_M [' + qs.einheiten['Normalspannung'] + ']')  
df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_max [' + qs.einheiten['Normalspannung'] + ']')  
df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_min [' + qs.einheiten['Normalspannung'] + ']')  
df_results_header_list[2].append('Ausn. druck') 

if include_tau:
    if werkstoff_parameter['Furnierebene']:
        df_results_header_list[2].append(GD.GREEK_UNICODE['tau'] + '_Fe' + ' [' + qs.einheiten['Schubspannung'] + ']') 
        df_results_header_list[2].append(GD.GREEK_UNICODE['tau'] + '_längs' + ' [' + qs.einheiten['Schubspannung'] + ']')
        df_results_header_list[2].append('Ausn. schub Furnier') 
        df_results_header_list[2].append('Ausn. schub längs')  
    else:
        df_results_header_list[2].append(GD.GREEK_UNICODE['tau'] + '_Qy' + ' [' + qs.einheiten['Schubspannung'] + ']') 
        df_results_header_list[2].append(GD.GREEK_UNICODE['tau'] + '_Mx' + ' [' + qs.einheiten['Schubspannung'] + ']') 
        df_results_header_list[2].append(GD.GREEK_UNICODE['tau'] + '_xy' + ' [' + qs.einheiten['Schubspannung'] + ']') 
        df_results_header_list[2].append(GD.GREEK_UNICODE['tau'] + '_tor' + ' [' + qs.einheiten['Schubspannung'] + ']') 
        df_results_header_list[2].append('Ausn. schub brutto') 
        df_results_header_list[2].append('Ausn. schub netto') 
        df_results_header_list[2].append('Ausn. schub torsion') 

if include_reibung:
    df_results_header_list[2].append('nxy_Qy' + ' [kN/m]') 
    df_results_header_list[2].append('nxy_Mx' + ' [kNm/m]') 
    df_results_header_list[2].append('nxy_P,Rd' + ' [kN/m]') 
    df_results_header_list[2].append('Ausn. reibung') 

for segment in range(1,qs.n_ebenen):
    seg = 'seg_0-' + str(segment)
    df_bauzustand_header_list[2].append(seg)
if include_bauzustand_sgr:
    df_bauzustand_header_list[3].extend(['Mz [MNm]', 'N [MN]', 'Q [MN]','Mx [MNm]', 'Reibung [MN]', 'Ausn. reibung'])
if include_bauzustand_spannung:
    df_bauzustand_header_list[3].append(GD.GREEK_UNICODE['sigma'] + '_max [' + qs.einheiten['Normalspannung'] + ']') 
    df_bauzustand_header_list[3].append(GD.GREEK_UNICODE['sigma'] + '_min [' + qs.einheiten['Normalspannung'] + ']') 
if include_bauzustand_sgr_m:
    df_bauzustand_header_list[3].append('nz_max [kN/m]')
    df_bauzustand_header_list[3].append('nz_min [kN/m]')

# Header in pd Format
df_results_header = pd.MultiIndex.from_product(df_results_header_list, names=['Nabenhöhe', 'Querschnitt', 'Segment'])
df_bauzustand_header = pd.MultiIndex.from_product(df_bauzustand_header_list, names=['Nabenhöhe', 'Querschnitt', 'Segment', 'Komponente'])
df_einwirkungs_header = pd.MultiIndex.from_product(df_einwirkung_header_list, names=['Nabenhöhe', 'Querschnitt', 'Ew. Dauer', 'Komponente', ])
df_einwirkungs_header_typ = pd.MultiIndex.from_product(df_einwirkung_typ_header_list, names=['Nabenhöhe', 'Querschnitt', 'Lasttyp','Komponente', ])
df_vorspannung_header = pd.MultiIndex.from_product(df_vorspannung_header_list, names=['Nabenhöhe', 'Querschnitt', 'Segment'])
# Dataframes leer 
results_df = {'max_druck':pd.DataFrame(columns=df_results_header),'max_schub':pd.DataFrame(columns=df_results_header)}#pd.DataFrame(columns=df_results_header)#
einwirkungs_parameter_df = {}
einwirkungs_df = {'max_druck':pd.DataFrame(columns=df_einwirkungs_header),'max_schub':pd.DataFrame(columns=df_einwirkungs_header)}
einwirkungs_df_typ = {'max_druck':pd.DataFrame(columns=df_einwirkungs_header_typ), 'max_schub':pd.DataFrame(columns=df_einwirkungs_header_typ)}

# Unabhängig von der Belastung
if include_bauzustand_sgr:
    bauzustands_df = pd.DataFrame(columns=df_bauzustand_header)
vorpsannungs_df = pd.DataFrame(columns=df_vorspannung_header)

'''
ÄUßERE BELASTUNG
'''

einwirkungsdauer = ['ständig', 'kurz', 'egal','spannkraft']
sicherheitsbeiwerte = {'dlc':1.35, 'wind':1.35, 'g':1.35, 'vorspannung':1.0, 'günstig':1.0} # dlc entsprechend der Kopflasten q diesem entsprechend und g TODO ansich komplexer im IEC

# TODO lastfall anzahl variabel halten siehe schon bei der erstellung der dfs sehr händisch und abhängig von key wörtern verschiedene!!
lastfälle = {'max_druck': '@max_Fxy','max_schub':'@max_all'}#'@max_Mz', 'max_all':} 

# sind schon designlasten Cosy noch IEA
kopf_lasten_IEA = { '@max_Fxy':{'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':5.98E+06, 'My':6.81E+06, 'Mz':2.31E+06},# NOTE My = Mxy #Masse charakt. 2.62E+06, Fx = Fxy 
                    '@max_Mz':{'Fx':419925,'Fy':-24468, 'Fz':-3650438, 'Mx':6294993, 'My':-3439939, 'Mz':10811885}, # Fx = Fxy, My = Mxy
                    '@max_all':{'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':5.98E+06, 'My':6.81E+06, 'Mz':10811885}} # bisher max Fx und max Mz 

qs.berechnungs_hinweise.append('   - Kopflast Fy (beam cosy) entspricht der resultierenden Fxy aus den IEA Lasten')
qs.berechnungs_hinweise.append('   - Kopflast Mz (beam cosy) entspricht der resultierenden Mxy aus den IEA Lasten')

skalierung_IEA_kleinwind = 707 / 13273 *2# Rotorflächen x2
kopf_lasten_IEA = utils.scale_dict_of_dicts(kopf_lasten_IEA, skalierung_IEA_kleinwind)

qs.berechnungs_hinweise.append('   - Kopflasten IEA skaliert anhand der Rotorflächen mit  (und dann x 2)' + str(round(skalierung_IEA_kleinwind,5)))

kopf_masse_design = - parameters_init['nacelle_mass'] * GD.GRAVITY * sicherheitsbeiwerte['g'] # -2628197

import inputs.DIN_Windlasten as wind_DIN

terrain_kategorie = 'II'
windzone = 2
# 25 m/s = cutout windspeed
basis_windgeschwindigkeit = wind_DIN.vb_von_v_nabenhöhe(25, terrain_kategorie, qs.nabenhöhe) #17

vb_bauzustad = wind_DIN.VB_WINDZONEN[windzone] # Windzone 2 vb = 25 m/s TODO heir extra windlast berechne und für den bauzustand verwenden

nachweis_parameter_dict = { GD.GREEK_UNICODE['gamma']+'_m': [nachweis_parameter['gamma_m']],
                            'kmod Ew. kurz': [nachweis_parameter['k_mod']['kurz']],
                            'kmod Ew. ständig': [nachweis_parameter['k_mod']['lang']],
                            'ksys': [nachweis_parameter['k_sys']],
                            'Reibung; ' + GD.GREEK_UNICODE['mu']: [nachweis_parameter['mu']]
}

nachweis_parameter_df = pd.DataFrame(nachweis_parameter_dict)

print ('Windbelastung aus:')
print ('    vb 10 m:', round(basis_windgeschwindigkeit,2))
print ('    terrain:', terrain_kategorie)
#wind_DIN.plot_DIN_all(basis_windgeschwindigkeit, categories=[terrain_kategorie])

max_ausnutzung = {}
for lastfall in lastfälle:
    last_at = lastfälle[lastfall]
    print ('\nKopflast von IEA:\n  ', last_at)

    querschnitte[0].berechnungs_hinweise.append('   - Last am Kopf aus IEA ' + last_at)
    einwirkungs_parameter = {'vb':round(basis_windgeschwindigkeit,2),'Terrain Kategorie':terrain_kategorie, 'cd': cd_zylinder, 
                            'Kopflast':last_at, 'DLC':'1.3_seed2_9ms', 'Kosy.':'Balken'}

    einwirkungs_parameter.update(sicherheitsbeiwerte)

    einwirkungs_parameter_df[lastfall] = pd.DataFrame(einwirkungs_parameter, index=[0])

    kopf_lasten_beam = utils.convert_coordinate_system_and_consider_einwirkungsdauer(kopf_lasten_IEA[last_at], n_nodes, sicherheitsbeiwerte['g'], 
                                                                                    kopf_masse = kopf_masse_design, dimension=parameters_init['dimension']) 

    wind_kraft_z = {}
    lasten_files, lasten_dicts_dauer, lasten_dicts_typ = {}, {}, {}
    for qs in querschnitte:
        QS_label = qs.name
        lasten_dict_base = {fi:np.zeros(n_nodes) for fi in GD.FORCES[parameters_init['dimension']]}

        # Ergebniss Daten vorbereiten -> leer initialisieren nur wenn noch nicht leer 
        if QS_label not in lasten_files:
            lasten_files[QS_label] = {}
            lasten_dicts_dauer[QS_label], lasten_dicts_typ[QS_label] = {}, {}

        if qs.nabenhöhe not in lasten_files[QS_label]:
            lasten_files[QS_label][qs.nabenhöhe] = {}
            lasten_dicts_dauer[QS_label][qs.nabenhöhe], lasten_dicts_typ[QS_label][qs.nabenhöhe] = {}, {}

        # Lastberechnung für den spezifischen QS
        v_z, Iv_z, qp_z, z = wind_DIN.DIN_potenz_profil(basis_windgeschwindigkeit, terrain_kategorie, qs.nabenhöhe)
        
        knoten_wind_kraft_z = {}
        for knoten_typ in qs.d_achse:
            wind_kraft_z, z_coords = wind_DIN.wind_kraft(vb=basis_windgeschwindigkeit, category=terrain_kategorie, height=qs.section_absolute_heights[knoten_typ], cd = qs.cd, Aref=qs.d_achse[knoten_typ]) # ergibt N/m
            knoten_wind_kraft_z[knoten_typ] = utils.linien_last_to_knoten_last(wind_kraft_z, qs.section_absolute_heights[knoten_typ], gamma_q = sicherheitsbeiwerte['wind']) # NOTE gamma_q nach IEC DLC 1.3
            #ebenen_wind_kraft_z = utils.linien_last_to_knoten_last(wind_kraft_z, qs.section_absolute_heights, gamma_q = sicherheitsbeiwerte['wind']) # NOTE gamma_q nach IEC 
            
        #postprocess.plot_along_height(knoten_wind_kraft_z, z_coords, label='wind_kraft [N]')

        gewichtskraft_design = {'FE':{k: sicherheitsbeiwerte['g']*v for k,v in qs.gewichtskraft['FE'].items()}} # im QS ist die charakteristische gespeichert
        gewichtskraft_design['Ebenen'] = {k: sicherheitsbeiwerte['g']*v for k,v in qs.gewichtskraft['Ebenen'].items()} 

        #_____Lasten sortiert nach Typ der Last
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['kopflast'] = kopf_lasten_beam['egal']
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['windkraft'] = knoten_wind_kraft_z['FE']
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['eigengewicht'] = gewichtskraft_design['FE']
        F_h_imp_ers = utils.horizontale_ersatzlast([gewichtskraft_design['FE'], kopf_lasten_beam['egal']], theta_x = parameters_init['imperfektion'])
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['F_h,ers,imp'] = F_h_imp_ers
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['bauzustand'] = utils.lasten_dict_bauzustand(lasten_dict_base, knoten_wind_kraft_z, qs.gewichtskraft['Ebenen'], sicherheitsbeiwerte,
                                                                                            QS_obj=qs, abminderung_wind=0.5, parameter = parameters_init)

        # _____Sortieren nach Einwirkungsdauer
        lasten_dicts_dauer[QS_label][qs.nabenhöhe]['egal']  = utils.update_lasten_dict(lasten_dict_base, [knoten_wind_kraft_z['FE'], F_h_imp_ers, 
                                                                                                         gewichtskraft_design['FE'], kopf_lasten_beam['egal']]) 
        lasten_dicts_dauer[QS_label][qs.nabenhöhe]['kurz'] = utils.update_lasten_dict(lasten_dict_base, [knoten_wind_kraft_z['FE'], kopf_lasten_beam['kurz']])
        # TODO Ein Teil es Kopfmoments Mz (beam) ist auch ständig!
        lasten_dicts_dauer[QS_label][qs.nabenhöhe]['ständig']  = utils.update_lasten_dict(lasten_dict_base, [gewichtskraft_design['FE'], F_h_imp_ers, kopf_lasten_beam['ständig']])
        # das ist der Lastfall spannkraft -> ist für die Spannkraftberechnung (gamma_m Eigengewicht = 1,0)
        lasten_dicts_dauer[QS_label][qs.nabenhöhe]['spannkraft']  = utils.update_lasten_dict(lasten_dict_base, [knoten_wind_kraft_z['FE'], qs.gewichtskraft['FE'], F_h_imp_ers, kopf_lasten_beam['spannkraft']])
    
        # ____________ LASTEN DATEI GENERIEREN - nur nach dauer sortiert da dies relevant für ausnutzung ist ___________________________
        for dauer in einwirkungsdauer:
            filename = 'K-IEA'+ last_at + '_W-v' +str(round(basis_windgeschwindigkeit,1)) + 'cd' + str(qs.cd) + '_D-' + dauer
            lasten_files[QS_label][qs.nabenhöhe][dauer] = utils.generate_lasten_file(n_nodes, lasten_dicts_dauer[QS_label][qs.nabenhöhe][dauer], 
                                                                                     file_base_name=filename, dimension=parameters_init['dimension'])
        
        # Händische Konstante Torsion dazu wird auch bei den schnittgrößen seperat betrachtet
        if parameters_init['dimension'] == '2D':
            qs.berechnungs_hinweise.append('   - Torsion konstant über die Höhe nur aus Kopflast')
            lasten_dicts_dauer[QS_label][qs.nabenhöhe]['egal']['Mx']  = np.append(np.zeros(parameters_init['n_elements']), kopf_lasten_IEA[last_at]['Mz'])
            lasten_dicts_dauer[QS_label][qs.nabenhöhe]['kurz']['Mx']  = np.append(np.zeros(parameters_init['n_elements']), kopf_lasten_IEA[last_at]['Mz'])
            lasten_dicts_typ[QS_label][qs.nabenhöhe]['kopflast']['Mx'] = np.append(np.zeros(parameters_init['n_elements']), kopf_lasten_IEA[last_at]['Mz'])

        # LASTEN NACH DAUER SORTIERT
        for dauer in einwirkungsdauer:
            if dauer == 'spannkraft':
                continue
            for komponente in lasten_dicts_dauer[QS_label][qs.nabenhöhe][dauer]:
                kategorie = GD.FORCE_CATEGORY[komponente]
                
                einwirkungs_df[lastfall].loc[:,(qs.nabenhöhe, QS_label, dauer, komponente + ' [' + einheiten_out[kategorie] + ']')] =\
                    np.around(lasten_dicts_dauer[QS_label][qs.nabenhöhe][dauer][komponente] *\
                        utils.unit_conversion(qs.einheiten[kategorie],einheiten_out[kategorie]),2)

        # NACH LASTTYP SORTIERT
        for typ in lasten_dicts_typ[QS_label][qs.nabenhöhe]:
            if typ == 'bauzustand':
                qs.berechnungs_hinweise.append('   - Für den Bauzustand wird die Windlast mit einem Abminderungsbeiwert von 0.5 berücksichtigt (DIN 1055)')
                qs.berechnungs_hinweise.append('   - Für den Bauzustand wird das Eigengewicht gamma_q = 1.0 berücksichtigt (Für Berechnung der Schiefstellung mit 1.35)')
                lasten_files[QS_label][qs.nabenhöhe][typ] = {}
                for segment in lasten_dicts_typ[QS_label][qs.nabenhöhe][typ]:
                    filename = 'K-IEA'+ last_at + '_W-v' +str(round(basis_windgeschwindigkeit,1)) +'cd' + str(qs.cd) + '_D-' + typ + '_' + segment
                    lasten_files[QS_label][qs.nabenhöhe][typ][segment] = utils.generate_lasten_file(qs.n_ebenen, lasten_dicts_typ[QS_label][qs.nabenhöhe][typ][segment], 
                                                                                                    file_base_name=filename, dimension=parameters_init['dimension'])
            
            else:
                filename = 'K-IEA'+ last_at + '_W-v' +str(round(basis_windgeschwindigkeit,1)) + 'cd' + str(qs.cd) + '_D-' + typ
                lasten_files[QS_label][qs.nabenhöhe][typ] = utils.generate_lasten_file(n_nodes, lasten_dicts_typ[QS_label][qs.nabenhöhe][typ], 
                                                                                        file_base_name=filename, dimension=parameters_init['dimension'])
                for komponente in lasten_dicts_typ[QS_label][qs.nabenhöhe][typ]:
                    kategorie = GD.FORCE_CATEGORY[komponente]

                    einwirkungs_df_typ[lastfall].loc[:,(qs.nabenhöhe, QS_label, typ, komponente + ' [' + einheiten_out[kategorie] + ']')] =\
                        np.around(lasten_dicts_typ[QS_label][qs.nabenhöhe][typ][komponente] *\
                            utils.unit_conversion(qs.einheiten[kategorie],einheiten_out[kategorie]),2)

        if qs.nabenhöhe == 30:
            postprocess.plot_dict_subplots(lasten_dicts_typ[QS_label][qs.nabenhöhe]['windkraft'], qs.x_FE, title='Windkraft ' + QS_label, unit='kN')
            postprocess.plot_dict_subplots(lasten_dicts_dauer[QS_label][qs.nabenhöhe]['egal'], qs.x_FE, title='Windkraft + Kopflasten gesamt ' + QS_label, unit = 'kN')

        plot_nach_dauer = False
        if plot_nach_dauer:
            for dauer in einwirkungsdauer:
                postprocess.plot_dict_subplots(lasten_dicts_dauer[QS_label][qs.nabenhöhe][dauer], qs.x_FE, title='Last ' + dauer +' ' + QS_label, unit = 'kN')


    '''
    SCHNITTGRÖßEN
    '''
    # Werden getrennt an den Stellen der Ebenen und den FE Knoten ausgewertet
    schnittgrößen_design = {'Ebenen':{},'FE':{},'bauzustand':{}}
    bauzustands_kram = {}

    for querschnitt in querschnitte:
        QS_label = querschnitt.name
        nabenhöhe = querschnitt.nabenhöhe

        if QS_label not in bauzustands_kram:
            bauzustands_kram[QS_label] = {}
        if nabenhöhe not in bauzustands_kram[QS_label]:
            bauzustands_kram[QS_label][nabenhöhe] = {}

        for knoten_typ in schnittgrößen_design:
            if QS_label not in schnittgrößen_design[knoten_typ]:
                schnittgrößen_design[knoten_typ][QS_label] = {}
            if nabenhöhe not in schnittgrößen_design[knoten_typ][QS_label]:
                schnittgrößen_design[knoten_typ][QS_label][nabenhöhe] = {}

        section_properties = querschnitt.section_parameters 

        parameters = {'FE': utils.add_model_data_from_dict(section_properties['FE'], parameters_init)}
        parameters['Ebenen'] = utils.add_model_data_from_dict(section_properties['Ebenen'], parameters_init)
        parameters['Ebenen']['n_elements'] = querschnitt.n_sections
        if parameters_init['type_of_bc'] == 'Feder':
            querschnitt.berechnungs_hinweise.append('   - Einspannung am Fuß mit Federsteifigkeiten u & Drehfder gamma: ' +  str(parameters_init['spring_stiffness']))
        elif parameters_init['type_of_bc'] == 'Eingespannt':
            querschnitt.berechnungs_hinweise.append('   - Einspannung am Fuß starr')

        beam = BeamModel(parameters['FE'], adjust_mass_density_for_total = False, optimize_frequencies_init=False , apply_k_geo=False)
        beam_ebenen = BeamModel(parameters['Ebenen'], adjust_mass_density_for_total = False, optimize_frequencies_init=False , apply_k_geo=False)
        
        print ('\nN_nodes', querschnitt.n_ebenen)
        print (QS_label, nabenhöhe, 'm')
        print ('     Gesamt Volumen des Querschnitts [m³]:', round(beam.volume,2))
        print ('     Gesamt Gewichtskraft Turm am Fuß [MN]:',round(sum(beam.eigengewicht * GD.UNIT_SCALE['MN']),2))
        print ('     Frequenzen [Hz]:', [round(beam.eigenfrequencies[i],3) for i in range(3)])

        lasten_nach_dauer = {'ständig':{'torsion':False},
                             'egal':{'torsion':kopf_lasten_IEA[last_at]['Mz']},
                             'kurz':{'torsion':kopf_lasten_IEA[last_at]['Mz']},
                             'spannkraft':{'torsion':kopf_lasten_IEA[last_at]['Mz']}}

        for dauer in einwirkungsdauer:
            # NOTE Torsion ist nicht in der Steifigkeitsmatrix, deswegen kann diese brechnung nur händisch ergänzt werden
            lasten_file = lasten_files[QS_label][nabenhöhe][dauer]
            schnittgrößen_design['FE'][QS_label][nabenhöhe][dauer] = {}
            schnittgrößen_design['Ebenen'][QS_label][nabenhöhe][dauer] = {}

            schnittgrößen_design['FE'][QS_label][nabenhöhe][dauer] = beam.static_analysis_solve(load_vector_file=lasten_file, 
                                                                        constant_torsion = lasten_nach_dauer[dauer]['torsion'])

            # SGR werte auch an Ebenen knoten bekommen -> Linear interpolieren 
            for key in schnittgrößen_design['FE'][QS_label][nabenhöhe][dauer]:
                schnittgrößen_design['Ebenen'][QS_label][nabenhöhe][dauer][key] = utils.interpolate(z_soll = querschnitt.section_absolute_heights['Ebenen'], 
                                                                                                    z_ist=beam.nodal_coordinates['x0'], 
                                                                                                    werte = schnittgrößen_design['FE'][QS_label][nabenhöhe][dauer][key])

            lasten_label = os.path.splitext(os.path.basename(lasten_files[QS_label][nabenhöhe][dauer]))[0]                                                                                        
                    
        print ('     Maximales Moment [MNm]:', round(max(abs(schnittgrößen_design['Ebenen'][QS_label][nabenhöhe]['egal']['g'])) * utils.unit_conversion('Nm', 'MNm'),2))   
        print ('     Maximales Moment ohne Kopflast [MNm] ~=', round((max(abs(schnittgrößen_design['Ebenen'][QS_label][nabenhöhe]['egal']['g'])) - kopf_lasten_IEA[last_at]['Fx']*nabenhöhe) * utils.unit_conversion('Nm', 'MNm'),2))   

        # BAUZUSTAND 
        if include_bauzustand_sgr:
            for segment in lasten_files[QS_label][nabenhöhe]['bauzustand']:
                schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]= {}
                bauzustands_kram[QS_label][nabenhöhe][segment] = {}
                lasten_file = lasten_files[QS_label][nabenhöhe]['bauzustand'][segment]
                schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment] = beam_ebenen.static_analysis_solve(load_vector_file=lasten_file, return_result=True)
                bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung']  = -1 * utils.remove_small_values_from_array(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['x']) * nachweis_parameter['mu']
                bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung Ausn.'] = np.divide(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['y'],
                                                                                            bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung'],
                                                                                            out=np.zeros_like(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['y']), 
                                                                                            where=bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung']!=0)

        plot_sgr = False                                                            
        if plot_sgr:
            lasten_label = os.path.splitext(os.path.basename(lasten_files[QS_label][nabenhöhe]['egal']))[0]
            postprocess.plot_static_result_forces(beam, 'internal_forces', ['x','y','g','a'], unit='MN', title_suffix='Last:' + lasten_label, figsize_scale = 1.1)
            postprocess.plot_static_result(beam, 'deformation', ['y'], unit='cm',title_suffix='Last:' + lasten_label,)
            #postprocess.plot_static_result(beam, 'deformation', ['g'], unit='cm',title_suffix='Last:' + lasten_label,)

    '''
    SPANNUNGEN UND AUSNUTZUNGEN
    NOTE: Spannungen werden quasi 2D berechnet wenn für die entpsrechenden Lasten die resultierenden Angegebn sind passt das auch so 
    '''    
    ausgabe_an_knoten = 'FE'# 'Ebenen'#
    r_1, r_2, r_3 = 1,2,3 # Rundung: ziffern zahl
    for querschnitt in querschnitte:
        QS_label = querschnitt.name
        nabenhöhe = querschnitt.nabenhöhe
        t = round(querschnitt.wand_stärke * utils.unit_conversion(querschnitt.einheiten['Länge'], 'cm')) # cm
        tX = querschnitt.t_laengslagen * utils.unit_conversion(querschnitt.einheiten['Länge'], 'cm')
        tY = querschnitt.t_querlagen * utils.unit_conversion(querschnitt.einheiten['Länge'], 'cm')       

        schnittgrößen = utils.parse_schnittgrößen_labels(schnittgrößen_design)

        # TODO: Achtung: Reihenfolge der hier aufgerufenen Funktionen ist wichtig
        querschnitt.calculate_normalspannung(schnittgrößen, add_vorspannkraft_grob = False, ist_bauzustand = True)
        sigma_bauzustand = {'sigma_max':querschnitt.sigma_zug,'sigma_min':querschnitt.sigma_zug,
                            'nz_max':querschnitt.nz_max, 'nz_min':querschnitt.nz_min}

        if lastfall == 'max_druck':
            querschnitt.spannkraft_berechnung(schnittgrößen, Spannglieder.suspa_draht_ex['Stahlparameter'], verluste_pauschal = spannkraft_verlust_pauschal,  unit = 'MN') # NOTE muss bisher vor ausnutzung bestimmt werden 
            n_ext_erf, n_int_pro_segment, n_summe_int, P_ist_fuge = querschnitt.get_spannglied_staffelung(n_oberste_segmente= vorpsannungs_params['n_segmente_ext'], 
                                                                                            Pd_ext = vorpsannungs_params['Pd_ext'], 
                                                                                            Pd_int = vorpsannungs_params['Pd_int'],
                                                                                            last = lastfall)

            #if last == 'max_druck': #TODO das ist nur ein sehr hässlicher weg das schnell zu  machen -> für schub den ruck aus max druck ansetzen
            P_ist_fuge_max_druck = P_ist_fuge
            querschnitt.berechnungs_hinweise.append('   - Vorspannung nur im Lastfall für maximalen Zug/Druck bestimmt')

        querschnitt.calculate_ausnutzung_normalspannung(schnittgrößen, add_vorspannkraft_grob = True, plot_spannungsverlauf=False)
        querschnitt.reibungs_nachweis_horizontalfugen(schnittgrößen, nachweis_parameter['mu'], P = P_ist_fuge_max_druck)
        #TODO das ist nicht gut hier da scheibenschub berechnung zwingend nach allen anderen gemacht werden muss da Datenstruktur der Querschnittswerte verändert wird
        # -> vermutlich behoben 
        querschnitt.calculate_ausnutzung_scheibenschub(schnittgrößen, lamellenbreite=0.13)
        querschnitt.compute_effektive_festigkeiten_design('kurz')

        # Für die Ausgabe (NOTE ist etwas verwirrend, da die Ausnutzen der Drcukspannung die Vorspannung intern beachtet und nicht als externe Schnittgröße)
        for knoten_typ in knoten_typen:
            for dauer in ['egal','ständig']:
                schnittgrößen_design[knoten_typ][QS_label][nabenhöhe][dauer]['G'] = np.copy(schnittgrößen_design[knoten_typ][QS_label][nabenhöhe][dauer]['x'])
                schnittgrößen_design[knoten_typ][QS_label][nabenhöhe][dauer]['x'] -= P_ist_fuge[knoten_typ]     

        # NOTE Syntax: results_df.loc[:, (Level1, Level2, Level3)] = daten
        results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Höhe [' + querschnitt.einheiten['Länge'] + ']')] = np.around(querschnitt.section_absolute_heights[ausgabe_an_knoten],r_2)
        results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'd_achse [' + querschnitt.einheiten['Länge'] + ']')] = np.around(querschnitt.d_achse[ausgabe_an_knoten],r_2)
        if include_sgr:
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Mz [MNm]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['g']* utils.unit_conversion('Nm', 'MNm') ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'N [MN]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['x']* utils.unit_conversion('N', 'MN') ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'G [MN]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['G']* utils.unit_conversion('N', 'MN') ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Q [MN]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['y']* utils.unit_conversion('N', 'MN') ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Mx [MNm]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['a']* utils.unit_conversion('Nm', 'MNm') ,r_2)

            results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_max [' + querschnitt.einheiten['Normalspannung'] + ']')] = np.around(querschnitt.sigma_zug_design[ausgabe_an_knoten],r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_min [' + querschnitt.einheiten['Normalspannung'] + ']')] = np.around(querschnitt.sigma_druck_design[ausgabe_an_knoten],r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. druck')] = np.around(querschnitt.ausnutzung_druck[ausgabe_an_knoten],r_3)
            #max_ausnutzung[QS_label][str(nabenhöhe) + ' m'] = max(querschnitt.ausnutzung_druck[ausgabe_an_knoten])

            #results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'P_erf [MN]')] = np.around(querschnitt.P_erf * utils.unit_conversion('N', 'MN') ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'P_erf,end [MN]')] = np.around(querschnitt.P_m0[ausgabe_an_knoten] * utils.unit_conversion('N', 'MN') ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'P_ist [MN]')] = np.around(querschnitt.P_ist_fuge[ausgabe_an_knoten] * utils.unit_conversion('N', 'MN') ,r_2)

        if include_tau:
            if werkstoff_parameter['Furnierebene']:
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Fe' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_Fe_design[ausgabe_an_knoten],r_2)    
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_längs' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_längs_design[ausgabe_an_knoten],r_2)
                    
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub Furnier')] = np.around(querschnitt.ausnutzung_schub_Fe[ausgabe_an_knoten],r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub längs')] = np.around(querschnitt.ausnutzung_schub_längs[ausgabe_an_knoten],r_3)

                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_Fe')] = np.around(querschnitt.fvFed,r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_längs')] = np.around(querschnitt.fvd_brutto,r_3)
            else:
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Qy' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_Qy_design[ausgabe_an_knoten],r_2)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Mx' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_Mx_design[ausgabe_an_knoten],r_2)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_xy' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_xy_design[ausgabe_an_knoten],r_2)    
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_tor' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_vtor_design[ausgabe_an_knoten],r_2)
                    
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub brutto')] = np.around(querschnitt.ausnutzung_schub_brutto[ausgabe_an_knoten],r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub netto')] = np.around(querschnitt.ausnutzung_schub_netto[ausgabe_an_knoten],r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub torsion')] = np.around(querschnitt.ausnutzung_schub_torsion[ausgabe_an_knoten],r_3)

                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_brutto')] = np.around(querschnitt.fvd_brutto,r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_netto')] = np.around(querschnitt.fvxyd_netto,r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_tor')] = np.around(querschnitt.fvtord,r_3)

        if include_reibung:
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'nxy_Qy' + ' [kN/m]')] = np.around(querschnitt.nxy_Qy[ausgabe_an_knoten] * utils.unit_conversion(querschnitt.einheiten['Grundschnittgrößen_f'], 'kN/m'),r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'nxy_Mx' + ' [kNm/m]')] = np.around(querschnitt.nxy_Mx[ausgabe_an_knoten]* utils.unit_conversion(querschnitt.einheiten['Grundschnittgrößen_m'], 'kNm/m'),r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'nxy_P,Rd' + ' [kN/m]')] = np.around(querschnitt.n_Rd[ausgabe_an_knoten]* utils.unit_conversion(querschnitt.einheiten['Grundschnittgrößen_f'], 'kN/m'),r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. reibung')] = np.around(querschnitt.ausnutzung_reibung[ausgabe_an_knoten],r_3)

        if include_sigma_M_N:
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_N [' + querschnitt.einheiten['Normalspannung'] + ']')] = querschnitt.sigma_N_design[ausgabe_an_knoten]
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_M [' + querschnitt.einheiten['Normalspannung'] + ']')] = querschnitt.sigma_M_design[ausgabe_an_knoten]

        # maximalen Ausnutzungen sammeln
        max_ausnutzung = utils.collect_max_ausn(results_df, max_ausnutzung, lastfall)

    if lastfall == 'max_druck':
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'Höhe [m]')] = np.around(querschnitt.section_absolute_heights['Ebenen'],r_1)
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'd_achse [m]')] = np.around(querschnitt.d_achse['Ebenen'],r_1)
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'P_erf,end [MN]')] = np.around(querschnitt.P_m0['Ebenen']* utils.unit_conversion('N', 'MN') ,r_1)
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n EX-'+str(n_draht_suspa_ex))] = n_ext_erf
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n Mono-'+str(n_mono_litzen))] = n_int_pro_segment
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n int. Summe')] = n_summe_int
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'P_ist Fuge [MN]')] = np.around(P_ist_fuge['Ebenen'] * utils.unit_conversion('N', 'MN') ,r_1)

# TODO der hier müsste mal dahingehend gechekct werden welche einfluss die Kopflast hat
for segment in schnittgrößen_design['bauzustand'][QS_label][nabenhöhe]:
    if include_bauzustand_sgr:
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Mz [MNm]')] = np.around(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['g']* utils.unit_conversion('Nm', 'MNm') ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'N [MN]')] = np.around(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['x']* utils.unit_conversion('N', 'MN') ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Q [MN]')] = np.around(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['y']* utils.unit_conversion('N', 'MN') ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Mx [MNm]')] = np.around(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['y']* utils.unit_conversion('Nm', 'MNm') ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Reibung [MN]')] = np.around(bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung'] * utils.unit_conversion('N', 'MN') ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Ausn. reibung')] = np.around(bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung Ausn.'],r_2)
    if include_bauzustand_spannung:
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, GD.GREEK_UNICODE['sigma'] + '_max [' + querschnitt.einheiten['Normalspannung'] + ']')] = np.around(sigma_bauzustand['sigma_max']['Ebenen'][segment] * querschnitt.einheiten_umrechnung['Normalspannung'], r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, GD.GREEK_UNICODE['sigma'] + '_min [' + querschnitt.einheiten['Normalspannung'] + ']')] = np.around(sigma_bauzustand['sigma_min']['Ebenen'][segment]* querschnitt.einheiten_umrechnung['Normalspannung'], r_2)
    if include_bauzustand_sgr_m:
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'nz_max [kN/m]')] = np.around(sigma_bauzustand['nz_max']['Ebenen'][segment] * utils.unit_conversion('N/m', 'kN/m') ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'nz_min [kN/m]')] = np.around(sigma_bauzustand['nz_min']['Ebenen'][segment] * utils.unit_conversion('N/m', 'kN/m') ,r_2)

print ('\nAnmerkungen zu Annahmen bei der Berechnung:')
for anmerkung in list(set(querschnitt.berechnungs_hinweise)):
    print (anmerkung)
anmerkungs_df = pd.DataFrame({'Anmerkungen zu Annahmen bei der Berechnung':list(set(querschnitt.berechnungs_hinweise))})

max_results_df = pd.DataFrame(max_ausnutzung)
print ('\nMaximale Ausnutzungen', QS_label)
print (max_results_df)

#______________________ ALLES IN EXCEL SPEICHERN _________________________________________________________________________________________

with pd.ExcelWriter(results_excel, mode= 'w', engine="openpyxl") as writer:# --> mode = 'a', if_sheet_exists='overlay'
    nrows = nachweis_parameter_df.shape[0]
    nz_max_cols = utils.get_spalten_nach_name(bauzustands_df, df_bauzustand_header_list, start_row = 0, name = 'nz_max')
    reibung_cols = utils.get_spalten_nach_name(bauzustands_df, df_bauzustand_header_list, start_row = 0, name = 'Ausn. reibung')

    holz_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=0, startcol=0, index=False)#
    holz_units_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=3, startcol=0, index=False)#
    qs_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=7, startcol=0)
    lagen_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=qs_df.shape[0] + 13, startcol=0)

    for lastfall in lastfälle:
        einwirkungs_parameter_df[lastfall].to_excel(writer, sheet_name= 'Ed_' + lastfall, startrow=0, startcol=0, index=False)
        einwirkungs_df[lastfall].to_excel(writer, sheet_name= 'Ed_' + lastfall, startrow=4, startcol=0)
        einwirkungs_df_typ[lastfall].to_excel(writer, sheet_name= 'Ed_' + lastfall, startrow=einwirkungs_df[lastfall].shape[0]+12, startcol=0)

    anmerkungs_df.to_excel(writer, sheet_name= 'Hinweise_Berechnung', startrow=1, startcol=0, index=False)    

    vorpsannungs_df.to_excel(writer, sheet_name= 'Vorspannung_Staffelung', startrow=5, startcol=0)
    spannglieder_werte_df.to_excel(writer, sheet_name= 'Vorspannung_Staffelung', startrow=0, startcol=0, index = False)

    for lastfall in lastfälle:
        results_df[lastfall].to_excel(writer, sheet_name= 'Ergebnisse_' +lastfall, startrow=nrows+2, startcol=0, index=True)#, float_format='%.4f')
        nachweis_parameter_df.to_excel(writer, sheet_name= 'Ergebnisse_' +lastfall, startrow=0, startcol=0, index=False)
        ausnutzungs_cols = utils.get_spalten_ausnutzung(results_df[lastfall], df_results_header_list, start_row = nrows+2)
    max_results_df.to_excel(writer, sheet_name= 'Ausnutzungen_max', startrow=1, index=False)

    bauzustands_df.to_excel(writer, sheet_name= 'Bauzustand', startcol=0, index=True)


utils.zellen_groeße_formatieren(results_excel, worksheet= 'Vorspannung_Staffelung', cell_width=16, n_cols=len(spannglieder_werte_df.columns)+1)
utils.zellen_groeße_formatieren(results_excel, worksheet= 'QS_Werte', cell_width=15, n_cols=max(len(qs_df.columns), len(holz_df.columns))+1)
for lastfall in lastfälle:
    utils.zellen_groeße_formatieren(results_excel, worksheet= 'Ed_'+lastfall, cell_width=15, n_cols=len(einwirkungs_df[lastfall].columns)+1)
    utils.zellen_groeße_formatieren(results_excel, worksheet= 'Ergebnisse_'+lastfall, cell_width=18, n_cols=len(results_df[lastfall].columns)+1)
    utils.add_databar_color(results_excel, worksheet = 'Ergebnisse_'+lastfall, columns = ausnutzungs_cols)

utils.zellen_groeße_formatieren(results_excel, worksheet= 'Ausnutzungen_max', cell_width=18, n_cols=len(max_results_df.columns)+1)
utils.zellen_groeße_formatieren(results_excel, worksheet= 'Hinweise_Berechnung', cell_width=100, n_cols=1, start_col='A')
utils.zellen_groeße_formatieren(results_excel, worksheet= 'Bauzustand', cell_width=15, n_cols=len(bauzustands_df.columns)+1, start_col='A')

utils.add_color_rule(results_excel, worksheet = 'Bauzustand', columns = nz_max_cols, regel = 'greaterThan', value = 0)
utils.add_color_rule(results_excel, worksheet = 'Bauzustand', columns = reibung_cols, regel = 'greaterThan', value = 1)

#utils.zelle_beschriften(results_excel, 'Ausnutzungen_max', 'B' + str(start_row_max_ti[t]), 
#                       't, tX, tY [cm] ' + ', '.join([str(int(t)), str(round(tX,1)), str(round(tY,1))]) ,'B' + str(start_row_max_ti[t])+ ':E'+ str(start_row_max_ti[t]))
# TODO das hier auf max_druck schub unterscheidung anpassen -> QS Werte werden bei der Ermüdungsberechnung gebraucht
dataframes_doc = {'QS_Werte_holz':holz_df,'QS_Werte_geometrie':qs_df,'QS_Werte_lagen':lagen_df,
                  'Einwirkungs_parameter':einwirkungs_parameter_df,'Einwirkung_nach_dauer':einwirkungs_df,'Einwirkung_nach_typ':einwirkungs_df_typ,
                  'Vorspannung':vorpsannungs_df,
                  'Anmerkungen':anmerkungs_df, 'Nachweis_parameter':nachweis_parameter_df,
                  'SGR_Spannungen':results_df,}
#utils.save_dataframes_to_pkl(dataframes_doc, ausgabe_an_knoten)

from datetime import datetime
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M")
print ('Ergebnisse in ', results_excel, 'geschrieben -', dt_string)

    

