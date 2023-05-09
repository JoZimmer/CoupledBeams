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
import json

from source.model import BeamModel
import source.postprocess as postprocess
import source.utilities.utilities as utils
from source.Querschnitte import KreisRing
import source.utilities.global_definitions as GD
from inputs import holz
from inputs.spannglieder import Spannglieder
import inputs.parameter_H130 as parameter_input

# with open('turm_parameter.json', 'r') as parameter_file:
#     parameters = json.loads(parameter_file.read())
# json Format unpraktisch weil man keine Kommentare hinzufügen kann
parameters = parameter_input.params_dict

params_beam = utils.add_defaults(parameters["model_parameter"])
params_material = holz.charakteristische_werte[parameters["material"]["holzguete"]]
params_lagen = parameters['lagenaufbau']
params_geometrie = parameters['turm_geometrie']
params_vorspannung = parameters['vorspannung']
params_nachweise = parameters['nachweise']
params_lasten = parameters['lasten']
params_output = parameters['output']
digits = params_output['Nachkommastellen']

# zusätzliches zur vorspannung was erst noch gerechnet werden muss
n_draht_suspa_ex = params_vorspannung['n_draehte_suspa_ex'] 
n_mono_litzen = params_vorspannung['n_monolitzen']  # Nabenhöhe > 100: 5
params_vorspannung.update({
                        'Pd_ext':Spannglieder.suspa_draht_ex[params_vorspannung['stahl_suspa_ex']]['Pm0_n'][n_draht_suspa_ex], 
                        'Pd_int':Spannglieder.monolitzen[params_vorspannung['stahl_monolitzen']]['Pmax_n'][n_mono_litzen]})

# wählen welche gerechnet werden sollen
lastfälle = params_lasten["kopflasten"]['lastfälle']

# sind schon designlasten Cosy noch IEA
kopf_lasten_IEA = { '@max_Fxy':{'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':5.98E+06, 'My':6.81E+06, 'Mz':2.31E+06},# NOTE My = Mxy #Masse charakt. 2.62E+06, Fx = Fxy 
                    '@max_Mz':{'Fx':419925,'Fy':-24468, 'Fz':-3650438, 'Mx':6294993, 'My':-3439939, 'Mz':10811885}, # Fx = Fxy, My = Mxy
                    '@max_all':{'Fx':1.17E+06,'Fy':4.80E+04, 'Fz':-3.64E+06, 'Mx':5.98E+06, 'My':6.81E+06, 'Mz':10811885}} # bisher max Fx und max Mz 

params_beam['material_density'] = params_material['rhok']
params_beam['E_Modul'] = params_material['E0mean']

n_nodes = params_beam['n_elements']+1

'''
------------------------------------------------------------------------------------------------------
QUERSCHNITTS DEFINITION
------------------------------------------------------------------------------------------------------
'''
knoten_typen = ['Ebenen','FE']
querschnitts_dateien_laden = False
if not querschnitts_dateien_laden:
    querschnitte,qs_werte, qs_werte_FE, qs_labels  = [],[],[],[]

    label = 'Gerade'
    for t in params_lagen['dicken']:

        lagenaufbauten = holz.Lagenaufbauten(params_lagen)
        lagen_aufbau = lagenaufbauten.get()
        
        kreis_ring = KreisRing(cd = params_lasten['turmlasten']['cd_zylinder'], cp_max=1.8, lagen_aufbau=lagen_aufbau,
                                holz_parameter = params_material, 
                                nachweis_parameter = params_nachweise,
                                hoehen_parameter= params_geometrie, einheiten=parameters['einheiten_input'],
                                FE_elements=params_beam['n_elements'])

        kreis_ring.name += ' ' + label + ' t' + str(t)
        qs_label_openfast = kreis_ring.name + '_h' + str(params_geometrie['nabenhöhe']) + '_du' + str(params_geometrie['d_unten_oben'][0]) +\
                        '_do' + str(params_geometrie['d_unten_oben'][1])
        qs_labels.append(kreis_ring.name)
        for values in kreis_ring.querschnitts_werte['Ebenen'].values():
            qs_werte.append(np.around(values, digits['QS_werte']))
        for values in kreis_ring.querschnitts_werte['FE'].values():
            qs_werte_FE.append(np.around(values, digits['QS_werte']))

        querschnitte.append(kreis_ring)
   
    qs_header = pd.MultiIndex.from_product([qs_labels, list(kreis_ring.querschnitts_werte['Ebenen'].keys())], names=['', 'Ebene'])
    qs_header_FE = pd.MultiIndex.from_product([qs_labels, list(kreis_ring.querschnitts_werte['FE'].keys())], names=['', 'Knoten'])
    #lagen_header = pd.MultiIndex.from_product([['Lagenaufbau [m]'], list(lagen_aufbau[0].keys())])

    qs_werte=np.array(qs_werte).transpose()
    qs_werte_FE=np.array(qs_werte_FE).transpose()
    qs_df = pd.DataFrame(qs_werte, columns=qs_header)
    qs_FE_df = pd.DataFrame(qs_werte_FE, columns=qs_header_FE)
    holz_df = pd.DataFrame(params_material, index=[0])
    holz_units_df = pd.DataFrame(holz.charakteristische_werte['units'])
    lagen_df = pd.DataFrame(lagen_aufbau,  index=list(range(1,len(lagen_aufbau)+1)))
    lagen_df.rename({'ti': 'ti [' + kreis_ring.einheiten['Länge'] +']', 'di': 'di [' + kreis_ring.einheiten['Länge']+']'}, axis='columns',inplace=True)
    lagen_df.index.name = 'Lage'
    
    # Querschnittswerte speichern für openfast einlesen -> NOTE gehe hier davon aus das gerade nur ein QS
    qs_FE_df.to_pickle(os_join(*['output','Turmwerte_OpenFAST', qs_label_openfast + '.pkl']))
    print ('  Querschnittswerte in dictionary als pickle gespeichert in:', os_join(*['output','Turmwerte_O_spenFAST', qs_label_openfast + '.pkl']))

'''
------------------------------------------------------------------------------------------------------
DATAFRAME
------------------------------------------------------------------------------------------------------
'''
# https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy

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
    df_results_header_list[2].append('My [MNm]') 
    df_results_header_list[2].append('N [MN]') 
    df_results_header_list[2].append('G [MN]') 
    df_results_header_list[2].append('Q [MN]') 
    df_results_header_list[2].append('Mz [MNm]') #Torsion
if include_sigma_M_N:
    df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_N [' + qs.einheiten['Normalspannung'] + ']') 
    df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_M [' + qs.einheiten['Normalspannung'] + ']')  
df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_max [' + qs.einheiten['Normalspannung'] + ']')  
df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_min [' + qs.einheiten['Normalspannung'] + ']')  
df_results_header_list[2].append('Ausn. druck') 

if include_tau:
    if params_material['Furnierebene']:
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
    df_results_header_list[2].append('nxy_Mz' + ' [kNm/m]') 
    df_results_header_list[2].append('nxy_P,Rd' + ' [kN/m]') 
    df_results_header_list[2].append('Ausn. reibung') 

for segment in range(1,qs.n_ebenen):
    seg = 'seg_0-' + str(segment)
    df_bauzustand_header_list[2].append(seg)
if include_bauzustand_sgr:
    df_bauzustand_header_list[3].extend(['My [MNm]', 'N [MN]', 'Q [MN]','Mz [MNm]', 'Reibung [MN]', 'Ausn. reibung'])
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
results_df = {last:pd.DataFrame(columns=df_results_header) for last in lastfälle}#pd.DataFrame(columns=df_results_header)#
einwirkungs_parameter_df = {}
einwirkungs_df = {last:pd.DataFrame(columns=df_einwirkungs_header) for last in lastfälle}
einwirkungs_df_typ = {last:pd.DataFrame(columns=df_einwirkungs_header_typ) for last in lastfälle}

# Unabhängig von der Belastung
if include_bauzustand_sgr:
    bauzustands_df = pd.DataFrame(columns=df_bauzustand_header)
vorpsannungs_df = pd.DataFrame(columns=df_vorspannung_header)

spannglieder_werte_df = pd.DataFrame(Spannglieder.parameter_dict(params_vorspannung))

'''
------------------------------------------------------------------------------------------------------
ÄUßERE BELASTUNG
------------------------------------------------------------------------------------------------------
'''

qs.berechnungs_hinweise.append('   - Kopflast Fx (beam cosy neu) entspricht der resultierenden Fxy aus den IEA Lasten')
qs.berechnungs_hinweise.append('   - Kopflast My (beam cosy neu) entspricht der resultierenden Mxy aus den IEA Lasten')

sk = params_lasten['kopflasten']['skalierungs_rotor_flaechen']
skalierung = sk[0] / sk[1] * sk[2]
if skalierung != 1.0:
    kopf_lasten_IEA = utils.scale_dict_of_dicts(kopf_lasten_IEA, skalierung)
    qs.berechnungs_hinweise.append('   - Kopflasten IEA skaliert anhand der Rotorflächen mit  (und dann x 2)' + str(round(skalierung,5)))

kopf_masse_design = - params_beam['nacelle_mass'] * GD.GRAVITY * params_nachweise['sicherheitsbeiwerte']['g'] # -2628197

import inputs.DIN_Windlasten as wind_DIN
terrain_kategorie = params_lasten['turmlasten']['terrain_kategorie']
basis_windgeschwindigkeit = wind_DIN.vb_von_v_nabenhöhe(params_lasten['turmlasten']['wind_nabenhoehe'], terrain_kategorie, qs.nabenhöhe) #17

vb_bauzustad = wind_DIN.VB_WINDZONEN[params_lasten['turmlasten']['windzone']] # Windzone 2 vb = 25 m/s TODO heir extra windlast berechne und für den bauzustand verwenden

params_nachweise_dict = { GD.GREEK_UNICODE['gamma']+'_M': [params_nachweise['gamma_M']],
                            'kmod Ew. kurz': [params_nachweise['k_mod']['kurz']],
                            'kmod Ew. ständig': [params_nachweise['k_mod']['lang']],
                            'ksys': [params_nachweise['k_sys']],
                            'Reibung; ' + GD.GREEK_UNICODE['mu']: [params_material['mu']]
}

params_nachweise_df = pd.DataFrame(params_nachweise_dict)

print ('Windbelastung aus:')
print ('    vb 10 m:', round(basis_windgeschwindigkeit,2))
print ('    terrain:', terrain_kategorie)
#wind_DIN.plot_DIN_all(basis_windgeschwindigkeit, categories=[terrain_kategorie])

max_ausnutzung = {}
for last_i, lastfall in enumerate(lastfälle):
    last_at = lastfälle[lastfall]
    print ('\nKopflast von IEA:\n  ', last_at)

    querschnitte[0].berechnungs_hinweise.append('   - Last am Kopf aus IEA ' + last_at)
    einwirkungs_parameter = {'vb':round(basis_windgeschwindigkeit,2),'Terrain Kategorie':terrain_kategorie, 
                            'cd': params_lasten['turmlasten']['cd_zylinder'], 
                            'Kopflast':last_at, 'DLC':'1.3_seed2_9ms', 'Kosy.':'Balken neu'}

    einwirkungs_parameter.update(params_nachweise['sicherheitsbeiwerte'])

    einwirkungs_parameter_df[lastfall] = pd.DataFrame(einwirkungs_parameter, index=[0])

    kopf_lasten_beam = utils.convert_coordinate_system_and_consider_einwirkungsdauer(kopf_lasten_IEA[last_at], n_nodes, params_nachweise['sicherheitsbeiwerte']['g'], 
                                                                                    kopf_masse = kopf_masse_design, dimension=params_beam['dimension']) 

    wind_kraft_z = {}
    lasten_files, lasten_dicts_dauer, lasten_dicts_typ = {}, {}, {}
    for qs in querschnitte:
        QS_label = qs.name
        lasten_dict_base = {fi:np.zeros(n_nodes) for fi in GD.FORCES[params_beam['dimension']]}

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
            knoten_wind_kraft_z[knoten_typ] = utils.linien_last_to_knoten_last_design(wind_kraft_z, qs.section_absolute_heights[knoten_typ], gamma_q = params_nachweise['sicherheitsbeiwerte']['wind']) # NOTE gamma_q nach IEC DLC 1.3
            #ebenen_wind_kraft_z = utils.linien_last_to_knoten_last(wind_kraft_z, qs.section_absolute_heights, gamma_q = params_nachweise['sicherheitsbeiwerte']['wind']) # NOTE gamma_q nach IEC 
            
        #postprocess.plot_along_height(knoten_wind_kraft_z, z_coords, label='wind_kraft [N]')

        gewichtskraft_design = {'FE':{k: params_nachweise['sicherheitsbeiwerte']['g']*v for k,v in qs.gewichtskraft['FE'].items()}} # im QS ist die charakteristische gespeichert
        gewichtskraft_design['Ebenen'] = {k: params_nachweise['sicherheitsbeiwerte']['g']*v for k,v in qs.gewichtskraft['Ebenen'].items()} 

        #_____Lasten sortiert nach Typ der Last
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['kopflast'] = kopf_lasten_beam['egal']
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['windkraft'] = knoten_wind_kraft_z['FE']
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['eigengewicht'] = gewichtskraft_design['FE']
        F_h_imp_ers = utils.horizontale_ersatzlast([gewichtskraft_design['FE'], kopf_lasten_beam['egal']], theta_x = params_beam['imperfektion'])
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['F_h,ers,imp'] = F_h_imp_ers
        lasten_dicts_typ[QS_label][qs.nabenhöhe]['bauzustand'] = utils.lasten_dict_bauzustand(lasten_dict_base, knoten_wind_kraft_z, qs.gewichtskraft['Ebenen'], params_nachweise['sicherheitsbeiwerte'],
                                                                                            QS_obj=qs, abminderung_wind=0.5, parameter = params_beam)

        # _____Sortieren nach params_nachweise['Einwirkungsdauer']
        lasten_dicts_dauer[QS_label][qs.nabenhöhe]['egal']  = utils.update_lasten_dict(lasten_dict_base, [knoten_wind_kraft_z['FE'], F_h_imp_ers, 
                                                                                                         gewichtskraft_design['FE'], kopf_lasten_beam['egal']]) 
        lasten_dicts_dauer[QS_label][qs.nabenhöhe]['kurz'] = utils.update_lasten_dict(lasten_dict_base, [knoten_wind_kraft_z['FE'], kopf_lasten_beam['kurz']])
        # TODO Ein Teil es Kopfmoments My (beam neu) ist auch ständig!
        lasten_dicts_dauer[QS_label][qs.nabenhöhe]['ständig']  = utils.update_lasten_dict(lasten_dict_base, [gewichtskraft_design['FE'], F_h_imp_ers, kopf_lasten_beam['ständig']])
        # das ist der Lastfall spannkraft -> ist für die Spannkraftberechnung (gamma_m Eigengewicht = 1,0)
        lasten_dicts_dauer[QS_label][qs.nabenhöhe]['spannkraft']  = utils.update_lasten_dict(lasten_dict_base, [knoten_wind_kraft_z['FE'], qs.gewichtskraft['FE'], F_h_imp_ers, kopf_lasten_beam['spannkraft']])
    
        # ____________ LASTEN DATEI GENERIEREN - nur nach dauer sortiert da dies relevant für ausnutzung ist ___________________________
        for dauer in params_nachweise['einwirkungsdauer']:
            filename = 'K-IEA'+ last_at + '_W-v' +str(round(basis_windgeschwindigkeit,1)) + 'cd' + str(qs.cd) + '_D-' + dauer
            lasten_files[QS_label][qs.nabenhöhe][dauer] = utils.generate_lasten_file(n_nodes, lasten_dicts_dauer[QS_label][qs.nabenhöhe][dauer], 
                                                                                     file_base_name=filename, dimension=params_beam['dimension'])
        
        # Händische Konstante Torsion dazu wird auch bei den schnittgrößen seperat betrachtet
        if params_beam['dimension'] == '2D':
            qs.berechnungs_hinweise.append('   - Torsion konstant über die Höhe nur aus Kopflast')
            lasten_dicts_dauer[QS_label][qs.nabenhöhe]['egal']['Mz']  = np.append(np.zeros(params_beam['n_elements']), kopf_lasten_IEA[last_at]['Mz'])
            lasten_dicts_dauer[QS_label][qs.nabenhöhe]['kurz']['Mz']  = np.append(np.zeros(params_beam['n_elements']), kopf_lasten_IEA[last_at]['Mz'])
            lasten_dicts_typ[QS_label][qs.nabenhöhe]['kopflast']['Mz'] = np.append(np.zeros(params_beam['n_elements']), kopf_lasten_IEA[last_at]['Mz'])

        einheiten_out = params_output['einheiten']['lasten']
        # LASTEN NACH DAUER SORTIERT
        for dauer in params_nachweise['einwirkungsdauer']:
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
                                                                                                    file_base_name=filename, dimension=params_beam['dimension'])
            
            else:
                filename = 'K-IEA'+ last_at + '_W-v' +str(round(basis_windgeschwindigkeit,1)) + 'cd' + str(qs.cd) + '_D-' + typ
                lasten_files[QS_label][qs.nabenhöhe][typ] = utils.generate_lasten_file(n_nodes, lasten_dicts_typ[QS_label][qs.nabenhöhe][typ], 
                                                                                        file_base_name=filename, dimension=params_beam['dimension'])
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
            for dauer in params_nachweise['einwirkungsdauer']:
                postprocess.plot_dict_subplots(lasten_dicts_dauer[QS_label][qs.nabenhöhe][dauer], qs.x_FE, title='Last ' + dauer +' ' + QS_label, unit = 'kN')


    '''
    ------------------------------------------------------------------------------------------------------
    SCHNITTGRÖßEN
    ------------------------------------------------------------------------------------------------------
    '''
    # Werden getrennt an den Stellen der Ebenen und den FE Knoten ausgewertet
    SGR_design = {'Ebenen':{},'FE':{},'bauzustand':{}}
    bauzustands_kram = {}

    for querschnitt in querschnitte:
        QS_label = querschnitt.name
        nabenhöhe = querschnitt.nabenhöhe

        if QS_label not in bauzustands_kram:
            bauzustands_kram[QS_label] = {}
        if nabenhöhe not in bauzustands_kram[QS_label]:
            bauzustands_kram[QS_label][nabenhöhe] = {}

        for knoten_typ in SGR_design:
            if QS_label not in SGR_design[knoten_typ]:
                SGR_design[knoten_typ][QS_label] = {}
            if nabenhöhe not in SGR_design[knoten_typ][QS_label]:
                SGR_design[knoten_typ][QS_label][nabenhöhe] = {}

        section_properties = querschnitt.section_parameters 

        parameters = {'FE': utils.add_model_data_from_dict(section_properties['FE'], params_beam)}
        parameters['Ebenen'] = utils.add_model_data_from_dict(section_properties['Ebenen'], params_beam)
        parameters['Ebenen']['n_elements'] = querschnitt.n_sections
        if params_beam['type_of_bc'] == 'Feder':
            querschnitt.berechnungs_hinweise.append('   - Einspannung am Fuß mit Federsteifigkeiten u & Drehfder gamma: ' +  str(params_beam['spring_stiffness']))
        elif params_beam['type_of_bc'] == 'Eingespannt':
            querschnitt.berechnungs_hinweise.append('   - Einspannung am Fuß starr')

        beam = BeamModel(parameters['FE'], adjust_mass_density_for_total = False, optimize_frequencies_init=False , apply_k_geo=False)
        beam_ebenen = BeamModel(parameters['Ebenen'], adjust_mass_density_for_total = False, optimize_frequencies_init=False , apply_k_geo=False)
        
        if last_i == 0:
            print ('\nN_nodes', querschnitt.n_ebenen)
            print (QS_label, nabenhöhe, 'm')
            print ('     Gesamt Volumen des Querschnitts [m³]:', round(beam.volume,2))
            print ('     Gesamt Gewichtskraft Turm am Fuß [MN]:',round(sum(beam.eigengewicht * GD.UNIT_SCALE['MN']),2))
            print ('     Frequenzen [Hz]:', [round(beam.eigenfrequencies[i],3) for i in range(3)])

        lasten_nach_dauer = {'ständig':{'torsion':False},
                             'egal':{'torsion':kopf_lasten_IEA[last_at]['Mz']},
                             'kurz':{'torsion':kopf_lasten_IEA[last_at]['Mz']},
                             'spannkraft':{'torsion':kopf_lasten_IEA[last_at]['Mz']}}

        for dauer in params_nachweise['einwirkungsdauer']:
            # NOTE Torsion ist nicht in der Steifigkeitsmatrix, deswegen kann diese brechnung nur händisch ergänzt werden
            lasten_file = lasten_files[QS_label][nabenhöhe][dauer]
            SGR_design['FE'][QS_label][nabenhöhe][dauer] = {}
            SGR_design['Ebenen'][QS_label][nabenhöhe][dauer] = {}

            SGR_design['FE'][QS_label][nabenhöhe][dauer] = beam.static_analysis_solve(load_vector_file=lasten_file, 
                                                                        constant_torsion = lasten_nach_dauer[dauer]['torsion'])

            # SGR werte auch an Ebenen knoten bekommen -> Linear interpolieren 
            for key in SGR_design['FE'][QS_label][nabenhöhe][dauer]:
                SGR_design['Ebenen'][QS_label][nabenhöhe][dauer][key] = utils.interpolate(z_soll = querschnitt.section_absolute_heights['Ebenen'], 
                                                                                                    z_ist=beam.nodal_coordinates['x0'], 
                                                                                                    werte = SGR_design['FE'][QS_label][nabenhöhe][dauer][key])

            lasten_label = os.path.splitext(os.path.basename(lasten_files[QS_label][nabenhöhe][dauer]))[0]                                                                                        
        if last_i == 0:          
            print ('     Maximales Moment [MNm]:', round(max(abs(SGR_design['Ebenen'][QS_label][nabenhöhe]['egal']['g'])) * utils.unit_conversion('Nm', 'MNm'),2))   
            print ('     Maximales Moment ohne Kopflast [MNm] ~=', round((max(abs(SGR_design['Ebenen'][QS_label][nabenhöhe]['egal']['g'])) - kopf_lasten_IEA[last_at]['Fx']*nabenhöhe) * utils.unit_conversion('Nm', 'MNm'),2))   

        # BAUZUSTAND 
        if include_bauzustand_sgr:
            for segment in lasten_files[QS_label][nabenhöhe]['bauzustand']:
                SGR_design['bauzustand'][QS_label][nabenhöhe][segment]= {}
                bauzustands_kram[QS_label][nabenhöhe][segment] = {}
                lasten_file = lasten_files[QS_label][nabenhöhe]['bauzustand'][segment]
                SGR_design['bauzustand'][QS_label][nabenhöhe][segment] = beam_ebenen.static_analysis_solve(load_vector_file=lasten_file, return_result=True)
                bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung']  = -1 * utils.remove_small_values_from_array(SGR_design['bauzustand'][QS_label][nabenhöhe][segment]['x']) * params_material['mu']
                bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung Ausn.'] = np.divide(SGR_design['bauzustand'][QS_label][nabenhöhe][segment]['y'],
                                                                                            bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung'],
                                                                                            out=np.zeros_like(SGR_design['bauzustand'][QS_label][nabenhöhe][segment]['y']), 
                                                                                            where=bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung']!=0)

        plot_sgr = False                                                            
        if plot_sgr:
            lasten_label = os.path.splitext(os.path.basename(lasten_files[QS_label][nabenhöhe]['egal']))[0]
            postprocess.plot_static_result_forces(beam, 'internal_forces', ['x','y','g','a'], unit='MN', title_suffix='Last:' + lasten_label, figsize_scale = 1.1)
            postprocess.plot_static_result(beam, 'deformation', ['y'], unit='cm',title_suffix='Last:' + lasten_label,)
            #postprocess.plot_static_result(beam, 'deformation', ['g'], unit='cm',title_suffix='Last:' + lasten_label,)

    '''
    ------------------------------------------------------------------------------------------------------
    SPANNUNGEN UND AUSNUTZUNGEN
    NOTE: Spannungen werden quasi 2D berechnet wenn für die entpsrechenden Lasten die resultierenden Angegebn sind passt das auch so 
    ------------------------------------------------------------------------------------------------------
    '''    
    ausgabe_an_knoten = params_output['ausgabe_an_knoten']
    e_o_s = params_output["einheiten"]['schnittgrößen']
    e_o = params_output["einheiten"] # alles wo man nicht zwischen moment und kraft unterschieden muss
    r_1, r_2, r_3 = digits['Maße'], digits['default'], digits['Ausnutzung'] # Rundung: ziffern zahl
    for qs in querschnitte:
        QS_label = qs.name
        nabenhöhe = qs.nabenhöhe
        t = round(qs.wand_stärke * utils.unit_conversion(qs.einheiten['Länge'], 'cm')) # cm
        tX = qs.t_laengslagen * utils.unit_conversion(qs.einheiten['Länge'], 'cm')
        tY = qs.t_querlagen * utils.unit_conversion(qs.einheiten['Länge'], 'cm')       

        schnittgrößen = utils.parse_schnittgrößen_labels(SGR_design)

        # TODO: Achtung: Reihenfolge der hier aufgerufenen Funktionen ist wichtig
        qs.calculate_normalspannung(schnittgrößen, add_vorspannkraft_grob = False, ist_bauzustand = True)
        sigma_bauzustand = {'sigma_max':qs.sigma_zug,'sigma_min':qs.sigma_zug,
                            'nz_max':qs.nz_max, 'nz_min':qs.nz_min}

        # berechnung der Vporspannung nim Lastfall max_druck (sollte gleich max zug sein)
        if next(iter(lastfälle)) != 'max_druck':
            raise ('Achtung im Lastfall dict muss der erste Lastfall der max_druck Fall sein um die Spannkraft zu berechnen')

        if lastfall == 'max_druck':
            qs.spannkraft_berechnung(schnittgrößen, params = params_vorspannung) # NOTE muss vor der staffelung aufgerufen werden
            qs.get_spannglied_staffelung(  n_oberste_segmente= params_vorspannung['n_segmente_extern'], 
                                                    Pd_ext = params_vorspannung['Pd_ext'], 
                                                    Pd_int = params_vorspannung['Pd_int'],
                                                    spannkanal = params_vorspannung['masse_spannkanal'])
            for knoten_typ in knoten_typen:
                qs.Iy[knoten_typ] -= qs.Iy_spannkanal[knoten_typ]


            #TODO das ist nur ein sehr hässlicher weg das schnell zu  machen -> für schub den ruck aus max druck ansetzen
            P_ist_fuge_max_druck = qs.P_ist_fuge
            qs.berechnungs_hinweise.append('   - Vorspannung nur im Lastfall für maximalen Zug/Druck bestimmt')

        qs.calculate_ausnutzung_normalspannung(schnittgrößen, add_vorspannkraft_grob = True, plot_spannungsverlauf=False)
        qs.reibungs_nachweis_horizontalfugen(schnittgrößen, params_material['mu'], P = P_ist_fuge_max_druck)
        #TODO das ist nicht gut hier da scheibenschub berechnung zwingend nach allen anderen gemacht werden muss da Datenstruktur der Querschnittswerte verändert wird
        # -> vermutlich behoben 
        qs.calculate_ausnutzung_scheibenschub(schnittgrößen, lamellenbreite=0.13)
        qs.compute_effektive_festigkeiten_design('kurz')

        # Für die Ausgabe (NOTE ist etwas verwirrend, da die Ausnutzen der Drcukspannung die Vorspannung intern beachtet und nicht als externe Schnittgröße)
        for knoten_typ in knoten_typen:
            for dauer in ['egal','ständig']:
                SGR_design[knoten_typ][QS_label][nabenhöhe][dauer]['G'] = np.copy(SGR_design[knoten_typ][QS_label][nabenhöhe][dauer]['x'])
                SGR_design[knoten_typ][QS_label][nabenhöhe][dauer]['x'] -= qs.P_ist_fuge[knoten_typ]     

        # NOTE Syntax: results_df.loc[:, (Level1, Level2, Level3)] = daten
        results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Höhe [' + qs.einheiten['Länge'] + ']')] = np.around(qs.section_absolute_heights[ausgabe_an_knoten],r_2)
        results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'd_achse [' + qs.einheiten['Länge'] + ']')] = np.around(qs.d_achse[ausgabe_an_knoten],r_2)
        if include_sgr:
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'My [MNm]')] = np.around(SGR_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['g']* utils.unit_conversion(qs.einheiten['Moment'], e_o_s['Moment']) ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'N [MN]')] = np.around(SGR_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['x']* utils.unit_conversion(qs.einheiten['Kraft'], e_o_s['Kraft']) ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'G [MN]')] = np.around(SGR_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['G']* utils.unit_conversion(qs.einheiten['Kraft'], e_o_s['Kraft']) ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Q [MN]')] = np.around(SGR_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['y']* utils.unit_conversion(qs.einheiten['Kraft'], e_o_s['Kraft']) ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Mz [MNm]')] = np.around(SGR_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['a']* utils.unit_conversion(qs.einheiten['Moment'], e_o_s['Moment']) ,r_2)

            results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_max [' + qs.einheiten['Normalspannung'] + ']')] = np.around(qs.sigma_zug_design[ausgabe_an_knoten],r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_min [' + qs.einheiten['Normalspannung'] + ']')] = np.around(qs.sigma_druck_design[ausgabe_an_knoten],r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. druck')] = np.around(qs.ausnutzung_druck[ausgabe_an_knoten],r_3)
            #max_ausnutzung[QS_label][str(nabenhöhe) + ' m'] = max(qs.ausnutzung_druck[ausgabe_an_knoten])

            #results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'P_erf [MN]')] = np.around(qs.P_erf * utils.unit_conversion(qs.einheiten['Kraft'], 'MN') ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'P_erf,end [MN]')] = np.around(qs.P_m0[ausgabe_an_knoten] * utils.unit_conversion(qs.einheiten['Kraft'], e_o['Vorspannung']) ,r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'P_ist [MN]')] = np.around(qs.P_ist_fuge[ausgabe_an_knoten] * utils.unit_conversion(qs.einheiten['Kraft'], e_o['Vorspannung']) ,r_2)

        if include_tau:
            if params_material['Furnierebene']:
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Fe' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(qs.tau_Fe_design[ausgabe_an_knoten],r_2)    
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_längs' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(qs.tau_längs_design[ausgabe_an_knoten],r_2)
                    
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub Furnier')] = np.around(qs.ausnutzung_schub_Fe[ausgabe_an_knoten],r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub längs')] = np.around(qs.ausnutzung_schub_längs[ausgabe_an_knoten],r_3)

                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_Fe')] = np.around(qs.fvFed,r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_längs')] = np.around(qs.fvd_brutto,r_3)
            else:
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Qx' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(qs.tau_Qx_design[ausgabe_an_knoten],r_2)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Mz' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(qs.tau_Mz_design[ausgabe_an_knoten],r_2)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_xy' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(qs.tau_xy_design[ausgabe_an_knoten],r_2)    
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_tor' + ' [' + qs.einheiten['Schubspannung'] + ']')] = np.around(qs.tau_vtor_design[ausgabe_an_knoten],r_2)
                    
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub brutto')] = np.around(qs.ausnutzung_schub_brutto[ausgabe_an_knoten],r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub netto')] = np.around(qs.ausnutzung_schub_netto[ausgabe_an_knoten],r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. schub torsion')] = np.around(qs.ausnutzung_schub_torsion[ausgabe_an_knoten],r_3)

                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_brutto')] = np.around(qs.fvd_brutto,r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_netto')] = np.around(qs.fvxyd_netto,r_3)
                results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'fvd_tor')] = np.around(qs.fvtord,r_3)

        if include_reibung:
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'nxy_Qy' + ' [kN/m]')] = np.around(qs.nxy_Qx[ausgabe_an_knoten] * utils.unit_conversion(qs.einheiten['Kraft'], 'kN'),r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'nxy_Mz' + ' [kNm/m]')] = np.around(qs.nxy_Mz[ausgabe_an_knoten]* utils.unit_conversion(qs.einheiten['Moment'], 'kNm'),r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'nxy_P,Rd' + ' [kN/m]')] = np.around(qs.n_Rd[ausgabe_an_knoten]* utils.unit_conversion(qs.einheiten['Kraft'], 'kN'),r_2)
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, 'Ausn. reibung')] = np.around(qs.ausnutzung_reibung[ausgabe_an_knoten],r_3)

        if include_sigma_M_N:
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_N [' + qs.einheiten['Normalspannung'] + ']')] = qs.sigma_N_design[ausgabe_an_knoten]
            results_df[lastfall].loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_M [' + qs.einheiten['Normalspannung'] + ']')] = qs.sigma_M_design[ausgabe_an_knoten]

        # maximalen Ausnutzungen sammeln
        max_ausnutzung = utils.collect_max_ausn(results_df, max_ausnutzung, lastfall)

    if lastfall == 'max_druck':
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'Höhe [m]')] = np.around(qs.section_absolute_heights['Ebenen'],r_1)
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'd_achse [m]')] = np.around(qs.d_achse['Ebenen'],r_1)
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'P_erf,end [MN]')] = np.around(qs.P_m0['Ebenen']* utils.unit_conversion(qs.einheiten['Kraft'], e_o['Vorspannung']) ,r_1)
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n EX-'+str(n_draht_suspa_ex))] = qs.n_ext_erf_list
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n Mono-'+str(n_mono_litzen))] = qs.n_int_pro_segment
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n int. Summe')] = qs.n_int_summe
        vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'P_ist Fuge [MN]')] = np.around(qs.P_ist_fuge['Ebenen'] * utils.unit_conversion(qs.einheiten['Kraft'], e_o['Vorspannung']) ,r_1)

qs_df.loc[:,(QS_label, 'Iy_löcher [m^4]')] = np.around(qs.Iy['Ebenen'],4)

# TODO der hier müsste mal dahingehend gechekct werden welche einfluss die Kopflast hat
for segment in SGR_design['bauzustand'][QS_label][nabenhöhe]:
    if include_bauzustand_sgr:
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'My [MNm]')] = np.around(SGR_design['bauzustand'][QS_label][nabenhöhe][segment]['g']* utils.unit_conversion(qs.einheiten['Moment'], e_o_s['Moment']) ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'N [MN]')] = np.around(SGR_design['bauzustand'][QS_label][nabenhöhe][segment]['x']* utils.unit_conversion(qs.einheiten['Kraft'], e_o_s['Kraft']) ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Q [MN]')] = np.around(SGR_design['bauzustand'][QS_label][nabenhöhe][segment]['y']* utils.unit_conversion(qs.einheiten['Kraft'], e_o_s['Kraft']) ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Mz [MNm]')] = np.around(SGR_design['bauzustand'][QS_label][nabenhöhe][segment]['a']* utils.unit_conversion(qs.einheiten['Moment'], e_o_s['Moment']) ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Reibung [MN]')] = np.around(bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung'] * utils.unit_conversion(qs.einheiten['Kraft'], e_o_s['Kraft']) ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Ausn. reibung')] = np.around(bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung Ausn.'],r_2)
    if include_bauzustand_spannung:
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, GD.GREEK_UNICODE['sigma'] + '_max [' + qs.einheiten['Normalspannung'] + ']')] = np.around(sigma_bauzustand['sigma_max']['Ebenen'][segment] * qs.einheiten_umrechnung['Normalspannung'], r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, GD.GREEK_UNICODE['sigma'] + '_min [' + qs.einheiten['Normalspannung'] + ']')] = np.around(sigma_bauzustand['sigma_min']['Ebenen'][segment]* qs.einheiten_umrechnung['Normalspannung'], r_2)
    if include_bauzustand_sgr_m:
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'nz_max [kN/m]')] = np.around(sigma_bauzustand['nz_max']['Ebenen'][segment] * utils.unit_conversion(qs.einheiten['Kraft'], 'kN') ,r_2)
        bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'nz_min [kN/m]')] = np.around(sigma_bauzustand['nz_min']['Ebenen'][segment] * utils.unit_conversion(qs.einheiten['Kraft'], 'kN') ,r_2)

print ('\nAnmerkungen zu Annahmen bei der Berechnung:')
for anmerkung in list(set(qs.berechnungs_hinweise)):
    print (anmerkung)
anmerkungs_df = pd.DataFrame({'Anmerkungen zu Annahmen bei der Berechnung':list(set(qs.berechnungs_hinweise))})

max_results_df = pd.DataFrame(max_ausnutzung)
print ('\nMaximale Ausnutzungen', QS_label)
print (max_results_df)

'''
------------------------------------------------------------------------------------------------------
EXCEL SPEICHERN
------------------------------------------------------------------------------------------------------
'''
results_excel = os_join(*params_output["Excel_Datei"])

with pd.ExcelWriter(results_excel, mode= 'w', engine="openpyxl") as writer:# --> mode = 'a', if_sheet_exists='overlay'
    nrows = params_nachweise_df.shape[0]
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
        params_nachweise_df.to_excel(writer, sheet_name= 'Ergebnisse_' +lastfall, startrow=0, startcol=0, index=False)
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
                  'Anmerkungen':anmerkungs_df, 'params_nachweise':params_nachweise_df,
                  'SGR_Spannungen':results_df,}
#utils.save_dataframes_to_pkl(dataframes_doc, ausgabe_an_knoten)

from datetime import datetime
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M")
print ('Ergebnisse in ', results_excel, 'geschrieben -', dt_string)

    

