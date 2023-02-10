import pandas as pd
import numpy as np
import source.utilities.global_definitions as GD
import source.utilities.utilities as utils

def df_from_nested_dict(dictonary:dict, append_short_cols:bool = True) -> pd.DataFrame:
    '''
    soll aus mehrstufigem Dictonary ein multiindex dataframe erstellen -> muss rekursiv sein
    append_short_cols: in dem dictionary werden alle spalten auf die gleich länge gebracht mit None werden
    '''

    reform = {}
    for outerKey, innerDict in dictonary.items():
        for innerKey, values in innerDict.items():
            if isinstance(values, dict):
                for innerKey1, values1 in values.items():
                    if isinstance(values1, dict):
                        for innerKey2, values2 in values1.items():
                            if isinstance(values2, dict):
                                for innerKey3, values3 in values2.items():
                                    reform[(outerKey, innerKey, innerKey1, innerKey2, innerKey3)] = values3
                            else:
                                reform[(outerKey, innerKey, innerKey1, innerKey2)] = values2
                    else:
                        reform[(outerKey, innerKey, innerKey1)] = values1
            else:
                reform[(outerKey, innerKey)] = values

    if append_short_cols:
        max_len = max([len(x) for x in reform.values()])
        for key, vals in reform.items():
            if len(vals) < max_len:
                diff = int(max_len - len(vals))
                if isinstance(vals, np.ndarray):
                    reform[key] = np.append(vals, [None] * diff)
                elif isinstance(vals,list):
                    reform[key].extend([None]*diff)

    return pd.DataFrame(reform)

class Writer(object):
    '''
    Soll anhand einiger bools die dataframes in dem ausgabe daten gesamelt werden erstellen.
    Dann sollten ergebnisse gegebn werden und die entsprechenden dfs gefüllt werden 
    Dann sollten diese in die Excel geschrieben werden
    TODO: generische umwandlung von mehr stufigen dictionaries in Multiindex dataframes:
    https://stackoverflow.com/questions/24988131/nested-dictionary-to-multiindex-dataframe-where-dictionary-keys-are-column-label 
    '''
    def __init__(self, includes:dict, querschnitte:list) -> None:

        '''
        Einheiten könnten hier mal noch zentraler gemanagetd werden und ggf auch checks durchgeführt
        '''
        
        self.include_tau = includes['include_tau']
        self.include_sigma_M_N = includes['include_sigma_M_N']
        self.include_sgr = includes['include_sgr']
        self.include_reibung = includes['include_reibung']
        self.include_bauzustand_sgr = includes['include_bauzustand_sgr']
        self.include_bauzustand_spannung = includes['include_bauzustand_spannung']
        self.include_bauzustand_sgr_m = includes['include_bauzustand_sgr_m'] # pro meter
        self.mit_furnier = includes['mit_funier']

        self.querschnitte = querschnitte


    def initialize_dataframes(self):

        df_results_header_list = [[],[], []] # 1. Liste = 1. Level usw...
        kraft_komponenten = ['Fx [kN]', 'Fy [kN]', 'Mz [kNm]','Mx [kNm]']
        df_einwirkung_header_list = [[],[],['ständig','kurz','egal'], kraft_komponenten]
        df_einwirkung_typ_header_list = [[],[],['Kopflast [kN]', 'windkraft [kN]', 'eigengewicht [kN]', 'F_h,ers,imp [kN]' ],kraft_komponenten]
        df_vorspannung_header_list = [[],[],['Höhe [m]', 'd_achse [m]', 'P_erf,end [MN]', 'n EX-84', 'n Mono-5', 'n int. Summe', 'P_ist Fuge [MN]']]
        df_bauzustand_header_list = [[],[],[],[]]

        for qs in self.querschnitte:
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
        #df_bauzustand_header_list[2].append('Höhe [' + qs.einheiten['Länge'] + ']') 
        df_results_header_list[2].append('d_achse [' + qs.einheiten['Länge'] + ']') 
        #df_results_header_list[2].append('P_erf [MN]') 
        df_results_header_list[2].append('P_erf,end [MN]') 
        df_results_header_list[2].append('P_ist [MN]') 

        if self.include_sgr:
            df_results_header_list[2].append('M [MNm]') 
            df_results_header_list[2].append('N [MN]') 
            df_results_header_list[2].append('G [MN]') 
            df_results_header_list[2].append('Q [MN]') 
        #df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_zug [' + qs.einheiten['Normalspannung'] + ']')
        if self.include_sigma_M_N:
            df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_N [' + qs.einheiten['Normalspannung'] + ']') 
            df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_M [' + qs.einheiten['Normalspannung'] + ']')  
        df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_max [' + qs.einheiten['Normalspannung'] + ']')  
        df_results_header_list[2].append(GD.GREEK_UNICODE['sigma'] + '_min [' + qs.einheiten['Normalspannung'] + ']')  
        df_results_header_list[2].append('Ausn. druck') 

        if self.include_tau:
            if self.mit_furnier:
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

        if self.include_reibung:
            df_results_header_list[2].append('nxy_Qy' + ' [kN/m]') 
            df_results_header_list[2].append('nxy_Mx' + ' [kNm/m]') 
            df_results_header_list[2].append('nxy_P,Rd' + ' [kN/m]') 
            df_results_header_list[2].append('Ausn. reibung') 

        # Bauzustand
        for segment in range(1,qs.n_ebenen):
            seg = 'seg_0-' + str(segment)
            df_bauzustand_header_list[2].append(seg)
        if self.include_bauzustand_sgr:
            df_bauzustand_header_list[3].extend(['M [MNm]', 'N [MN]', 'Q [MN]', 'Reibung [MN]', 'Ausn. reibung'])
        if self.include_bauzustand_spannung:
            df_bauzustand_header_list[3].append(GD.GREEK_UNICODE['sigma'] + '_max [' + qs.einheiten['Normalspannung'] + ']') 
            df_bauzustand_header_list[3].append(GD.GREEK_UNICODE['sigma'] + '_min [' + qs.einheiten['Normalspannung'] + ']') 
        if self.include_bauzustand_sgr_m:
            df_bauzustand_header_list[3].append('nz_max [kN/m]')
            df_bauzustand_header_list[3].append('nz_min [kN/m]')

        # Header in pd Format
        df_results_header = pd.MultiIndex.from_product(df_results_header_list, names=['Nabenhöhe', 'Querschnitt', 'Segment'])
        df_bauzustand_header = pd.MultiIndex.from_product(df_bauzustand_header_list, names=['Nabenhöhe', 'Querschnitt', 'Segment', 'Komponente'])
        df_einwirkungs_header = pd.MultiIndex.from_product(df_einwirkung_header_list, names=['Nabenhöhe', 'Querschnitt', 'Ew. Dauer', 'Komponente', ])
        df_einwirkungs_header_typ = pd.MultiIndex.from_product(df_einwirkung_typ_header_list, names=['Nabenhöhe', 'Querschnitt', 'Lasttyp','Komponente', ])
        df_vorspannung_header = pd.MultiIndex.from_product(df_vorspannung_header_list, names=['Nabenhöhe', 'Querschnitt', 'Segment'])
        # Dataframes leer 
        self.results_df = pd.DataFrame(columns=df_results_header)
        self.bauzustands_df = pd.DataFrame(columns=df_bauzustand_header)
        self.einwirkungs_df = pd.DataFrame(columns=df_einwirkungs_header)
        self.einwirkungs_df_typ = pd.DataFrame(columns=df_einwirkungs_header_typ)
        self.vorpsannungs_df = pd.DataFrame(columns=df_vorspannung_header)

    def fill_dataframes(self, ausgabe_an_knoten:str, schnittgrößen_design:dict):
        '''
        ausgabe_an_knoten = 'Ebenen' oder 'FE'# 

        '''
        knoten_typen = ['Ebenen', 'FE']
        max_ausnutzung = {}
        r_1, r_2, r_3 = 1,2,3 # Rundung: ziffern zahl

        for querschnitt in self.querschnitte:
            QS_label = querschnitt.name
            nabenhöhe = querschnitt.nabenhöhe
            # Für die Ausgabe (NOTE ist etwas verwirrend, da die Ausnutzen der Drcukspannung die Vorspannung intern beachtet und nicht als externe Schnittgröße)
            for knoten_typ in knoten_typen:
                for dauer in ['egal','ständig']:
                    schnittgrößen_design[knoten_typ][QS_label][nabenhöhe][dauer]['G'] = np.copy(schnittgrößen_design[knoten_typ][QS_label][nabenhöhe][dauer]['x'])
                    schnittgrößen_design[knoten_typ][QS_label][nabenhöhe][dauer]['x'] -= querschnitt.P_ist_fuge[knoten_typ]     

            # NOTE Syntax: results_df.loc[:, (Level1, Level2, Level3)] = daten
            self.results_df.loc[:,(nabenhöhe, QS_label, 'Höhe [' + querschnitt.einheiten['Länge'] + ']')] = np.around(querschnitt.section_absolute_heights[ausgabe_an_knoten],r_2)
            self.results_df.loc[:,(nabenhöhe, QS_label, 'd_achse [' + querschnitt.einheiten['Länge'] + ']')] = np.around(querschnitt.d_achse[ausgabe_an_knoten],r_2)
            if self.include_sgr:
                self.results_df.loc[:,(nabenhöhe, QS_label, 'M [MNm]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['g']* utils.unit_conversion('Nm', 'MNm') ,r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, 'N [MN]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['x']* utils.unit_conversion('N', 'MN') ,r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, 'G [MN]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['G']* utils.unit_conversion('N', 'MN') ,r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, 'Q [MN]')] = np.around(schnittgrößen_design[ausgabe_an_knoten][QS_label][nabenhöhe]['egal']['y']* utils.unit_conversion('N', 'MN') ,r_2)

                self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_max [' + querschnitt.einheiten['Normalspannung'] + ']')] =\
                    np.around(querschnitt.sigma_zug_design[ausgabe_an_knoten],r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_min [' + querschnitt.einheiten['Normalspannung'] + ']')] =\
                    np.around(querschnitt.sigma_druck_design[ausgabe_an_knoten],r_2)

                self.results_df.loc[:,(nabenhöhe, QS_label, 'Ausn. druck')] = np.around(querschnitt.ausnutzung_druck[ausgabe_an_knoten],r_3)
                max_ausnutzung[QS_label][str(nabenhöhe) + ' m'] = max(querschnitt.ausnutzung_druck[ausgabe_an_knoten])

                #self.results_df.loc[:,(nabenhöhe, QS_label, 'P_erf [MN]')] = np.around(querschnitt.P_erf * utils.unit_conversion('N', 'MN') ,r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, 'P_erf,end [MN]')] = np.around(querschnitt.P_m0[ausgabe_an_knoten] * utils.unit_conversion('N', 'MN') ,r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, 'P_ist [MN]')] = np.around(querschnitt.P_ist_fuge[ausgabe_an_knoten] * utils.unit_conversion('N', 'MN') ,r_2)

            if self.include_tau:
                if self.mit_furnier:
                    self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Fe' + ' [' + querschnitt.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_Fe_design[ausgabe_an_knoten],r_2)    
                    self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_längs' + ' [' + querschnitt.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_längs_design[ausgabe_an_knoten],r_2)
                        
                    self.results_df.loc[:,(nabenhöhe, QS_label, 'Ausn. schub Furnier')] = np.around(querschnitt.ausnutzung_schub_Fe[ausgabe_an_knoten],r_3)
                    self.results_df.loc[:,(nabenhöhe, QS_label, 'Ausn. schub längs')] = np.around(querschnitt.ausnutzung_schub_längs[ausgabe_an_knoten],r_3)

                    self.results_df.loc[:,(nabenhöhe, QS_label, 'fvd_Fe')] = np.around(querschnitt.fvFed,r_3)
                    self.results_df.loc[:,(nabenhöhe, QS_label, 'fvd_längs')] = np.around(querschnitt.fvd_brutto,r_3)
                else:
                    self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Qy' + ' [' + querschnitt.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_Qy_design[ausgabe_an_knoten],r_2)
                    self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_Mx' + ' [' + querschnitt.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_Mx_design[ausgabe_an_knoten],r_2)
                    self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_xy' + ' [' + querschnitt.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_xy_design[ausgabe_an_knoten],r_2)    
                    self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['tau'] + '_tor' + ' [' + querschnitt.einheiten['Schubspannung'] + ']')] = np.around(querschnitt.tau_vtor_design[ausgabe_an_knoten],r_2)
                        
                    self.results_df.loc[:,(nabenhöhe, QS_label, 'Ausn. schub brutto')] = np.around(querschnitt.ausnutzung_schub_brutto[ausgabe_an_knoten],r_3)
                    self.results_df.loc[:,(nabenhöhe, QS_label, 'Ausn. schub netto')] = np.around(querschnitt.ausnutzung_schub_netto[ausgabe_an_knoten],r_3)
                    self.results_df.loc[:,(nabenhöhe, QS_label, 'Ausn. schub torsion')] = np.around(querschnitt.ausnutzung_schub_torsion[ausgabe_an_knoten],r_3)

                    self.results_df.loc[:,(nabenhöhe, QS_label, 'fvd_brutto')] = np.around(querschnitt.fvd_brutto,r_3)
                    self.results_df.loc[:,(nabenhöhe, QS_label, 'fvd_netto')] = np.around(querschnitt.fvxyd_netto,r_3)
                    self.results_df.loc[:,(nabenhöhe, QS_label, 'fvd_tor')] = np.around(querschnitt.fvtord,r_3)

            if self.include_reibung:
                self.results_df.loc[:,(nabenhöhe, QS_label, 'nxy_Qy' + ' [kN/m]')] = np.around(querschnitt.nxy_Qy[ausgabe_an_knoten] * utils.unit_conversion(querschnitt.einheiten['Grundschnittgrößen_f'], 'kN/m'),r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, 'nxy_Mx' + ' [kNm/m]')] = np.around(querschnitt.nxy_Mx[ausgabe_an_knoten]* utils.unit_conversion(querschnitt.einheiten['Grundschnittgrößen_m'], 'kNm/m'),r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, 'nxy_P,Rd' + ' [kN/m]')] = np.around(querschnitt.n_Rd[ausgabe_an_knoten]* utils.unit_conversion(querschnitt.einheiten['Grundschnittgrößen_f'], 'kN/m'),r_2)
                self.results_df.loc[:,(nabenhöhe, QS_label, 'Ausn. reibung')] = np.around(querschnitt.ausnutzung_reibung[ausgabe_an_knoten],r_3)


            if self.include_sigma_M_N:
                self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_N [' + querschnitt.einheiten['Normalspannung'] + ']')] = querschnitt.sigma_N_design[ausgabe_an_knoten]
                self.results_df.loc[:,(nabenhöhe, QS_label, GD.GREEK_UNICODE['sigma'] + '_M [' + querschnitt.einheiten['Normalspannung'] + ']')] = querschnitt.sigma_M_design[ausgabe_an_knoten]
        
            self.vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'Höhe [m]')] = np.around(querschnitt.section_absolute_heights['Ebenen'],r_1)
            self.vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'd_achse [m]')] = np.around(querschnitt.d_achse['Ebenen'],r_1)
            self.vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'P_erf,end [MN]')] = np.around(querschnitt.P_m0['Ebenen']* utils.unit_conversion('N', 'MN') ,r_1)
            self.vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n EX-84')] = querschnitt.n_ext_erf
            self.vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n Mono-5')] = querschnitt.n_int_pro_segment
            self.vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'n int. Summe')] = querschnitt.n_summe_int
            self.vorpsannungs_df.loc[:,(nabenhöhe, QS_label, 'P_ist Fuge [MN]')] = np.around(querschnitt.P_ist_fuge['Ebenen'] * utils.unit_conversion('N', 'MN') ,r_1)

            for segment in schnittgrößen_design['bauzustand'][QS_label][nabenhöhe]:
                if self.include_bauzustand_sgr:
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'M [MNm]')] = np.around(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['g']* utils.unit_conversion('Nm', 'MNm') ,r_2)
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'N [MN]')] = np.around(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['x']* utils.unit_conversion('N', 'MN') ,r_2)
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Q [MN]')] = np.around(schnittgrößen_design['bauzustand'][QS_label][nabenhöhe][segment]['y']* utils.unit_conversion('N', 'MN') ,r_2)
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Reibung [MN]')] = np.around(bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung'] * utils.unit_conversion('N', 'MN') ,r_2)
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'Ausn. reibung')] = np.around(bauzustands_kram[QS_label][nabenhöhe][segment]['Reibung Ausn.'],r_2)
                if self.include_bauzustand_spannung:
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, GD.GREEK_UNICODE['sigma'] + '_max [' + querschnitt.einheiten['Normalspannung'] + ']')] = np.around(sigma_bauzustand['sigma_max']['Ebenen'][segment] * querschnitt.einheiten_umrechnung['Normalspannung'], r_2)
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, GD.GREEK_UNICODE['sigma'] + '_min [' + querschnitt.einheiten['Normalspannung'] + ']')] = np.around(sigma_bauzustand['sigma_min']['Ebenen'][segment]* querschnitt.einheiten_umrechnung['Normalspannung'], r_2)
                if self.include_bauzustand_sgr_m:
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'nz_max [kN/m]')] = np.around(sigma_bauzustand['nz_max']['Ebenen'][segment] * utils.unit_conversion('N/m', 'kN/m') ,r_2)
                    self.bauzustands_df.loc[:,(nabenhöhe, QS_label, segment, 'nz_min [kN/m]')] = np.around(sigma_bauzustand['nz_min']['Ebenen'][segment] * utils.unit_conversion('N/m', 'kN/m') ,r_2)

        print ('\nAnmerkungen zu Annahmen bei der Berechnung:')
        for anmerkung in list(set(querschnitt.berechnungs_hinweise)):
            print (anmerkung)
        anmerkungs_df = pd.DataFrame({'Anmerkungen zu Annahmen bei der Berechnung':list(set(querschnitt.berechnungs_hinweise))})

        max_results_df = pd.DataFrame(max_ausnutzung)
        print ('\nMaximale Ausnutzungen')
        print (max_results_df)

    def write_dfs_to_excel(self, results_excel:str):
        '''
        results_excel: pfad inklusive datei name einer excel Datei
        '''
        #______________________ ALLES IN EXCEL SPEICHERN _________________________________________________________________________________________

        with pd.ExcelWriter(results_excel, mode= 'w', engine="openpyxl") as writer:# --> mode = 'a', if_sheet_exists='overlay'
            nrows = nachweis_parameter_df.shape[0]
            ausnutzungs_cols = utils.get_spalten_ausnutzung(results_df, df_results_header_list, start_row = nrows+2)
            nz_max_cols = utils.get_spalten_nach_name(bauzustands_df, df_bauzustand_header_list, start_row = 0, name = 'nz_max')
            reibung_cols = utils.get_spalten_nach_name(bauzustands_df, df_bauzustand_header_list, start_row = 0, name = 'Ausn. reibung')

            holz_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=0, startcol=0, index=False)#
            qs_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=4, startcol=0)
            lagen_df.to_excel(writer, sheet_name= 'QS_Werte', startrow=qs_df.shape[0] + 10, startcol=0)

            einwirkungs_parameter_df.to_excel(writer, sheet_name= 'Einwirkung_design', startrow=0, startcol=0, index=False)
            einwirkungs_df.to_excel(writer, sheet_name= 'Einwirkung_design', startrow=4, startcol=0)
            einwirkungs_df_typ.to_excel(writer, sheet_name= 'Einwirkung_design', startrow=einwirkungs_df.shape[0]+12, startcol=0)

            vorpsannungs_df.to_excel(writer, sheet_name= 'Vorspannung_Staffelung', startrow=4, startcol=0)

            anmerkungs_df.to_excel(writer, sheet_name= 'Hinweise_Berechnung', startrow=1, startcol=0, index=False)

            results_df.to_excel(writer, sheet_name= 'Berechnungs_Ergebnisse', startrow=nrows+2, startcol=0, index=True)#, float_format='%.4f')
            nachweis_parameter_df.to_excel(writer, sheet_name= 'Berechnungs_Ergebnisse', startrow=0, startcol=0, index=False)
            max_results_df.to_excel(writer, sheet_name= 'Ausnutzungen_max', startrow=1, index=True)

            bauzustands_df.to_excel(writer, sheet_name= 'Bauzustand', startcol=0, index=True)

        utils.zellen_groeße_formatieren(results_excel, worksheet= 'Einwirkung_design', cell_width=15, n_cols=len(einwirkungs_df.columns)+1)
        utils.zellen_groeße_formatieren(results_excel, worksheet= 'Vorspannung_Staffelung', cell_width=15, n_cols=len(vorpsannungs_df.columns)+1)
        utils.zellen_groeße_formatieren(results_excel, worksheet= 'QS_Werte', cell_width=13, n_cols=max(len(qs_df.columns), len(holz_df.columns))+1)
        utils.zellen_groeße_formatieren(results_excel, worksheet= 'Berechnungs_Ergebnisse', cell_width=17, n_cols=len(results_df.columns)+1)
        utils.zellen_groeße_formatieren(results_excel, worksheet= 'Ausnutzungen_max', cell_width=20, n_cols=len(max_results_df.columns)+1)
        utils.zellen_groeße_formatieren(results_excel, worksheet= 'Hinweise_Berechnung', cell_width=100, n_cols=1, start_col='A')
        utils.zellen_groeße_formatieren(results_excel, worksheet= 'Bauzustand', cell_width=15, n_cols=len(bauzustands_df.columns)+1, start_col='A')

        utils.add_databar_color(results_excel, worksheet = 'Berechnungs_Ergebnisse', columns = ausnutzungs_cols)
        utils.add_color_rule(results_excel, worksheet = 'Bauzustand', columns = nz_max_cols, regel = 'greaterThan', value = 0)
        utils.add_color_rule(results_excel, worksheet = 'Bauzustand', columns = reibung_cols, regel = 'greaterThan', value = 1)

        #utils.zelle_beschriften(results_excel, 'Ausnutzungen_max', 'B' + str(start_row_max_ti[t]), 
        #                       't, tX, tY [cm] ' + ', '.join([str(int(t)), str(round(tX,1)), str(round(tY,1))]) ,'B' + str(start_row_max_ti[t])+ ':E'+ str(start_row_max_ti[t]))
        dataframes_doc = {'QS_Werte_holz':holz_df,'QS_Werte_geometrie':qs_df,'QS_Werte_lagen':lagen_df,
                        'Einwirkungs_parameter':einwirkungs_parameter_df,'Einwirkung_nach_dauer':einwirkungs_df,'Einwirkung_nach_typ':einwirkungs_df_typ,
                        'Vorspannung':vorpsannungs_df,
                        'Anmerkungen':anmerkungs_df, 'Nachweis_parameter':nachweis_parameter_df,
                        'SGR_Spannungen':results_df,}
        utils.save_dataframes_to_pkl(dataframes_doc, ausgabe_an_knoten)


