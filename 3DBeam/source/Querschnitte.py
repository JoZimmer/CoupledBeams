from copy import copy
from math import cos, radians, sin
import pickle 
import matplotlib.pyplot as plt
import numpy as np
import source.utilities.utilities as utils
import source.utilities.global_definitions as GD

'''

Nachweisführung 

   Plattenbeanspruchung 
       Spannungen und Nachweise für 
       - aus Biegung
           - in Haupttragrichtung 
           - in Nebentragrichtung  
       - aus Schub 
           - in Haupttragrichtung
               benötigt: 
               äquivalente Schubflächen für BSP
               statisches Moment des Brettsperrholzelements in Haupttragrichtung                 
               Rollschubfestigkeit in schwerpunktsnächsten Querlage maßgebend 
               Schubfestigkeit der Längslagen  
           - in Nebentragrichtung
               statisches Moment des Brettsperrholzelements in Nebntragrichtung
       - Torsion bei Plattenbeanspruchung 

  

   Scheibenbeanspruchung 
       Spannungen und Nachweise für
       - Schub 
           -Abscheren der Bretter entlang einer Fuge 
               - Haupttragrichtung/Nebntragrichtung 
                   benötigt: 
                   T/As(min{A0,net und A90net)
           - Schubversagen in der Klebefläche in den Kreuzungspunkten 
                   benötigt: 
                   polares Trägheitsmoment/Brettbreite/ Mt im Brettsperrholzelement 
                   Anzahl der Klebeflächen/Anzahl der Klebefugen/ Kreuzungsfelder
           - Schubversagen der gesamten Scheibe --> Überschreitung der Rollschubfestigkeit


       - Knicken? 
       - Beulen? 

''' 

class Querschnitt(object):

    '''
    Klasse geschlossener Querschnitte aus BSP definiert durch den durchmesser der Achse
    '''

    def __init__(self, lagen_aufbau, holz_parameter:dict, nachweis_parameter:dict, hoehen_parameter:dict, 
                        einheiten:dict, FE_elements:int):
        '''
        d_achse: Durchmesser des Turms in der Achse des Querschnitts
        holz_parameter: dict siehe global definitions von einer Bestimmten Holzklasse/Art
        hfract, absolute_height
        FE_elements: anzahl an Elementen mit denen dann später die Berechnung durchgeführt wird
        hoehn_parameter: parameter entlang der Höhe die die genaue Geometrie festlegen
        '''
        self.hoehen_parameter = hoehen_parameter
        self.nabenhöhe = hoehen_parameter['nabenhöhe']
        self.FE_nodes = FE_elements +1

        self.wand_stärke = 0
        self.lagen_aufbau = lagen_aufbau
        self.t_laengslagen = sum([lage['ti'] for lage in lagen_aufbau if lage['ortho'] == 'X'])
        self.t_querlagen = sum([lage['ti'] for lage in lagen_aufbau if lage['ortho'] == 'Y'])

        for lage in lagen_aufbau:
            self.wand_stärke += lage['ti']

        self.initalize_höhen_parameter()

        self.wichte = holz_parameter['rhok']
        
        self.holz_parameter = holz_parameter
        self.nachweis_parameter= nachweis_parameter
        self.einheiten = einheiten      

        self.initialize_einheiten_umrechnung()

        self.berechnungs_hinweise = []


    def initialize_einheiten_umrechnung(self):
        self.einheiten_umrechnung = {'Normalspannung':1, 'Schubspannung':1}

        einheit_spannung = self.einheiten['Kraft'] + '/' + self.einheiten['Länge'] + '²'

        if einheit_spannung != self.einheiten['Festigkeit']:
            self.einheiten_umrechnung['Normalspannung'] =  utils.unit_conversion(einheit_spannung, self.einheiten['Festigkeit'])
            self.einheiten_umrechnung['Schubspannung'] =  utils.unit_conversion(einheit_spannung, self.einheiten['Festigkeit'])
            self.einheiten['Normalspannung'] = self.einheiten['Festigkeit']
            self.einheiten['Schubspannung'] = self.einheiten['Festigkeit']
        
        else:
            self.einheiten['Normalspannung'] = einheit_spannung
            self.einheiten['Schubspannung'] = einheit_spannung

        self.einheiten['Grundschnittgrößen_f'] = self.einheiten['Kraft'] + '/' + self.einheiten['Länge'] 
        self.einheiten['Grundschnittgrößen_m'] = self.einheiten['Kraft'] + self.einheiten['Länge'] + '/' + self.einheiten['Länge'] 

    def initalize_höhen_parameter(self):
        '''
        hoehen_parameter:
            'd_knick' = Durchmesser des Knicks (obere für Freigang relevanter. 
                        None: überhaupt kein Knick wird angeordnet. Durchmesser ändert sich linear von d_unten zu d_oben
            'h_knick_von_oben': Abstand des Knicks von der OK Turm (gemessen von oben) -> es Wird der Knoten genommen der am nächsten dran oben drüber liegt
            'höhe_sektionen': range in meter in der die Höhe der einzelnen Sektionen liegen soll -> wird anhand der Nabenhöhe gewählt
            'd_unten_oben': (Anfangs) Werte der Durchmesser am Fuß und am Kopf
            'd_knick_übergang': 
                2. Durchmesser wenn ein Übergnag stattfinden soll 
                automatisch: setzt den Durchmesser der an der stelle vorhanden ist fest und ändert den Fußdurchmesser
                None: macht gar nichts
            'n_sektionen_übergang': Anzahl an Sektionen unterhalb des Knicks in denen ein 2. Radius Verlauf stattfinden soll 
            'd_unten_angepasst': im Fall des 2. Knicks kann der untere Radius angepasst werden 
        NOTE: das ganze ist nicht sooo super variable -> Bisher her nur für den fall gedacht mit den 2 kicken, einer in die iene Richtung einer in die andere
        '''
        from scipy.optimize import minimize_scalar
        from functools import partial

        def func_hi(h, hi):
            return h/hi - int(h/hi)

        opt_func = partial(func_hi, self.nabenhöhe)
        bounds_sektionen_höhen = tuple(self.hoehen_parameter['höhe_sektionen'])
        minimization_result = minimize_scalar(opt_func,
                                            method='Bounded',
                                            bounds=bounds_sektionen_höhen)

        self.section_heights = {}                                   

        self.section_heights['Ebenen'] = minimization_result.x
        self.n_sections = int(self.nabenhöhe/self.section_heights['Ebenen'])
        self.n_ebenen = self.n_sections + 1
        self.section_absolute_heights = {'Ebenen':np.round(np.arange(0, self.nabenhöhe, self.section_heights['Ebenen']),2)}
        #self.section_absolute_heights = np.linspace(0, self.nabenhöhe, self.n_ebenen)
        self.hfract = {'Ebenen':self.section_absolute_heights['Ebenen']/self.section_absolute_heights['Ebenen'][-1]}

        print ('Für Nabenhöhe', self.nabenhöhe, 'ergeben sich', self.n_sections, 'Sektionen mit einer Höhe von je', 
                round(self.section_heights['Ebenen'],2), 'm')

        d_u = self.hoehen_parameter['d_unten_oben'][0]
        d_o = self.hoehen_parameter['d_unten_oben'][1]

        # _________KEIN KNICK

        if self.hoehen_parameter['d_knick'] == None:
            self.d_achse = {'Ebenen':np.linspace(d_u, d_o, self.n_ebenen)}

        # ________ KNICK OBEN
        else:
            h_knick_von_unten = self.section_absolute_heights[-1] - self.hoehen_parameter['h_knick_von_oben']
            diff_knoten_i = self.section_absolute_heights - h_knick_von_unten
            for i, dif in enumerate(diff_knoten_i):
                if dif < 0:
                    continue
                if dif > 0:
                    ebene_knick = i
                    break
            d_achse_o = np.linspace(self.hoehen_parameter['d_knick'], d_o, self.n_ebenen - ebene_knick)    
            d_achse_u = np.linspace(d_u, self.hoehen_parameter['d_knick'], ebene_knick+1)
            d_knick = self.hoehen_parameter['d_knick']

            self.hoehen_parameter['h_knick_von_unten'] = [self.section_absolute_heights[ebene_knick]]
            knick_ebenen = [ebene_knick]
            print ('Bei h=', round(self.section_absolute_heights[ebene_knick],2), 'm ist ein Knick angeordnet.')

            self.d_achse['Ebenen'] = np.concatenate((d_achse_u, d_achse_o[1:]))

        # ______ ÜBERGNAGS KNICK
            if self.hoehen_parameter['d_knick_übergang'] != None and self.hoehen_parameter['d_unten_angepasst'] != d_u:
                if self.hoehen_parameter['d_knick_übergang'] == 'automatisch':
                    d_knick_übergang = self.d_achse[ebene_knick - self.hoehen_parameter['n_sektionen_übergang']]
                else:
                    d_knick_übergang = self.hoehen_parameter['d_knick_übergang']

                self.hoehen_parameter['h_knick_von_unten'].append(self.section_absolute_heights[ebene_knick-self.hoehen_parameter['n_sektionen_übergang']])
                knick_ebenen.append(ebene_knick-self.hoehen_parameter['n_sektionen_übergang'])
                print ('Bei h=', round(self.section_absolute_heights[ebene_knick-self.hoehen_parameter['n_sektionen_übergang']],2), 'm ist ein 2. Knick angeordnet.')
                
                d_achse_u = np.linspace(self.hoehen_parameter['d_unten_angepasst'], d_knick_übergang, ebene_knick - self.hoehen_parameter['n_sektionen_übergang']+1)
                d_achse_übergang = np.linspace(d_knick_übergang, d_knick,  self.hoehen_parameter['n_sektionen_übergang']+1)
                self.d_achse['Ebenen'] = np.concatenate((d_achse_u[:-1], d_achse_übergang, d_achse_o[1:]))

            
            #for d_knick in self.hoehen_parameter['d_knick']:
            #knick_ebenen = [ebene_knick-self.hoehen_parameter['n_sektionen_übergang'], ebene_knick]
            self.hoehen_parameter['knick_ebenen'] = knick_ebenen

        # ________ JETZT DAS GANZE AUF EBENE DER FE KNOTEN BESTIMMEN
        self.section_absolute_heights['FE'] = np.round(np.linspace(0,self.section_absolute_heights['Ebenen'][-1], self.FE_nodes),2)
        self.section_heights['FE'] = np.diff(self.section_absolute_heights['FE'])[0]
        self.hfract['FE'] = self.section_absolute_heights['FE']/self.section_absolute_heights['FE'][-1]
        
        d_achse_FE = []

        for x in self.section_absolute_heights['FE']:
            for i, ebene in enumerate(self.d_achse['Ebenen'][:-1]):
                x1, x2 = self.section_absolute_heights['Ebenen'][i], self.section_absolute_heights['Ebenen'][i+1]
                if x1 <= x < x2:
                    y1 = self.d_achse['Ebenen'][i]
                    y2 = self.d_achse['Ebenen'][i+1]
                    val_x = np.interp(x, (x1,x2), (y1,y2))

                    d_achse_FE.append(val_x)

        d_achse_FE.append(self.d_achse['Ebenen'][-1])
        self.d_achse['FE'] = np.asarray(d_achse_FE)

        # ______ VON D_ACHSE ABHÄNIGE PARAMETER
        self.d_außen, self.d_innen = {}, {}
        for knoten in self.d_achse:
            self.d_außen[knoten] = self.d_achse[knoten] + self.wand_stärke
            self.d_innen[knoten] = self.d_achse[knoten] - self.wand_stärke



# _________ EFFEKTIVE Steifigkeitswerte PLATTENBEANSPRUCHUNG _________________________ 

    def compute_effektive_moment_of_inertia_platte_y(self):
        t_laengslagen= 0
        width= 1
        eigenanteil = 0
        steineranteil = 0
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'X':
                t_laengslagen+= lage['ti']
                eigenanteil += (lage['ti']**3)*width*(1/12)
                steineranteil += lage['ti']*width*(lage['a']**2)
        self.Iz_eff_platte= eigenanteil+steineranteil
        
    def compute_effektive_moment_of_inertia_platte_z(self):
        t_laengslagen= 0
        width= 1
        eigenanteil = 0
        steineranteil = 0
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'Y':
                t_laengslagen+= lage['ti']
                eigenanteil += (lage['ti']**3)*width*(1/12)
                steineranteil += lage['ti']*width*lage['a']**2
        self.Iz_eff_platte= eigenanteil+steineranteil

    def compute_effektive_moment_of_inertia_platte_Sxz(self):
        pass

    def compute_effektive_moment_of_inertia_platte_Syz(self):
        pass

 
# ________ FESTIGKEITEN/ WIDERTÄNDE ____________________________

    def compute_effektive_festigkeiten_charakteristisch(self):
        '''
        die Festigkeiten werden Abgemindert durch den Lagenaufbau 
        -> für den Nachweis muss daher die im gesamten QS wirkende Spannung dieser Abgeminderten Festigkeit gegenübergestellt werden
        NOTE ACHTUNG im Fall der Scheibenbeanspruchung ist das bisher andersrum eingebaut
        '''
        # effektive festigkeiten laengs 
        t_querlagen = 0 
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'Y':
                t_querlagen+= lage['ti']

        self.fmk_eff_laengs = (self.wand_stärke - t_querlagen)/self.wand_stärke * self.holz_parameter['fmk']
        self.ft0k_eff_laengs = (self.wand_stärke - t_querlagen)/self.wand_stärke * self.holz_parameter['ft0k']
        self.fc0k_eff_laengs = (self.wand_stärke - t_querlagen)/self.wand_stärke * self.holz_parameter['fc0k']
        self.fvk_eff_laengs = (self.wand_stärke - t_querlagen)/self.wand_stärke * self.holz_parameter['fvk']
    
        # effektive Festigkeiten quer
        t_laengslagen = 0 
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'X':
                t_laengslagen+= lage['ti']

        self.fmk_eff_quer = (self.wand_stärke - t_laengslagen)/self.wand_stärke * self.holz_parameter['fmk']
        self.ft0k_eff_quer = (self.wand_stärke - t_laengslagen)/self.wand_stärke * self.holz_parameter['ft0k']
        self.fc0k_eff_quer = (self.wand_stärke - t_laengslagen)/self.wand_stärke * self.holz_parameter['fc0k']
        self.fvk_eff_quer = (self.wand_stärke - t_laengslagen)/self.wand_stärke * self.holz_parameter['fvk']

        # Scheibfestigkeiten
        
        t_i_längs = [lage['ti'] for lage in self.lagen_aufbau if lage['ortho'] == 'X']
        Axnet = 0.8*(t_i_längs[0] + t_i_längs[-1]) + sum(t_i_längs[1:-1])
        Aynet = self.t_querlagen
        Axynet = min(Axnet, Aynet) # "schwächere" Scheibenebene Versagt

        self.fvk_brutto = self.holz_parameter['fvk']
        self.fvxyk_netto = self.holz_parameter['fvxyk'] * (Axynet / self.wand_stärke) # In norm und VO wird diese Abminderung auf Einwirkungsseite vorgenommen
        self.fvtork = self.holz_parameter['ftornodek'] * (Axynet / self.wand_stärke)

        # Schubfestigkeit einer Furnierebene
        if 'fvFek' in self.holz_parameter:
            self.fvFek = self.holz_parameter['fvFek']

    def compute_effektive_festigkeiten_design(self, einwirkungsdauer):
        '''
        TODO: Datenstruktur der Festigkeitskennwerte wäre vmtl auch in dict besser:
        f = {'charakterisitsch':{'v':..., 'c0':...},
             'design':{'v':..., 'c0':...}}
        '''


        self.compute_effektive_festigkeiten_charakteristisch()

        sicherheits_faktor = self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']

        #effektive Festigkeite laengs
        self.fmd_eff_laengs = self.fmk_eff_laengs * sicherheits_faktor
        self.ft0d_eff_laengs= self.ft0k_eff_laengs * sicherheits_faktor
        self.fc0d_eff_laengs= self.fc0k_eff_laengs * sicherheits_faktor
        self.fvd_eff_laengs= self.fvk_eff_laengs * sicherheits_faktor
    
        #effektive Festigkeiten Quer
        self.fmd_eff_quer = self.fmk_eff_quer * sicherheits_faktor
        self.ft0d_eff_quer= self.ft0k_eff_quer * sicherheits_faktor
        self.fc0d_eff_quer= self.fc0k_eff_quer * sicherheits_faktor
        self.fvd_eff_quer= self.fvk_eff_quer * sicherheits_faktor

        # Scheibefestigkeiten
        self.fvd_brutto = self.fvk_brutto * sicherheits_faktor
        self.fvxyd_netto = self.fvxyk_netto * sicherheits_faktor 
        self.fvtord = self.fvtork * sicherheits_faktor

        if 'fvFek' in self.holz_parameter:
            self.fvFed = self.fvFek * sicherheits_faktor
    
#__________ STATISCHES MOMENT______________________

    def compute_static_moment(self):
        '''
        wird von der jeweiligen child klasse gemacht
        '''
        pass

    def compute_Wx(self):
        '''
        wird von der jeweiligen child klasse gemacht
        '''
        pass

# _________ NORMALSPANNUNGEN nx - im Balken cosy nz__________________________________
    
    def calculate_normalspannung(self, SGR_design, add_vorspannkraft_grob = False, plot = False, ist_bauzustand = False):
        '''
        ergibt dictionaries mit der einwirkungsdauers Dauer als key:
            sigma_druck: Spannung aus Normalkraft und moment
            sigma_zug: Spannung aus moment abzüglich Normalkraft (Annahme immer Druckkraft)
            sigma_N: Spannung nur mit Normalkraft berechnet
            sigma_M: Spannung nur mit Moment berechnet
        '''
        self.sigma_druck, self.sigma_zug, self.sigma_N, self.sigma_M, self.sigma_druck_ohne_vorsp = {},{},{},{},{}
        self.sigma_druck_innen, self.sigma_zug_innen, self.sigma_M_innen , self.sigma_M_neg= {},{}, {}, {}
        self.nz_min, self.nz_max = {},{}

        knoten_typ = list(self.d_achse.keys())
        knoten_typ_sgr = list(self.d_achse.keys())
        if ist_bauzustand:
            knoten_typ = ['Ebenen']
            knoten_typ_sgr = ['bauzustand']

        for i, knoten in enumerate(knoten_typ):
            e = self.d_achse[knoten]/2
            self.berechnungs_hinweise.append('   - Normalspannungen werden als Mittelwert mit Hebelarm e = Radius_achse berechnet')
            e_innen = self.d_innen[knoten]/2

            self.sigma_druck[knoten], self.sigma_zug[knoten], self.sigma_N[knoten], self.sigma_M[knoten], self.sigma_druck_ohne_vorsp[knoten] = {},{},{},{},{}
            self.sigma_druck_innen[knoten], self.sigma_zug_innen[knoten], self.sigma_M_innen[knoten] , self.sigma_M_neg[knoten] = {},{}, {}, {}
            if ist_bauzustand:
                self.nz_min[knoten], self.nz_max[knoten] = {}, {}

            SGR_design_current = SGR_design[knoten_typ_sgr[i]][self.name][self.nabenhöhe]

            for dauer in SGR_design_current:
                # Wirkungsrichtung von M neutralisieren um es dann als Druck und zug aufbringen zu können

                SGR_design_current[dauer]['Mz'] = abs(SGR_design_current[dauer]['Mz'])

                self.sigma_druck[knoten][dauer] = -(SGR_design_current[dauer]['Mz'])/self.Iz[knoten] * e + SGR_design_current[dauer]['Nx'] / self.A[knoten]
                self.sigma_zug[knoten][dauer] = (SGR_design_current[dauer]['Mz'])/self.Iz[knoten] * e + SGR_design_current[dauer]['Nx'] / self.A[knoten]

                self.sigma_druck_ohne_vorsp[knoten][dauer] = -(SGR_design_current[dauer]['Mz'])/self.Iz[knoten] * e + SGR_design_current[dauer]['Nx'] / self.A[knoten]
                self.sigma_druck_innen[knoten][dauer] = -(SGR_design_current[dauer]['Mz'])/self.Iz[knoten] * e_innen + SGR_design_current[dauer]['Nx'] / self.A[knoten]
                self.sigma_zug_innen[knoten][dauer] = (SGR_design_current[dauer]['Mz'])/self.Iz[knoten] * e_innen + SGR_design_current[dauer]['Nx'] / self.A[knoten]


                self.sigma_N[knoten][dauer] = SGR_design_current[dauer]['Nx'] / self.A[knoten]
                self.sigma_M[knoten][dauer] = (SGR_design_current[dauer]['Mz'])/self.Iz[knoten] * e 
                self.sigma_M_neg[knoten][dauer] = -(SGR_design_current[dauer]['Mz'])/self.Iz[knoten] * e 
                self.sigma_M_innen[knoten][dauer] = (SGR_design_current[dauer]['Mz'])/self.Iz[knoten] * e_innen 

                if ist_bauzustand:
                    self.nz_min[knoten][dauer] = self.sigma_druck[knoten][dauer] * self.wand_stärke
                    self.nz_max[knoten][dauer] = self.sigma_zug[knoten][dauer] * self.wand_stärke
        

            if add_vorspannkraft_grob:
                self.sigma_druck[knoten]['ständig'] -= self.sigma_druck_P[knoten] # für die berechnung der ausnutzung
                self.sigma_druck[knoten]['egal'] -= self.sigma_druck_P[knoten] # für die ergebniss ausgabe
                self.sigma_zug[knoten]['ständig'] -= self.sigma_druck_P[knoten] # für die berechnung der ausnutzung
                self.sigma_zug[knoten]['egal'] -= self.sigma_druck_P[knoten] # für die ergebniss ausgabe

            
        if plot:
            knoten_typ = 'Ebenen'            
            
            # TODO oder input an welchem Knoten und auslagern in postprocess mit eingang Spannung und knoten 
            knoten = [0, list(abs(self.sigma_druck[knoten_typ]['egal'])).index(max(abs(self.sigma_druck[knoten_typ]['egal']))),-1]

            fig, ax = plt.subplots(ncols=len(knoten), figsize= (15, 4), sharey=True)
            ef = utils.unit_conversion('N/m²', 'N/mm²')
            fig.suptitle('Normalspannung im Schnitt '  + self.name)
            import matplotlib.patches as patch
            efl = utils.unit_conversion(self.einheiten['Länge'],'cm')
            x = np.array([0, self.wand_stärke, 3*self.wand_stärke, 4*self.wand_stärke])* efl
            dx = self.wand_stärke * efl

            for i, knoten_i in enumerate(knoten):
                dy = 0.4#self.sigma_druck['egal'][0] * ef * -1/10
                y0_patch = -dy/2
                p_links = patch.Rectangle((x[0],y0_patch), dx, dy, fill= None, hatch= '///' )
                p_rechts = patch.Rectangle((x[2],y0_patch), dx, dy, fill= None, hatch= '///' )
                ax[i].set_title('Knoten ' +str(knoten_i) + ' h='+ str(round(self.section_absolute_heights[knoten_typ][knoten_i],2)) + 'm')
                ax[i].add_patch(p_links)
                ax[i].add_patch(p_rechts)
                
                y = [self.sigma_zug[knoten_typ]['egal'][knoten_i]*ef,
                    self.sigma_zug_innen[knoten_typ]['egal'][knoten_i]*ef,
                    self.sigma_druck_innen[knoten_typ]['egal'][knoten_i]*ef,
                    self.sigma_druck_ohne_vorsp[knoten_typ]['egal'][knoten_i]*ef]

                # TODO wenn Vorspannung dazu dann so ca.:
                # y = np.array(y) - self.sigma_zug['egal'][knoten_i]*ef

                colors = {1:'tab:red', -1:'tab:blue'}
                for si, spannung in enumerate(y):
                    ax[i].plot([x[si],x[si]], [0, spannung], color = colors[np.sign(spannung)])
                
                for xi, yi in zip(x,y):
                    ax[i].annotate(str(round(yi,2)), xy=(xi, yi))# - sign*0.4))
                
                ax[i].plot([x[0],x[1]], [y[0], y[1]], color = colors[np.sign(y[1])])
                ax[i].plot([x[2],x[3]], [y[2], y[3]], color = colors[np.sign(y[3])])

                ax[i].set_xticks([])
                ax[i].set_ylabel(GD.GREEK_UNICODE['sigma'] + ' [N/mm²]')
                ax[i].grid()

            plt.tight_layout()
            plt.show()

    def calculate_ausnutzung_normalspannung(self, SGR_design, add_vorspannkraft_grob = False, plot_spannungsverlauf=False):
        '''
        add_vorspannkraft_grob: Spannkraftverluste in %
        Berechnung der Ausnutzungen für:
            sigma_druck: Spannung aus Normalkraft und moment
            sigma_zug: Spannung aus moment abzüglich Normalkraft(Annahme immer Druckkraft)
            sigma_N: Spannung nur mit Normalkraft berechnet
            sigma_M: Spannung nur mit Moment berechnet

        Ergibt:
            die design Spannungen unabhängig von der einwirkungsdauer in der Einheit der Festigkeit 
            die Ausnutzungen in unter Berücksichtigung der Einwirkungsdauern
        Lasten kommen schon als Design
        Lasten = Schnittgrößen als dict sortiert nach einwirkungsdauer
        '''

        self.calculate_normalspannung(SGR_design, add_vorspannkraft_grob, plot=plot_spannungsverlauf)

        # erst mal Einheiten klar machen
        einheiten_faktor = self.einheiten_umrechnung['Normalspannung']
        self.sigma_druck_design, self.sigma_zug_design, self.sigma_N_design, self.sigma_M_design = {}, {}, {}, {}
        self.ausnutzung_druck, self.ausnutzung_zug, self.ausnutzung_N, self.ausnutzung_M = {}, {}, {}, {}
        for knoten in self.d_achse:
            self.sigma_druck_design[knoten] = self.sigma_druck[knoten]['egal'] * einheiten_faktor
            self.sigma_zug_design[knoten] = self.sigma_zug[knoten]['egal'] * einheiten_faktor
            self.sigma_N_design[knoten] = self.sigma_N[knoten]['egal'] * einheiten_faktor
            self.sigma_M_design[knoten] = self.sigma_M[knoten]['egal'] * einheiten_faktor

            self.ausnutzung_druck[knoten], self.ausnutzung_zug[knoten], self.ausnutzung_N[knoten], self.ausnutzung_M[knoten] = 0,0,0,0
            for dauer in SGR_design[knoten][self.name][self.nabenhöhe]:
                if dauer == 'egal' or dauer == 'spannkraft':
                    continue

                self.compute_effektive_festigkeiten_design(dauer)

                self.ausnutzung_druck[knoten] += abs(self.sigma_druck[knoten][dauer]* einheiten_faktor)/self.fc0d_eff_laengs 
                self.ausnutzung_zug[knoten] += abs(self.sigma_zug[knoten][dauer]* einheiten_faktor)/self.ft0d_eff_laengs 
                self.ausnutzung_N[knoten] += abs(self.sigma_N[knoten][dauer]* einheiten_faktor)/self.fc0d_eff_laengs # Annahme nur Normalkraft
                self.ausnutzung_M[knoten] += abs(self.sigma_M[knoten][dauer]* einheiten_faktor)/self.fc0d_eff_laengs         

    def spannkraft_berechnung(self, SGR_design, spannglied_parameter, n_drähte = 60, n_jahre=20, feuchteänderung = 0, dT = 0, verluste_pauschal = 0, unit = 'MN'):
        '''
        fester_wert = Verluste in %

        Zusammenfassung:
        T:\Projekte\Projekte 2022\(E) Entwicklung\22-E-001 - Studie Windraeder\Vorspannung -> docx Dokument
        Bemessungswert der Vorspannkraft ergibt sich auch der Mittleren Vorspannkraft abzüglich der Verluste infolge:
            - Kriechen 
            - Schwinden des Wekrstoffs (kann vermieden werden wenn holz nicht weiter austrocknet)
            - Langzeitrelaxation des Spannstahls 
            - Ankerschlupf
            - Temperatur
            --> Ausgabe der Verluste als dictionary sortiert nach diesen Arten + die gesamte = 'gesamt'

        Teilsicherheitsbeiwerte: 
            - Wirkt die Vorspannung günstig: γP = 1.0
            - Wirkt sie ungünstig: γP,unfav = 1.3  NA: γP,unfav = 1.0 (bei Th. II Ordnung eines extern Vorgespannten Bauteils)

        Spanndraht SUSPA:
            T:\Projekte\Projekte 2022\(E) Entwicklung\22-E-001 - Studie Windraeder\Literatur\Vorspannung\2021_Zulassung_Spannverfahren2 ABG für SUSPA Draht für Windenergieanlagen.pdf

        TODO:
            - NKL ?
            - F_quasi-Ständig nur Eigengewicht? oder auch eine Belastung aus z.B. aus einer Mittleren Windgeschwindigkeit 
        '''
        def get_belastungsniveau():
            '''
            fc,quasi-ständig / fc,mean
            bis dahin nimm einfach > 30% als ungünstigsten fall
            '''
            return 0
         
        self.calculate_normalspannung(SGR_design)


        # NOTE das hier ist nur auf Knoten_typ Ebene zu machen. Zwischen den Ebenen ist die Kraft aus Vorspannung konstant
        knoten = 'Ebenen'
        self.U_außen[knoten] = self.d_außen[knoten] * np.pi
        #'spannkraft': hat die SGR berechnet mit den für die Spannkraft berechnung relevante einwirkungskombination und sicherheitsbeiwerten
        self.berechnungs_hinweise.append('   - Erforderliche Spannkraft berechnet abzüglich des Eigengewichts mit gamma_m = 1.0')
        self.P_erf = self.sigma_zug[knoten]['spannkraft'] * self.wand_stärke * self.U_außen[knoten] # dauer 

        n_stunden = n_jahre * 365 * 24

        self.spannkraft_verluste = {}

        self.P_m0 = {'Ebenen':copy(self.P_erf)}

        if not verluste_pauschal:
            '''
            irgendwann hier die detaillierte Verlustberechnung
            '''
            # _____KRIECHEN
            belastungs_niveau = 30 # TODO Funktion die diese bestimmt und 30 , 15 oder 'sonst' zurückgibt
            k_def = self.nachweis_parameter['k_def']['NKL1'][belastungs_niveau]
            A_netto = self.A_längslagen[knoten]
            # Kraft in Richtung der Spannkraft (bisher ist Mz ständig nur aus Imperfektion)
            F_quasi_ständig = abs(SGR_design[knoten][self.name][self.nabenhöhe]['ständig']['Nx']) + abs(SGR_design[knoten][self.name][self.nabenhöhe]['ständig']['Mz']/self.d_achse[knoten])
            d_epsilon_kr = (self.P_erf + F_quasi_ständig) /(self.holz_parameter['E0mean'] * A_netto) * k_def

            self.spannkraft_verluste['kriechen'] = n_drähte * spannglied_parameter['Ap'] * spannglied_parameter['Ep'] * d_epsilon_kr

            # _______SCHWINDEN
            # TODO : Geringe Feuchteänderung sollte das ziehl sein
            d_epsilon_s = self.holz_parameter['alpha'] * feuchteänderung
            self.spannkraft_verluste['kriechen'] = n_drähte * spannglied_parameter['Ap'] * spannglied_parameter['Ep'] * d_epsilon_s

            # ______ANKERSCHLUPF
            # TODO: schlupf dehnung abhängig von der gesamt länge der jewiligen spannglieder
            d_epsilon_0 = spannglied_parameter['schlupf'] / (self.section_absolute_heights[knoten] * utils.unit_conversion(self.einheiten['Länge'], 'mm'))
            self.spannkraft_verluste['ankerschlupf'] = n_drähte * spannglied_parameter['Ap'] * spannglied_parameter['Ep'] * d_epsilon_0

            # _____RELAXATION
            nu_spannglied = 1.0 # Ausnutzung des Spannstahls 
            rho1000 = spannglied_parameter['rho1000']
            d = 0.66 * rho1000 * np.exp(9.1*nu_spannglied) * ((n_stunden/1000) ** (0.75*(1-nu_spannglied))) * 1E-05
            self.spannkraft_verluste['relaxation'] = d * self.P_erf

            # _____TEMPERATURERHÖHUNG
            d_epsilon_T = (spannglied_parameter['alphaT'] - self.holz_parameter['alphaT']) * (dT + 273.15) # in Kelvin
            self.spannkraft_verluste['temperatur'] = n_drähte * spannglied_parameter['Ap'] * spannglied_parameter['Ep'] * d_epsilon_T

            # ____REIBUNGSVERLUSTE

            #______SUMME
            dP_ges = sum(list(self.spannkraft_verluste.values()))
            self.spannkraft_verluste['gesamt'] = dP_ges

            #self.P_end = self.P_m0 - dP_ges
            # for i, fugenkraft in enumerate(self.P_erf):
            #     if self.P_end[i] <= fugenkraft:
            #         print ('Aufgrund von Spannkraftverlusten ist in Fuge', i, 'Mehr Spannkraft erforderlich')
            # TODO hier müssen die Verluste jetzt auf die ursprüngliche Kraft draufgerechnet werden und dann alles noch mal neu gemacht werden 
            # das solange bis sich die Verluste nicht mehr nennenswert verändern

        else:
            self.spannkraft_verluste['gesamt'] = verluste_pauschal/100 * np.asarray(self.P_m0['Ebenen'])

        # die absoluten Verluste werden auf die gesamt Kraft aufaddiert
        self.P_m0['Ebenen'] += self.spannkraft_verluste['gesamt']

        self.P_m0['FE'] = np.zeros(self.FE_nodes)

        for i, x in enumerate(self.section_absolute_heights['FE']):
            for j, x_ebene in enumerate(self.section_absolute_heights['Ebenen'][1:]):
                j += 1
                if x >= self.section_absolute_heights['Ebenen'][j-1] and x <= x_ebene:
                    self.P_m0['FE'][i] += self.P_m0['Ebenen'][j-1]

        # NOTE NormalSpannung durch Vorspannung wird erst in spannglied_staffelung berechnet
        #self.sigma_druck_P = {'Ebenen': self.P_m0 / self.wand_stärke / self.U_außen[knoten]}
        
    def get_spannglied_staffelung(self, n_oberste_segmente:int, Pd_ext, Pd_int):
        '''
        n_oberste_segmente: Anzahl der oberen Segmente die durch die externe Vorspannung abgedeckt werden sollen 
        Pd: werte des gewählten spannsystems   
        das wird hier alles in abhängigkeit der Ebenen definition gemacht und dann für die FE definition interpoliert 
        return: Werte über Spanngliedanzahlen werden nur an den Ebenen Knoten zurück gegebn. Spannkraft werte an beiden Knotenebenen 
        '''
        from math import ceil

        P_oben = self.P_m0['Ebenen'][-n_oberste_segmente]
        n_ext_erf = ceil(P_oben/Pd_ext)
        P_ist = n_ext_erf * Pd_ext
        n_int_pro_segment = [0] * n_oberste_segmente
        P_ist_fuge = {'Ebenen':[P_ist]*n_oberste_segmente}
        self.P_ist_fuge = {'Ebenen':{}}

        for P_erf in self.P_m0['Ebenen'][-n_oberste_segmente-1::-1]:
            n_int_erf = ceil((P_erf-P_ist)/Pd_int)
            n_int_pro_segment.append(n_int_erf)
            P_ist += n_int_erf * Pd_int
            P_ist_fuge['Ebenen'].append(P_ist)

        n_int_summe = []
        n_sum_ist = 0
        for n in n_int_pro_segment:
            n_sum_ist += n
            n_int_summe.append(n_sum_ist)

        self.P_ist_fuge['Ebenen'] = np.asarray(list(reversed(P_ist_fuge['Ebenen'])))

        n_ext_erf_list = [0]*self.n_sections
        n_ext_erf_list.append(n_ext_erf)

        # die differenzen zu Pm0 jetzt noch auf sigma_druck_P drauf rechene -> das passiert jetzt wieder knotentyp weiße

        P_diff = self.P_ist_fuge['Ebenen'] - self.P_m0['Ebenen']

        # TODO ist control gleich dem was es soll
        # self.sigma_druck_P['FE'] checken das das richtig gemacht wird
        # += P_diff / self.wand_stärke / self.U_außen['Ebenen']
        self.sigma_druck_P= {'Ebenen':self.P_ist_fuge['Ebenen']/ self.wand_stärke / self.U_außen['Ebenen']}

        self.P_ist_fuge['FE'] = np.zeros(self.FE_nodes)

        for i, x in enumerate(self.section_absolute_heights['FE']):
            for j, x_ebene in enumerate(self.section_absolute_heights['Ebenen'][1:]):
                j += 1
                if x >= self.section_absolute_heights['Ebenen'][j-1] and x <= x_ebene:
                    self.P_ist_fuge['FE'][i] += self.P_ist_fuge['Ebenen'][j-1]

        self.sigma_druck_P['FE'] = self.P_ist_fuge['FE'] / self.wand_stärke / self.U_außen['FE']

        return n_ext_erf_list, list(reversed(n_int_pro_segment)), list(reversed(n_int_summe)), self.P_ist_fuge

# _________ Schalenbeanspruchung NORMALSPANNUNGEN ny__________________________________

    def calculate_normalspannung_y(self, windbelastung):

        t_querlagen= 0.08
        radius= self.d_achse/2

        n_y= radius* windbelastung['Fy']

        self.sigma_y= n_y/t_querlagen 

    def nachweis_normalspannung_y(self, windbelastung):
        self.calculate_normalspannung_y(windbelastung)

        einwirkung= self.sigma_y

        einheiten_faktor = self.einheiten_umrechnung['Schubspannung']
        
        self.ausnutzung_sigma_ny = self.sigma_y* einheiten_faktor/self.fc0d_eff_quer

# _________ PLATTENBEANSPRUCHUNG BIEGUNG ny__________________________________


    def calculate_normalspannung_Plattenbeanspruchung_Nebentragrichtung(self, wind_max, segmentbreite): 
        '''Längsspannungen bei Plattenbeanspruchung aus Windsog und -Druck ohne Berücksichtigung von Schalentragwirkung!'''
        
        
        self.compute_effektive_moment_of_inertia_platte_z()
        M0_k= (wind_max['Fy']*segmentbreite**2)/8
        self.sigma_ny_biegung= M0_k/self.Iz_eff_platte * self.wand_stärke/2

    def ausnutzung_Plattenbeanspruchung_Nebentragrichtung(self, wind_max, segmentbreite):
        
        self.calculate_normalspannung_Plattenbeanspruchung_Nebentragrichtung(wind_max, segmentbreite)

        einwirkung= self.sigma_ny_biegung

        einheiten_faktor = self.einheiten_umrechnung['Normalspannung']
        
        self.ausnutzung_sigma_ny_biegung= self.sigma_y* einheiten_faktor/self.fc0d_eff_quer

        return self.ausnutzung_sigma_ny_biegung

    def calculate_normalspannung_Plattenbeanspruchung_Haupttragrichtung():
        return 0

    def calculate_schub_plattenbeanspruchung(self, wind_max):

        '''
        nach neuer Norm muss für die schubfestigkeit bei Biegung ein Krack Faktor mit berücksichtigt werden
        '''
        width_segment= 3
        #in kN/m
        self.Q_platte= 1/2*wind_max*width_segment
        self.tau_platte= self.Q_platte/self.wand_stärke
        return 0

    def calculate_schub_knick(self):
        d_unten=12
        d_knick=5.5
        d_oben=3.4
        h_turm=110
        h_knick=50
        winkel= np.arctan((d_unten-d_knick)/h_knick)-np.arctan(d_knick-d_oben/h_turm-h_knick)

    def calculate_schub_plattenbeanspruchung_rollschub(self):
        '''Rollschubfestigkeit der schwerpunktnächsten Querlage maßgebend'''
        A_t_net=0 
        I0_net=0
        S0R_net=0

    def compute_drillsteifigkeit(self):
        self.K_xy= np.sqrt(self.Iz_eff_platte*self.holz_parameter['E0_mean']*self.Iz_eff_platte*self.holz_parameter['E90_mean'])

    
# _________ SCHUBSPANNUNGEN SCHEIBENBEANSPRUCHUNG nxy__________________________________

    def calculate_schubspannung_scheibe(self, SGR_design, lamellenbreite=0.13):
        '''
        aus Vorlesung 7 und prEN 1995-1-1_v3:

            tau_v,xy,d = nxy,d/A_xy,net
            n_xy,d = Schubfluss als Integral der Schubspannung von -t/2 bis t/2 oder F/l (VO wenn einzel Last)
            -> Festigkeit wird analog reduziert 
      
            tau_tor,node,d = 3/2 * tau_v,xy,d * (tl/bl)
                mit: tl = dünnste Lagendicke (bezeichnung in Vorlesung mit d)
                     bl = Brettbreite/Lamellenbreite bzw. Schwindnutabstand (l steht für Lamelle)
            ODER mittels Meachnik:
                Schubkraft * Abstand zur betrachteten Fuge ist maßgebendes Moment das aufgenommen werden muss
                Anzahl der Kreuzungsfelder -->Hersteller  
                Anzahl der Klebefugen 
        '''
        self.tau_Qy, self.tau_Mx, self.tau_xy, self.tau_vtor, self.tau_Fe, self.tau_längs = {},{},{}, {},{}, {}
        self.tau_Qy_design, self.tau_Mx_design, self.tau_xy_design,  self.tau_vtor_design, self.tau_Fe_design, self.tau_längs_design = {},{},{},{},{},{}

        if self.holz_parameter['Furnierebene']:
            self.berechnungs_hinweise.append('   - Verwende Furnierebenen Ansatz für die Schubspannungsberechnung')

        for knoten in self.d_achse:

            SGR_design_current = SGR_design[knoten][self.name][self.nabenhöhe]

            self.tau_Qy[knoten], self.tau_Mx[knoten], self.tau_xy[knoten], self.tau_vtor[knoten], self.tau_Fe[knoten], self.tau_längs[knoten] = {},{},{}, {},{}, {}
            
            for dauer in SGR_design_current:
                if self.holz_parameter['Furnierebene']:
                    self.Sy_max = self.compute_static_moment(d_achse = self.d_achse[knoten], mit_furnier = self.holz_parameter['Furnierebene'])
                    self.Wx = self.compute_Wx(d_achse = self.d_achse[knoten], mit_furnier = self.holz_parameter['Furnierebene'])
                    self.Iz = self.compute_flächenträgheitsmoment(d_achse = self.d_achse[knoten], mit_furnier= True)

                    n = self.holz_parameter['G4545']/self.holz_parameter['G0mean']
                    for lage in self.Sy_max:
                        if lage == 'furnier':
                            t = self.t_querlagen 
                            kraft_aufteilung = n * t / (self.t_laengslagen + n*t)

                            Qy = SGR_design_current[dauer]['Qy'] * kraft_aufteilung
                            Mx = SGR_design_current[dauer]['Mx'] * kraft_aufteilung
                            tau_Qy = (Qy * self.Sy_max[lage]) / (self.Iz[lage] * t) * self.einheiten_umrechnung['Schubspannung']
                            tau_Mx  = Mx /self.Wx[lage]  * self.einheiten_umrechnung['Schubspannung']

                            self.tau_Fe[knoten][dauer] = tau_Mx + tau_Qy

                        if lage=='längslagen':
                            t = self.t_laengslagen 
                            kraft_aufteilung = self.t_laengslagen / (self.t_laengslagen + n*self.t_querlagen)
                            Qy = SGR_design_current[dauer]['Qy'] * kraft_aufteilung
                            Mx = SGR_design_current[dauer]['Mx'] * kraft_aufteilung
                            tau_Qy = (Qy * self.Sy_max[lage]) / (self.Iz[lage] * t) * self.einheiten_umrechnung['Schubspannung']
                            tau_Mx  = Mx /self.Wx[lage]  * self.einheiten_umrechnung['Schubspannung']

                            self.tau_längs[knoten][dauer] = tau_Mx + tau_Qy

                # _____ NORMALES BSP    
                else:
                    # TODO mit wand_stärke richtig? oder sollte das nur mit querlagen sein? -> das wird wieder auf widerstandsseite berücksichtigt?!
                    self.tau_Qy[knoten][dauer] = (SGR_design_current[dauer]['Qy'] * self.Sy_max[knoten]) / (self.Iz[knoten] * self.wand_stärke) * self.einheiten_umrechnung['Schubspannung']
                    self.tau_Mx[knoten][dauer]  = SGR_design_current[dauer]['Mx'] /self.Wx[knoten]  * self.einheiten_umrechnung['Schubspannung']

                    self.tau_xy[knoten][dauer] = self.tau_Mx[knoten][dauer] + self.tau_Qy[knoten][dauer]
                    
                    # TODO tl in abhänigkeit bestimmen was dei maßgebende Schicht ist? also wenn Axy,net = Ay,net ist dann nur die querlagen betrachten? 
                    # Maximale Lamellen dicke 6 cm -> 12 cm Lagendicke setzt sich demnach aus 6+6 zusammen
                    tl = min(max([lage['ti'] for lage in self.lagen_aufbau]), 0.04)

                    self.tau_vtor[knoten][dauer] = 1.5 * self.tau_xy[knoten][dauer] * (tl/lamellenbreite)
                    self.berechnungs_hinweise.append('   - Torsion im Kreuzungspunkt der Lammellen in Folge Scheibenschub mit Lamellendicke ' + str(tl) + ' und Lamellenbreite ' + str(lamellenbreite) + ' berechnet')

            if not self.holz_parameter['Furnierebene']:
                self.tau_Qy_design[knoten] = self.tau_Qy[knoten]['egal']
                self.tau_Mx_design[knoten] = self.tau_Mx[knoten]['egal']
                self.tau_xy_design[knoten] = self.tau_xy[knoten]['egal']
                self.tau_vtor_design[knoten] = self.tau_vtor[knoten]['egal']
            else:
                self.tau_Fe_design[knoten] = self.tau_Fe[knoten]['egal']
                self.tau_längs_design[knoten] = self.tau_längs[knoten]['egal']
    
    def calculate_ausnutzung_scheibenschub(self, SGR_design, lamellenbreite=0.08):
        '''
        SGR_design = Schnittgrößen
        NOTE/TODO Einwirkungsdauer ist hier manuell auf 'kurz' gesetz
        lamellenbreite: Standartwert 0.08 -> TODO teste auswirkung

        Nachweis beinhaltet 3 Versagensmechanismen
            - Versagen des Gesamtquerschnitts = Brutto_schub
            - abscheren der Bretter (= Netto-Schub, tau_v,xy,d)
            - Torsion im Kreuzungspunkt (tau_tor,node,d)

        Quellen:
          T:\Organisation\00_AUSTAUSCH\zzz_Ehemalige\AH\Leitfaden BSP\BSP-LITERATUR -> Vorlesungsunterlagen BSP 
          Entwurf der Norm prEN 1995-1-1_v3
        '''

        self.calculate_schubspannung_scheibe(SGR_design, lamellenbreite)
        
        self.ausnutzung_schub_brutto, self.ausnutzung_schub_netto, self.ausnutzung_schub_torsion = {},{},{}
        self.ausnutzung_schub_Fe, self.ausnutzung_schub_längs = {},{}
        for knoten in self.d_achse:
            self.ausnutzung_schub_brutto[knoten], self.ausnutzung_schub_netto[knoten], self.ausnutzung_schub_torsion[knoten] = 0,0,0
            self.ausnutzung_schub_Fe[knoten], self.ausnutzung_schub_längs[knoten] = 0,0
            
            for dauer in SGR_design[knoten][self.name][self.nabenhöhe]:
                if dauer == 'egal' or dauer == 'spannkraft':
                    continue

                self.compute_effektive_festigkeiten_design(dauer)

                if self.holz_parameter['Furnierebene']:
                    self.ausnutzung_schub_Fe[knoten] += self.tau_Fe[knoten][dauer] / self.fvFed
                    self.ausnutzung_schub_längs[knoten] += self.tau_längs[knoten][dauer] / self.fvd_brutto

                else:
                    self.ausnutzung_schub_brutto[knoten] += self.tau_xy[knoten][dauer] / self.fvd_brutto
                    self.ausnutzung_schub_netto[knoten] += self.tau_xy[knoten][dauer] / self.fvxyd_netto
                    self.ausnutzung_schub_torsion[knoten] += self.tau_vtor[knoten][dauer] / self.fvtord

        self.berechnungs_hinweise.append('   - Schub Brutto Nachweis nur mit Schmalseitenverklebten Aufbauten (laut VO7)!')

    def calculate_schubkerrekturbeiwert(self, anzahl_lagen):
        '''
        k_n Abhängig von der Anzahl der Schichten und das Verhältnis der Schubmodule längs und quer zur Faser 
        Verhältnis der Schubmoduli bei Fichte 0,1
        Berechnung nach DIN
        Einfluss auf Schubsteifigkeit
        '''
        n=anzahl_lagen
        self.k_n= (11/(2*(n-1)))*((n-1)/20+(n+1)/2)

    #def calculate_rollschubnachweis(self):


# _________ KOMBINIERTE AUSNUTZUNG __________________________________

    def reibungs_nachweis_horizontalfugen(self, SGR_design, mu):
        '''
        aus Beton Kalender 2011 Kapitel 4.7.4
        Als vertikale, die Reibung aktivierende Kraft wird die gesamt Vorspannkraft je Fuge angenommen = Pm0
        Einwirkung sind die Schnittgrößen gesamt unabhängig der Einwirkungsdauer
        TODO: Oberste Fuge ist der Übergang von Holz zu Stahl
        '''
        # Vorzeichen werden hier vernachlässigt, da Pm0 ansich negative wirkt aber als absoluter wert gespeichert ist

        self.n_Rd, self.nxy_Qy, self.nxy_Mx, self.ausnutzung_reibung = {},{},{},{}
        for knoten in self.d_achse:
            try:
                Pd = self.P_ist_fuge[knoten]
                self.berechnungs_hinweise.append('   - Gestaffelte Spannglied Führung für Reibungsberechnung angesetzt')
            except NameError:
                Pd = self.P_m0[knoten]
                self.berechnungs_hinweise.append('\n   Erforderliche Vorspannkraft für Reibungsberechnung angesetzt. Staffelung nicht berücksichtigt')

            self.n_Rd[knoten] = mu * Pd / (np.pi * self.d_achse[knoten])

            self.nxy_Qy[knoten] = SGR_design[knoten][self.name][self.nabenhöhe]['egal']['Qy'] / (self.d_achse[knoten]/2)
            self.nxy_Mx[knoten] = SGR_design[knoten][self.name][self.nabenhöhe]['egal']['Mx'] / (2*np.pi*(self.d_achse[knoten]/2)**2)

            self.ausnutzung_reibung[knoten] = (self.nxy_Qy[knoten] + self.nxy_Mx[knoten]) / self.n_Rd[knoten]

# _________ KOMBINIERTE AUSNUTZUNG __________________________________

    def calculate_ausnutzung_kombiniert(self, LF):

        sigma_d_normal = self.calculate_normalspannung(LF)
        tau_q= self.calculate_schubspannung_querkraft(LF)
        tau_Mt= self.calculate_schubspannung_torsion(LF)
    
        self.compute_effektive_festigkeiten_design('kurz')

        fd = self.fvd_eff_laengs * (utils.unit_conversion(self.einheiten['Festigkeit'], self.einheiten['Normalspannung']))
        einwirkung_d=  (sigma_d_normal[0]+tau_q[0]+tau_Mt[0])

        nu = einwirkung_d/fd
        self.ausnutzung_kombiniert = nu
       
# _________ OBJEKT FUNKTIONEN __________________________________

    def compute_sectional_mean(self, parameter):
        
        if not isinstance(parameter, list):
            parameter = list(parameter)
        sectional_sum = [j+i for i, j in zip(parameter[:-1], parameter[1:])]
        return np.array(sectional_sum)/2

    def plot_properties_along_height(self, properties:list, units:list=None, E_modul=12000, prop_to_compare=None, h_absolut=True, knoten_typ = 'Ebenen'):
        '''
        werte entlang der Höhe plotten. 
        gewünschte parameter als strings in einer liste geben
        knoten_typ: 'FE' oder 'Ebenen' an welchen der knoten das alles gemacht werden soll
        E_modul = optional wenn EI berechnet werden soll
        'Iz','EIz', 'd_achse', 'a_außen', 'a_innen', 'd_außen','M'
        '''
        prop_values = {'Iz':self.Iz[knoten_typ],'EIz':self.Iz[knoten_typ] * E_modul, 'd_achse':self.d_achse[knoten_typ], 'd_außen':self.d_außen[knoten_typ],
                        }
        prop_units_default = {'Iz':'m^4', 'EIz':'Nm²', 'd_achse':'m', 'a_außen':'m', 'a_innen':'m', 'd_außen':'m', 'M':'kg'}

        fg, ax = plt.subplots(ncols = len(properties), figsize=(2,5))
        if len(properties) == 1:
            ax = [ax]

        if h_absolut:
            y = self.section_absolute_heights[knoten_typ]
            y_label = 'h [m]'
        else:
            y = self.hfract[knoten_typ]
            y_label = 'Hfract [-]'
    
        for p_i, prop in enumerate(properties):
            ax[p_i].plot(prop_values[prop], y)
            if prop_to_compare != None:
                ax[p_i].plot(prop_to_compare, y)
            ax[p_i].grid()
            ax[p_i].set_xlabel(prop + ' [' + prop_units_default[prop] +']')
        ax[0].set_ylabel(y_label)

        plt.show()

    def save_section_parameters(self):

        # self.section_parameters_mean = {
        #     'n_sections':len(self.Iz) - 1,
        #     'n_ebenen':len(self.Iz),
        #     'Iz': self.compute_sectional_mean(self.Iz),
        #     'd_achse':self.compute_sectional_mean(self.d_achse),
        #     'A':self.compute_sectional_mean(self.A),
        #     'section_absolute_heights':self.section_absolute_heights,
        #     'section_heights':self.section_heights
        # }

        self.section_parameters = {}

        for knoten in self.d_achse:
            self.section_parameters[knoten] = {
                'n_sections':len(self.Iz[knoten]) - 1,
                'n_ebenen':len(self.Iz[knoten]),
                'Iz': self.Iz[knoten],
                'd_achse':self.d_achse[knoten],
                'A':self.A[knoten],
                'section_absolute_heights':self.section_absolute_heights[knoten],
                'section_heights':self.section_heights[knoten]
            }

    def save_QS_parameters_charakteristisch(self):
        '''
        WICHTIG: uinterschied zu save_section_parameters: die Werte an den Knoten werden gespeichert und nicht die section mittelwerte
        Gedacht um die werte die für einen Bestimmten QS charakterisitisch sind zu speichern um sie dann in ein excel zu packen
        NOTE Wandstärke bisher konstant
        '''
        if 'transportbreite_max'  in self.hoehen_parameter:
            s_max = self.hoehen_parameter['transportbreite_max']
            bogen_max = self.d_außen['Ebenen'] * np.arcsin(s_max/self.d_außen['Ebenen'])
            n_ist = self.U_außen['Ebenen'] / bogen_max
            s_kleiner = self.d_außen['Ebenen'] * np.sin(self.U_außen['Ebenen']/(np.ceil(n_ist) * self.d_außen['Ebenen']))
            n_kleiner = np.ceil(n_ist)
            s_größer = self.d_außen['Ebenen'] * np.sin(self.U_außen['Ebenen']/(np.floor(n_ist) * self.d_außen['Ebenen']))
            n_größer = np.floor(n_ist)


            n_segment_info = {'b_max [m]':bogen_max,
                              's <' + str(self.hoehen_parameter['transportbreite_max']) : s_kleiner,
                              'n Seg. max':n_kleiner,
                              's >' + str(self.hoehen_parameter['transportbreite_max']) : s_größer,
                              'n Seg. min':n_größer
                             }

        self.querschnitts_werte = {}
        self.querschnitts_werte_all = {}
        for knoten in self.d_achse:
            V_sections = np.append(self.V[knoten], 0)
            if self.holz_parameter['Furnierebene']:
                Iz_lagen = self.compute_flächenträgheitsmoment(self.d_achse[knoten], mit_furnier = True)
                Iz_eff = np.add(Iz_lagen['längslagen'], Iz_lagen['furnier'])
            else:
                Iz_lagen = self.compute_flächenträgheitsmoment(self.d_achse[knoten], mit_furnier = False, nur_längs = True)
                Iz_eff = Iz_lagen['längslagen']
            self.querschnitts_werte[knoten] = {
                'Höhe [m]':self.section_absolute_heights[knoten],
                'Hfract':self.hfract[knoten],
                'd_achse [m]':self.d_achse[knoten],
                'd_außen [m]':self.d_außen[knoten],
                'Iz [m^4]': self.Iz[knoten],
                'A [m²]':self.A[knoten],
                'U_achse [m]':self.U[knoten],
                'U_außen [m]':self.U_außen[knoten],
                't [m]':np.asarray([self.wand_stärke]*len(self.d_achse[knoten])),
                'V_sections [m³]':V_sections,
                'MassDen [kg/m]':self.mass_density[knoten],
                'EIz,eff [Nm^2]': self.holz_parameter['E0mean'] * Iz_eff,
            }

            self.querschnitts_werte_all[knoten] = {
                'Höhe [m]':self.section_absolute_heights[knoten],
                'Hfract':self.hfract[knoten],
                'd_achse [m]':self.d_achse[knoten],
                'd_außen [m]':self.d_außen[knoten],
                'Iz [m^4]': self.Iz[knoten],
                'A [m²]':self.A[knoten],
                'U_achse [m]':self.U[knoten],
                'U_außen [m]':self.U_außen[knoten],
                't [m]':np.asarray([self.wand_stärke]*len(self.d_achse[knoten])),
                'V [m³]':self.V[knoten],
                'MassDen [kg/m]':self.mass_density[knoten],
                'Eigengewicht [kg]':self.eigengewicht[knoten],
                'Gewichtskraft [N]':self.gewichtskraft[knoten],
                'EIz,eff [Nm^2]': self.holz_parameter['E0mean'] * Iz_eff,
                'Sy_max [m^3]': self.Sy_max[knoten],
                'cd': np.asarray([self.cd]*len(self.d_achse[knoten]))
            }

            if 'transportbreite_max'  in self.hoehen_parameter:
                if knoten == 'Ebenen':
                    self.querschnitts_werte[knoten].update(n_segment_info)

        return self.querschnitts_werte

    def export_object_to_pkl(self, dest_file):

        with open(dest_file, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print ('\nSaved object in', dest_file)

    # _________ FATIGUE SPEZIFISCH__________________________________

    def calculate_normalspannung_at_knoten(self, SGR, knoten, knoten_typ):
        '''
        HIER: SGR sind Zeitreihen 

        ergibt normalspannung an einem Bestimmten knoten i für SGR = {'Mz':float, 'Nx':float}
        return: sigma_druck, sigma_zug (ansich min und max beide druck da N immer größer)
        knoten_typ: Ebenen oder FE
        '''
   
        e = self.d_achse[knoten_typ][knoten]/2
        #self.berechnungs_hinweise.append('   - Normalspannungen werden als Mittelwert mit Hebelarm e = Radius_achse berechnet')

        # Wirkungsrichtung von M neutralisieren um es dann als Druck und zug aufbringen zu können

        SGR['Mz'] = abs(SGR['Mz'])

        sigma_druck = -(SGR['Mz'])/self.Iz[knoten_typ][knoten] * e + SGR['Nx'] / self.A[knoten_typ][knoten]
        sigma_zug = (SGR['Mz'])/self.Iz[knoten_typ][knoten] * e + SGR['Nx'] / self.A[knoten_typ][knoten]

        return sigma_druck#, sigma_zug
        
    def calculate_schubspannung_scheibe_at_knoten(self, SGR, knoten, lamellenbreite=0.13, knoten_typ='Ebenen'):
        '''
        HIER: SGR sind Zeitreihen 

        aus Vorlesung 7 und prEN 1995-1-1_v3:

            tau_v,xy,d = nxy,d/A_xy,net
            n_xy,d = Schubfluss als Integral der Schubspannung von -t/2 bis t/2 oder F/l (VO wenn einzel Last)
            -> Festigkeit wird analog reduziert 
      
            tau_tor,node,d = 3/2 * tau_v,xy,d * (tl/bl)
                mit: tl = dünnste Lagendicke (bezeichnung in Vorlesung mit d)
                     bl = Brettbreite/Lamellenbreite bzw. Schwindnutabstand (l steht für Lamelle)
            ODER mittels Meachnik:
                Schubkraft * Abstand zur betrachteten Fuge ist maßgebendes Moment das aufgenommen werden muss
                Anzahl der Kreuzungsfelder -->Hersteller  
                Anzahl der Klebefugen 
        '''
        self.tau_Qy, self.tau_Mx, self.tau_xy, self.tau_vtor, self.tau_Fe, self.tau_längs = {},{},{}, {},{}, {}
        self.tau_Qy_design, self.tau_Mx_design, self.tau_xy_design,  self.tau_vtor_design, self.tau_Fe_design, self.tau_längs_design = {},{},{},{},{},{}

        if self.holz_parameter['Furnierebene']:
            self.berechnungs_hinweise.append('   - Verwende Furnierebenen Ansatz für die Schubspannungsberechnung')

        #for knoten in self.d_achse:
        SGR_design_current = SGR

        self.tau_Qy[knoten_typ], self.tau_Mx[knoten_typ], self.tau_xy[knoten_typ], self.tau_vtor[knoten_typ], self.tau_Fe[knoten_typ], self.tau_längs[knoten_typ] = {},{},{}, {},{}, {}
        
        #for dauer in SGR_design_current:
        if self.holz_parameter['Furnierebene']:
            self.Sy_max = self.compute_static_moment(d_achse = self.d_achse[knoten_typ], mit_furnier = self.holz_parameter['Furnierebene'])
            self.Wx = self.compute_Wx(d_achse = self.d_achse[knoten_typ], mit_furnier = self.holz_parameter['Furnierebene'])
            self.Iz = self.compute_flächenträgheitsmoment(d_achse = self.d_achse[knoten_typ], mit_furnier= True)

            n = self.holz_parameter['G4545']/self.holz_parameter['G0mean']
            for lage in self.Sy_max:
                if lage == 'furnier':
                    t = self.t_querlagen 
                    kraft_aufteilung = n * t / (self.t_laengslagen + n*t)

                    Qy = SGR_design_current['Qy'] * kraft_aufteilung
                    Mx = SGR_design_current['Mx'] * kraft_aufteilung
                    tau_Qy = (Qy * self.Sy_max[lage][knoten]) / (self.Iz[lage][knoten] * t) #* self.einheiten_umrechnung['Schubspannung']
                    tau_Mx  = Mx /self.Wx[lage][knoten]  #* self.einheiten_umrechnung['Schubspannung']

                    self.tau_Fe[knoten_typ] = tau_Mx + tau_Qy

                if lage=='längslagen':
                    t = self.t_laengslagen #sum([lage['ti'] for lage in self.lagen_aufbau if lage['ortho'] == 'Y']) # ist schon in self.t_querlagen
                    kraft_aufteilung = self.t_laengslagen / (self.t_laengslagen + n*self.t_querlagen)
                    Qy = SGR_design_current['Qy'] * kraft_aufteilung
                    Mx = SGR_design_current['Mx'] * kraft_aufteilung
                    tau_Qy = (Qy * self.Sy_max[lage][knoten]) / (self.Iz[lage][knoten] * t) #* self.einheiten_umrechnung['Schubspannung']
                    tau_Mx  = Mx /self.Wx[lage][knoten]  #* self.einheiten_umrechnung['Schubspannung']

                    self.tau_längs[knoten_typ] = tau_Mx + tau_Qy

        # _____ NORMALES BSP    
        else:
            self.tau_Qy[knoten_typ] = (SGR_design_current['Qy'] * self.Sy_max[knoten_typ][knoten]) / (self.Iz[knoten_typ][knoten] * self.wand_stärke) #* self.einheiten_umrechnung['Schubspannung']
            self.tau_Mx[knoten_typ]  = SGR_design_current['Mx'] /self.Wx[knoten_typ][knoten]  #* self.einheiten_umrechnung['Schubspannung']

            self.tau_xy[knoten_typ] = self.tau_Mx[knoten_typ] + self.tau_Qy[knoten_typ]
            
            # TODO tl in abhänigkeit bestimmen was dei maßgebende Schicht ist? also wenn Axy,net = Ay,net ist dann nur die querlagen betrachten? 
            # Maximale Lamellen dicke 6 cm -> 12 cm Lagendicke setzt sich demnach aus 6+6 zusammen
            tl = min(max([lage['ti'] for lage in self.lagen_aufbau]), 0.04)

            self.tau_vtor[knoten_typ] = 1.5 * self.tau_xy[knoten_typ] * (tl/lamellenbreite)
            self.berechnungs_hinweise.append('   - Torsion im Kreuzungspunkt der Lammellen in Folge Scheibenschub mit Lamellendicke ' + str(tl) + ' und Lamellenbreite ' + str(lamellenbreite) + ' berechnet')

        if not self.holz_parameter['Furnierebene']:
            self.tau_Qy_design[knoten_typ] = self.tau_Qy[knoten_typ]
            self.tau_Mx_design[knoten_typ] = self.tau_Mx[knoten_typ]
            self.tau_xy_design[knoten_typ] = self.tau_xy[knoten_typ]
            self.tau_vtor_design[knoten_typ] = self.tau_vtor[knoten_typ]
        else:
            self.tau_Fe_design[knoten_typ] = self.tau_Fe[knoten_typ]
            self.tau_längs_design[knoten_typ] = self.tau_längs[knoten_typ]

class KreisRing(Querschnitt):

    def __init__(self, cd = 1.0, cp_max =1.8 , lagen_aufbau = None, holz_parameter = {}, nachweis_parameter= {}, 
                hoehen_parameter ={}, einheiten = {}, FE_elements=10) -> None:
        
        super().__init__(lagen_aufbau, holz_parameter, nachweis_parameter, hoehen_parameter, einheiten,FE_elements)


        self.name = 'Ring'
        self.nachweis_parameter= nachweis_parameter

        self.A, self.U,self.U_außen, self.V, self.Iz, self.Sy_max, self.Wx, self.mass_density, self.eigengewicht, self.gewichtskraft = {},{},{},{},{},{},{},{},{},{}
        for knoten in self.d_achse:

            self.A[knoten] = self.compute_area(self.d_außen[knoten]/2)-self.compute_area(self.d_innen[knoten]/2)
            self.U[knoten] = self.d_achse[knoten] * np.pi
            self.U_außen[knoten] = self.d_außen[knoten] * np.pi

            self.V[knoten] = self.compute_sectional_mean(self.A[knoten]) * self.section_heights[knoten]
            self.Iz[knoten] = self.compute_flächenträgheitsmoment(self.d_achse[knoten])
            self.Sy_max[knoten] = self.compute_static_moment(self.d_achse[knoten])
            self.Wx[knoten] = self.compute_Wx(self.d_achse[knoten])

            self.mass_density[knoten] = self.A[knoten] * self.wichte
            self.eigengewicht[knoten] =  - np.append(self.V[knoten] * self.wichte, 0)
            gewichtskraft = -self.V[knoten] * self.wichte * GD.GRAVITY

            if self.holz_parameter['Furnierebene']:# berücksichtigung der unterschiedlichen wichten
                self.compute_geometrie_werte_Furnier()
                gewichtskraft -= self.V_Fe[knoten] * (self.holz_parameter['rhok_Fe'] - self.wichte) * GD.GRAVITY

            self.gewichtskraft[knoten] = {'Fx':np.append(gewichtskraft, 0)} # gewichtskraft wird immer unten am Element/Ebene angegeben

        #self.initalize_geometrie_parameter_FE()
        self.compute_netto_flächen()
        self.cd = cd
        self.cp_max = cp_max
        
        self.save_section_parameters()
        self.save_QS_parameters_charakteristisch()


    def compute_flächenträgheitsmoment(self, d_achse, mit_furnier = False, nur_längs=False):

        '''
        bisher nur für konstate Wandstärke über die Höhe
        d_achse wird hier als argument gegeben um es über unterschiedliche Höghen, also längen des arrays berechnen zu können
        '''

        if mit_furnier:
            #if 'di' not in self.lagen_aufbau[0]:
            utils.get_d_achse_lagen(d_achse, self.lagen_aufbau)
            Iz = {'furnier':0, 'längslagen':0}
            for i, lage in enumerate(self.lagen_aufbau):
                d_außen = lage['di'] + lage['ti']/2
                d_innen = lage['di'] - lage['ti']/2
                if lage['ortho'] == 'X':
                    Iz['längslagen'] += np.pi / 64 * (d_außen**4 - d_innen**4)
                if lage['ortho'] == 'Y':
                    Iz['furnier'] += np.pi / 64 * (d_außen**4 - d_innen**4)
            return Iz
        
        elif nur_längs:
            #if 'di' not in self.lagen_aufbau[0]:
            utils.get_d_achse_lagen(d_achse, self.lagen_aufbau)
            Iz = {'quer':0, 'längslagen':0}
            for i, lage in enumerate(self.lagen_aufbau):
                d_außen = lage['di'] + lage['ti']/2
                d_innen = lage['di'] - lage['ti']/2
                if lage['ortho'] == 'X':
                    Iz['längslagen'] += np.pi / 64 * (d_außen**4 - d_innen**4)
                if lage['ortho'] == 'Y':
                    Iz['quer'] += 0 / 64 * (d_außen**4 - d_innen**4)
            return Iz
        
        else:
            d_außen = d_achse + self.wand_stärke
            d_innen = d_achse - self.wand_stärke
            Iz = np.pi / 64 * (d_außen**4 - d_innen**4)
            return Iz

    def compute_area(self, r):
        area= np.pi * r**2

        return area

    def compute_netto_flächen(self):

        self.A_längslagen, self.A_querlagen = {},{}
        for knoten in self.d_achse:
            self.A_längslagen[knoten] = 0
            self.A_querlagen[knoten] = 0

            for i, lage in enumerate(self.lagen_aufbau):
                t_von_innen = sum([li['ti'] for li in self.lagen_aufbau[:i]])
                r_i_innen = self.d_innen[knoten]/2 + t_von_innen
                r_i_außen = self.d_innen[knoten]/2 + t_von_innen + lage['ti']
                if lage['ortho'] == 'X':
                    self.A_längslagen[knoten] += self.compute_area(r_i_außen) - self.compute_area(r_i_innen)
                elif lage['ortho'] == 'Y':
                    self.A_längslagen[knoten] += self.compute_area(r_i_außen) - self.compute_area(r_i_innen)

    def compute_static_moment(self, d_achse, mit_furnier = False):
        '''
        Sy berechnene mit nur den Längslagen?!
        my Sy bei z=0 (Viertelkreis)
        Berechnung für Kreisring 
        TODO: Achtung, hier wird die datenstruktur von Klassenvariablen verändert ggf.
        '''
        if mit_furnier:
            utils.get_d_achse_lagen(d_achse, self.lagen_aufbau)
            Sy_max = {'furnier':0, 'längslagen':0}
            for i, lage in enumerate(self.lagen_aufbau):
                if lage['ortho'] == 'X':
                    Sy_max['längslagen'] += ((lage['di']/2)**2)*lage['ti']
                if lage['ortho'] == 'Y':
                    Sy_max['furnier'] += ((lage['di']/2)**2)*lage['ti']

            return Sy_max

        else:
            Sy_max= ((d_achse/2)**2)*self.wand_stärke 
            return Sy_max

    def compute_Wx(self, d_achse, mit_furnier = False):
        '''
        TODO: Achtung, hier wird die datenstruktur von Klassenvariablen verändert ggf.
        '''

        if mit_furnier:
            #if 'di' not in self.lagen_aufbau[0]:
            utils.get_d_achse_lagen(d_achse, self.lagen_aufbau)
            Wx = {'furnier':0, 'längslagen':0}
            for i, lage in enumerate(self.lagen_aufbau):
                if lage['ortho'] == 'X':
                    Wx['längslagen'] += 2 * self.compute_area(lage['di']/2) * lage['ti']
                if lage['ortho'] == 'Y':
                    Wx['furnier'] += 2 * self.compute_area(lage['di']/2) * lage['ti']

            return Wx
        else:
            A_m= self.compute_area(d_achse/2)
            Wx = 2 * A_m * self.wand_stärke
            return Wx

    def initalize_geometrie_parameter_FE(self):
        '''
        die Querschnittswerte anhand der höhen parameter auch auf die anzahl an FE Elementen hin linear interpolieren
        '''
        self.x_FE = np.linspace(0,self.section_absolute_heights[-1], self.FE_nodes)
        self.section_heights_FE = np.diff(self.x_FE)
        geometrie_parameter = [self.d_achse, self.A]
        self.d_achse_FE, self.A_FE = [],[]
        geometrie_parameter_FE = [self.d_achse_FE, self.A_FE]

        for param_i, parameter in enumerate(geometrie_parameter):
            for x in self.x_FE:
                for i, ebene in enumerate(parameter[:-1]):
                    x1, x2 = self.section_absolute_heights[i], self.section_absolute_heights[i+1]
                    if x1 <= x < x2:
                        y1 = parameter[i]
                        y2 = parameter[i+1]
                        val_x = np.interp(x, (x1,x2), (y1,y2))

                        geometrie_parameter_FE[param_i].append(val_x)

            geometrie_parameter_FE[param_i].append(parameter[-1])
            geometrie_parameter_FE[param_i] = np.asarray(geometrie_parameter_FE[param_i])

        self.V_FE = self.compute_sectional_mean(self.A_FE) * self.section_heights_FE
        self.eigengewicht_FE = - np.append(self.V_FE *  self.wichte, 0)
        self.gewichtskraft_FE = {'Fx': - np.append(self.V_FE *  self.wichte * GD.GRAVITY, 0)}

        if self.holz_parameter['Furnierebene']:
            utils.get_d_achse_lagen(np.asarray(self.d_achse_FE), self.lagen_aufbau)
            self.A_Fe_FE=0
            for i, lage in enumerate(self.lagen_aufbau):
                if lage['ortho'] == 'Y':
                    d_a = lage['di'] + lage['ti']/2
                    d_i = lage['di'] - lage['ti']/2
                    self.A_Fe_FE += self.compute_area(d_a/2) - self.compute_area(d_i/2)

            self.V_Fe_FE = self.compute_sectional_mean(self.A_Fe_FE) * self.section_heights_FE
            self.gewichtskraft_FE['Fx'][:-1] -= self.V_Fe_FE * (self.holz_parameter['rhok_Fe'] - self.wichte) * GD.GRAVITY

    def compute_geometrie_werte_Furnier(self):
        '''
        d_achse, A, V der Furniere 
        '''
        self.A_Fe, self.V_Fe = {},{}
        for knoten in self.d_achse:
            #if 'di' not in self.lagen_aufbau[0]:
            utils.get_d_achse_lagen(self.d_achse[knoten], self.lagen_aufbau)
            self.A_Fe[knoten] = 0
            for i, lage in enumerate(self.lagen_aufbau):
                if lage['ortho'] == 'Y':
                    d_a = lage['di'] + lage['ti']/2
                    d_i = lage['di'] - lage['ti']/2
                    self.A_Fe[knoten] += self.compute_area(d_a/2) - self.compute_area(d_i/2)
            self.V_Fe[knoten] = self.compute_sectional_mean(self.A_Fe[knoten]) * self.section_heights[knoten]

class nEck(Querschnitt):

    def __init__(self, n_ecken, cd = 1.5, cp_max= 2.1, lagen_aufbau = None, holz_parameter = {}, nachweis_parameter= {}, hoehen_parameter ={}, einheiten = {}):
        '''
        werte können für einzelne sections oder als arrays gegeben werden
        Geometrie: https://de.wikipedia.org/wiki/Regelm%C3%A4%C3%9Figes_Polygon  
        '''

        super().__init__(lagen_aufbau, holz_parameter, nachweis_parameter, hoehen_parameter, einheiten)

        self.n_ecken = n_ecken
        self.name = str(self.n_ecken) + '-Eck'
        self.alpha = 360/self.n_ecken
        self.nachweis_parameter= nachweis_parameter
        self.a_außen = sin(radians(self.alpha/2)) * self.d_außen
        self.a_innen = sin(radians(self.alpha/2)) * self.d_innen

        self.compute_winkel_term()

        Iz_außen = self.compute_flächenträgheitsmoment_neck(self.a_außen)
        Iz_innen = self.compute_flächenträgheitsmoment_neck(self.a_innen)
        
        self.Iz = Iz_außen - Iz_innen

        # if self.lagen_aufbau:
        #     self.compute_effective_moment_of_inertia()

        self.A_außen = self.compute_area_neck(self.d_außen/2)
        self.A_innen = self.compute_area_neck(self.d_innen/2)
        self.A_m= self.compute_area_neck(self.d_achse/2)

        self.A = self.A_außen-self.A_innen

        self.cd = cd
        self.cp_max=cp_max
        self.save_section_parameters()
        self.save_QS_parameters_charakteristisch()


    def compute_winkel_term(self):
        
        self.winkel_term = ((2 + cos(radians(self.alpha))) / (1-cos(radians(self.alpha)))**2 ) * sin(radians(self.alpha))

    def compute_area_neck(self, r_u):
        ''' 
        r_u: Umkreis radius des n-ecks
        '''

        A = self.n_ecken/2 *  r_u**2 * sin(2*np.pi/self.n_ecken)

        return A

    def compute_flächenträgheitsmoment_neck(self, a):

        self.Iz = self.n_ecken/96 * a**4 * self.winkel_term
        return self.Iz

    def compute_sectional_properties(self):

        self.section_volume = self.compute_sectional_mean(self.A_außen) * self.section_heights

        mass_total = np.sum(self.section_mass)
        volume_total = np.sum(self.section_volume)
        self.mass_density = mass_total/volume_total

    def compute_effective_moment_of_inertia(self):
        '''
        Flächenträgheitsmoment der Längslagen
        '''
        for i, lage in enumerate(self.lagen_aufbau):
            nEck_section = nEck(self.n_ecken, lage['di'], lage['ti'], holz_parameter= self.holz_parameter, hoehen_parameter=self.hoehen_parameter,
                                einheiten =self.einheiten)
            lage['Iz'] = nEck_section.Iz

        self.Iz_eff = 0
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'X':
                self.Iz_eff += lage['Iz']

    def compute_static_moment(self):
        e=self.n_ecken/4   # da maximales statisches moment bei einem viertel des querschnitts
        r=self.d_achse/2
        Sy_max=0
        z_before= r         #beginnend an der Ecke 
        length= 2*r*sin(np.pi/self.n_ecken)
        for i in range(int(e)): 
            alpha_i= (2*np.pi/self.n_ecken)*(i+1) 
            z_e= (cos(alpha_i)*r) 
            Sy_max+= (z_e*length+((z_before-z_e)/2)*length)*self.wand_stärke
            z_before = z_e 
        self.Sy_max=Sy_max   
        return Sy_max


def plot_properties_along_height_list(geometric_objects:list, properties:list, units:list=None, E_modul=12000, prop_to_compare=None):
    '''
    liste an nEck Objekten
    werte entlang der Höhe plotten. 
    gewünschte parameter als strings in einer liste geben
    E_modul = optional wenn EI berechnet werden soll
    hfract muss hier für definiert werden

    '''
    fg, ax = plt.subplots(ncols = len(properties))
    if len(properties) == 1:
        ax = [ax]

    legend = ['gesamt QS', 'Längslagen']

    factors = [geometric_objects[-1].Iz_eff/geometric_objects[0].Iz_eff, geometric_objects[-1].Iz_eff/geometric_objects[1].Iz_eff, geometric_objects[1].Iz_eff/geometric_objects[0].Iz_eff]
    factors_label = ['EI12/EI8','EI12/EI10','EI10/EI8']

    for f_i , factor in enumerate(factors):
        print (factors_label[f_i] + ':', round(factor[0],2))

    for e_i, cross_section in enumerate(geometric_objects):
        prop_values = {'Iz':cross_section.Iz,'EIz':cross_section.Iz * E_modul, 'Iz_eff':cross_section.Iz_eff,'EIz_eff':cross_section.Iz_eff * E_modul,
                       'a_außen':cross_section.a_außen, 'a_innen':cross_section.a_innen, 'd_außen':cross_section.d_außen,
                        }
        prop_units_default = {'Iz':'m^4', 'EIz':'Nm²', 'Iz_eff':'m^4', 'EIz_eff':'Nm²', 'a_außen':'m', 'a_innen':'m', 'd_außen':'m', 'M':'kg'}

        for p_i, prop in enumerate(properties):
            if isinstance(cross_section, nEck):
                ax[p_i].plot(prop_values[prop], cross_section.hfract, label = str(cross_section.n_ecken) + 'Ecken')
            else:
                ax[p_i].plot(prop_values[prop], cross_section.hfract, label = 'Kreisring')
            #label = factors_label[e_i]
            #fa = factors[e_i]
            #ax[p_i].plot(factors[e_i], cross_section.hfract, label = factors_label[e_i])
            if prop_to_compare != None:
                ax[p_i].plot(prop_to_compare, cross_section.hfract)
                
            ax[p_i].grid()
            ax[p_i].set_xlabel(prop + ' [' + prop_units_default[prop] +']')
        ax[0].set_ylabel('Hfract [-]')
        ax[0].grid()
    
    #f = round(geometric_objects[1].Iz[0]/geometric_objects[0].Iz[0],2)
    #ax[0].plot(0,0, label = 'factor '+ str(f))

    plt.grid()
    plt.legend()
    plt.show()
