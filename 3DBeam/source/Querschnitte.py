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
        self.initalize_höhen_parameter()

        self.wand_stärke = 0
        self.lagen_aufbau = lagen_aufbau
        self.t_laengslagen = sum([lage['ti'] for lage in lagen_aufbau if lage['ortho'] == 'X'])
        self.t_querlagen = sum([lage['ti'] for lage in lagen_aufbau if lage['ortho'] == 'Y'])

        for lage in lagen_aufbau:
            self.wand_stärke += lage['ti']
        
        self.wichte = holz_parameter['rhok']

        self.d_außen = self.d_achse + self.wand_stärke
        self.d_innen = self.d_achse - self.wand_stärke

        self.FE_nodes = FE_elements +1

        self.holz_parameter = holz_parameter
        self.nachweis_parameter= nachweis_parameter
        self.einheiten = einheiten      

        self.initialize_einheiten_umrechnung()


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

        self.section_heights = minimization_result.x
        self.n_sections = int(self.nabenhöhe/self.section_heights)
        self.n_ebenen = self.n_sections + 1
        self.section_absolute_heights = np.arange(0, self.nabenhöhe, self.section_heights)
        #self.section_absolute_heights = np.linspace(0, self.nabenhöhe, self.n_ebenen)
        self.hfract = self.section_absolute_heights/self.section_absolute_heights[-1]

        print ('Für Nabenhöhe', self.nabenhöhe, 'ergeben sich', self.n_sections, 'Sektionen mit einer Höhe von je', 
                round(self.section_heights,2), 'm')

        d_u = self.hoehen_parameter['d_unten_oben'][0]
        d_o = self.hoehen_parameter['d_unten_oben'][1]

        # _________KEIN KNICK

        if self.hoehen_parameter['d_knick'] == None:
            self.d_achse = np.linspace(d_u, d_o, self.n_ebenen)

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

            self.d_achse = np.concatenate((d_achse_u, d_achse_o[1:]))

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
                self.d_achse= np.concatenate((d_achse_u[:-1], d_achse_übergang, d_achse_o[1:]))

            
            #for d_knick in self.hoehen_parameter['d_knick']:
            #knick_ebenen = [ebene_knick-self.hoehen_parameter['n_sektionen_übergang'], ebene_knick]
            self.hoehen_parameter['knick_ebenen'] = knick_ebenen

        print()


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
                #print(lage['ti'])
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
                #print(lage['ti'])
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
        self.fvxyk_netto = self.holz_parameter['fvxyk'] * (Axynet / self.wand_stärke) # In norm und VO wird diese Abmonderung auf Einwirkungsseite vorgenommen
        self.fvtork = self.holz_parameter['ftornodek'] * (Axynet / self.wand_stärke)

    def compute_effektive_festigkeiten_design(self, einwirkungsdauer):


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
    
#__________ STATISCHES MOMENT______________________

    def compute_static_moment(self):
        '''
        wird von der jeweiligen child klasse gemacht
        '''
        pass

# _________ NORMALSPANNUNGEN nx__________________________________
    
    def calculate_normalspannung(self, SGR_design, add_vorspannkraft_grob = False, plot = False):
        '''
        ergibt dictionaries mit der einwirkungsdauers Dauer als key:
            sigma_druck: Spannung aus Normalkraft und moment
            sigma_zug: Spannung aus moment abzüglich Normalkraft (Annahme immer Druckkraft)
            sigma_N: Spannung nur mit Normalkraft berechnet
            sigma_M: Spannung nur mit Moment berechnet
        '''
        
        e = self.d_außen/2
        e_innen = self.d_innen/2#, self.d_achse/2]

        self.sigma_druck, self.sigma_zug, self.sigma_N, self.sigma_M, self.sigma_druck_ohne_vorsp = {},{},{},{},{}
        self.sigma_druck_innen, self.sigma_zug_innen, self.sigma_M_innen , self.sigma_M_neg= {},{}, {}, {}

        for dauer in SGR_design:
            # Wirkungsrichtung von M neutralisieren um es dann als Druck und zug aufbringen zu können

            SGR_design[dauer]['Mz'] = abs(SGR_design[dauer]['Mz'])

            self.sigma_druck[dauer] = -(SGR_design[dauer]['Mz'])/self.Iz * e + SGR_design[dauer]['Nx'] / self.A
            self.sigma_zug[dauer] = (SGR_design[dauer]['Mz'])/self.Iz * e + SGR_design[dauer]['Nx'] / self.A

            self.sigma_druck_ohne_vorsp[dauer] = -(SGR_design[dauer]['Mz'])/self.Iz * e + SGR_design[dauer]['Nx'] / self.A
            self.sigma_druck_innen[dauer] = -(SGR_design[dauer]['Mz'])/self.Iz * e_innen + SGR_design[dauer]['Nx'] / self.A
            self.sigma_zug_innen[dauer] = (SGR_design[dauer]['Mz'])/self.Iz * e_innen + SGR_design[dauer]['Nx'] / self.A


            self.sigma_N[dauer] = SGR_design[dauer]['Nx'] / self.A
            self.sigma_M[dauer] = (SGR_design[dauer]['Mz'])/self.Iz * e 
            self.sigma_M_neg[dauer] = -(SGR_design[dauer]['Mz'])/self.Iz * e 
            self.sigma_M_innen[dauer] = (SGR_design[dauer]['Mz'])/self.Iz * e_innen 

        if add_vorspannkraft_grob:
            self.sigma_druck['ständig'] -= self.sigma_druck_P # np.negative(self.sigma_zug['egal']) * (1+add_vorspannkraft_grob/100) # für die berechnung der ausnutzung
            self.sigma_druck['egal'] -= self.sigma_druck_P # np.negative(self.sigma_zug['egal']) * (1+add_vorspannkraft_grob/100) # für die ergebniss ausgabe

            self.U_außen = self.d_außen * np.pi

            self.P_control = np.negative(self.sigma_zug['egal']) * (1+add_vorspannkraft_grob/100) * self.wand_stärke * self.U_außen # * utils.unit_conversion(self.einheiten['Kraft'], unit)

        if plot:            
            
            # TODO oder input an welchem Knoten und auslagern in postprocess mit eingang Spannung und knoten 
            knoten = [0, list(abs(self.sigma_druck['egal'])).index(max(abs(self.sigma_druck['egal']))),-1]

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
                ax[i].set_title('Knoten ' +str(knoten_i) + ' h='+ str(round(self.section_absolute_heights[knoten_i],2)) + 'm')
                ax[i].add_patch(p_links)
                ax[i].add_patch(p_rechts)
                
                y = [self.sigma_zug['egal'][knoten_i]*ef,
                    self.sigma_zug_innen['egal'][knoten_i]*ef,
                    self.sigma_druck_innen['egal'][knoten_i]*ef,
                    self.sigma_druck_ohne_vorsp['egal'][knoten_i]*ef]

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

    def calculate_ausnutzung_normalspannung(self, lasten_design, add_vorspannkraft_grob = False, plot_spannungsverlauf=False):
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

        self.calculate_normalspannung(lasten_design, add_vorspannkraft_grob, plot=plot_spannungsverlauf)

        # erst mal Einheiten klar machen
        einheiten_faktor = self.einheiten_umrechnung['Normalspannung']
        self.sigma_druck_design = self.sigma_druck['egal'] * einheiten_faktor
        self.sigma_zug_design = self.sigma_zug['egal'] * einheiten_faktor
        self.sigma_N_design = self.sigma_N['egal'] * einheiten_faktor
        self.sigma_M_design = self.sigma_M['egal'] * einheiten_faktor

        self.ausnutzung_druck, self.ausnutzung_zug, self.ausnutzung_N, self.ausnutzung_M = 0,0,0,0
        for dauer in lasten_design:
            if dauer == 'egal' or dauer == 'spannkraft':
                continue

            self.compute_effektive_festigkeiten_design(dauer)

            self.ausnutzung_druck += abs(self.sigma_druck[dauer]* einheiten_faktor)/self.fc0d_eff_laengs 
            self.ausnutzung_zug += abs(self.sigma_zug[dauer]* einheiten_faktor)/self.ft0d_eff_laengs 
            self.ausnutzung_N += abs(self.sigma_N[dauer]* einheiten_faktor)/self.fc0d_eff_laengs # Annahme nur Normalkraft
            self.ausnutzung_M += abs(self.sigma_M[dauer]* einheiten_faktor)/self.fc0d_eff_laengs         

    def spannkraft_berechnung(self, lasten_design, spannglied_parameter, n_drähte = 60, n_jahre=20, feuchteänderung = 0, dT = 0, verluste_pauschal = 0, unit = 'MN'):
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
        # eigengewicht entfernen um zugspannugen nicht 
        self.calculate_normalspannung(lasten_design)

        self.U_außen = self.d_außen * np.pi

        self.P_erf = self.sigma_zug['spannkraft'] * self.wand_stärke * self.U_außen # dauer 'spannkraft': hat die SGR berechnet mit den für die Spannkraft berechnung relevante einwirkungskombination

        n_stunden = n_jahre * 365 * 24

        self.spannkraft_verluste = {}

        #self.P_m0 = self.get_spannglied_staffelung() # TODO fürs erste -> an dieser stelle sollte eine genau spangleidwahl und staffelung erstellt werden
        self.P_m0 = copy(self.P_erf)

        if not verluste_pauschal:
            '''
            irgendwann hier die detaillierte Verlustberechnung
            '''
            # _____KRIECHEN
            belastungs_niveau = 30 # TODO Funktion die diese bestimmt und 30 , 15 oder 'sonst' zurückgibt
            k_def = self.nachweis_parameter['k_def']['NKL1'][belastungs_niveau]
            A_netto = self.A_längslagen
            # Kraft in Richtung der Spannkraft (bisher ist Mz ständig nur aus Imperfektion)
            F_quasi_ständig = abs(lasten_design['ständig']['Nx']) + abs(lasten_design['ständig']['Mz']/self.d_achse)
            d_epsilon_kr = (self.P_erf + F_quasi_ständig) /(self.holz_parameter['E0mean'] * A_netto) * k_def

            self.spannkraft_verluste['kriechen'] = n_drähte * spannglied_parameter['Ap'] * spannglied_parameter['Ep'] * d_epsilon_kr

            # _______SCHWINDEN
            # TODO : Geringe Feuchteänderung sollte das ziehl sein
            d_epsilon_s = self.holz_parameter['alpha'] * feuchteänderung
            self.spannkraft_verluste['kriechen'] = n_drähte * spannglied_parameter['Ap'] * spannglied_parameter['Ep'] * d_epsilon_s

            # ______ANKERSCHLUPF
            # TODO: schlupf dehnung abhängig von der gesamt länge der jewiligen spannglieder
            d_epsilon_0 = spannglied_parameter['schlupf'] / (self.section_absolute_heights * utils.unit_conversion(self.einheiten['Länge'], 'mm'))
            self.spannkraft_verluste['ankerschlupf'] = n_drähte * spannglied_parameter['Ap'] * spannglied_parameter['Ep'] * d_epsilon_0

            # _____RELAXATION
            nu_spannglied = 1.0 # Ausnutzung des Spannstahls 
            rho1000 = spannglied_parameter['rho1000']
            d = 0.66 * rho1000 * np.exp(9.1*nu_spannglied) * ((n_stunden/1000) ** (0.75*(1-nu_spannglied))) * 1E-05
            self.spannkraft_verluste['relaxation'] = d * self.P_erf

            # _____TEMPERATURERHÖHUNG
            d_epsilon_T = (spannglied_parameter['alphaT'] - self.holz_parameter['alphaT']) * (dT + 273.15) # in Kelvin
            self.spannkraft_verluste['temperatur'] = n_drähte * spannglied_parameter['Ap'] * spannglied_parameter['Ep'] * d_epsilon_T

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
            self.spannkraft_verluste['gesamt'] = verluste_pauschal/100 * np.asarray(self.P_m0)

        # die absoluten Verluste werden auf die gesamt Kraft aufaddiert
        self.P_m0 += self.spannkraft_verluste['gesamt']

        self.sigma_druck_P = self.P_m0 / self.wand_stärke / self.U_außen

        # TODO wieder in Spannungen umrechenen und dem Nachweis hinzufügen
        
    def get_spannglied_staffelung(self, überschuss_toleranz = 2):
        '''
        P_ext -> externe Vorpsannkraft die alle Fugen von oben bis zum Knick überdrücken muss
        überschuss_toleranz = Differenz anhand derer die Staffelung der internen vorspannung gewählt wird [MN]. 
                              ist die Differenz der erfoderderlichen Überdrückkraft zwischen 2 Fugen größer der toleranz wird gestaffelt
        TODO stimmt noch nicht ganz
        '''
        P_m0 = [0] # soll quasi von oben nach unten gefüllt werden und am ende solang wie die anzahl an steffelungen sein
        if self.hoehen_parameter['d_knick'][0] != None:
            staffelungs_ebenen = copy(self.hoehen_parameter['knick_ebenen'])
            #staffelungs_ebenen.insert(0, self.n_ebenen)
            for ebene in staffelungs_ebenen:
                P_ext = self.P_erf[ebene] - sum(P_m0) # in dieser Liste kommt erst der obere dann der übergangsknick
                P_m0.append(P_ext)
        
        # Gerader Turm
        else:
            staffelungs_ebenen = [self.n_ebenen]

        staffelungs_ebenen.append(staffelungs_ebenen[-1]-1)
        # # KOMPLIZIERT
        # 
        # while staffelungs_ebenen[-1] != 0:
        #     P_erf_curr = self.P_erf[staffelungs_ebenen[-1]] - sum(P_m0)
        #     dif = self.P_erf[staffelungs_ebenen[-1]-1] - P_erf_curr
        #     if self.P_erf[staffelungs_ebenen[-1]-1] - P_erf_curr >= überschuss_toleranz *utils.unit_conversion('MN',self.einheiten['Kraft']):
        #         staffelungs_ebenen.append(staffelungs_ebenen[-1]-1)
        #         P_m0.append(P_erf_curr)
        #         su = sum(P_m0)
        #     else:
        #         staffelungs_ebenen[-1] -= 1

        # # JEDE FUGE 

        while staffelungs_ebenen[-1] != 0:
            P_erf_curr = self.P_erf[staffelungs_ebenen[-1]] - sum(P_m0)
            P_m0.append(P_erf_curr)
            staffelungs_ebenen.append(staffelungs_ebenen[-1]-1)
        
        sum_P = sum(P_m0)
        
        return P_m0

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
        #trägt nicht über die lange Seite (12m) ab.. Außerdem Schalenwirkung

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

    def calculate_schubspannung_scheibe(self, lasten_design, lamellenbreite=0.08):
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

        self.compute_static_moment()
        self.Wx = 2 * self.A_m * self.wand_stärke

        self.tau_Qy, self.tau_Mx, self.tau_xy, self.tau_vtor = {},{},{}, {}
        for dauer in lasten_design:

            self.tau_Qy[dauer] = (lasten_design[dauer]['Qy'] * self.Sy_max) / (self.Iz * self.wand_stärke) * self.einheiten_umrechnung['Schubspannung']
            self.tau_Mx[dauer]  = lasten_design[dauer]['Mx'] /self.Wx  * self.einheiten_umrechnung['Schubspannung']

            self.tau_xy[dauer] = self.tau_Mx[dauer] + self.tau_Qy[dauer]
            
            # TODO tl in abhänigkeit bestimmen was dei maßgebende Schicht ist? also wenn Axy,net = Ay,net ist dann nur die querlagen betrachten? 
            # Maximale Lamellen dicke 6 cm -> 12 cm Lagendicke setzt sich demnach aus 6+6 zusammen
            tl = min(max([lage['ti'] for lage in self.lagen_aufbau]), 0.06)

            self.tau_vtor[dauer] = 1.5 * self.tau_xy[dauer] * (tl/lamellenbreite)

        self.tau_Qy_design = self.tau_Qy['egal']
        self.tau_Mx_design = self.tau_Mx['egal']
        self.tau_xy_design = self.tau_xy['egal']
        self.tau_vtor_design = self.tau_vtor['egal']
    

        # # Alt bzw. die "mechanische Alternative zum Trosionsversagen der Klebefuge"
        # n_k= 10*50
        # lamellenbreite= 0.2
        # self.Ip_Brett= lamellenbreite**4/6
        # Mt_fuge=lamellenbreite/2*(self.tau_Mx+self.tau_Qy)/self.wand_stärke

        # tau_Mt_fuge= (3*Mt_fuge/n_k+lamellenbreite**3)

    def calculate_ausnutzung_scheibenschub(self, lasten_design, lamellenbreite=0.08):
        '''
        lasten_design = Schnittgrößen
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

        self.calculate_schubspannung_scheibe(lasten_design, lamellenbreite)
        
        self.ausnutzung_schub_brutto, self.ausnutzung_schub_netto, self.ausnutzung_schub_torsion = 0,0,0
        
        for dauer in lasten_design:
            if dauer == 'egal' or dauer == 'spannkraft':
                continue

            self.compute_effektive_festigkeiten_design(dauer)

            self.ausnutzung_schub_brutto += self.tau_xy[dauer] / self.fvd_brutto
            self.ausnutzung_schub_netto += self.tau_xy[dauer] / self.fvxyd_netto
            self.ausnutzung_schub_torsion += self.tau_vtor[dauer] / self.fvtord

        print ('\nSchub Brutto Nachweis nur mit Schmalseitenverklebten Aufbauten (laut VO7)??!')

        
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

    def plot_properties_along_height(self, properties:list, units:list=None, E_modul=12000, prop_to_compare=None, h_absolut=False):
        '''
        werte entlang der Höhe plotten. 
        gewünschte parameter als strings in einer liste geben
        E_modul = optional wenn EI berechnet werden soll
        'Iz','EIz', 'd_achse', 'a_außen', 'a_innen', 'd_außen','M'
        '''
        prop_values = {'Iz':self.Iz,'EIz':self.Iz * E_modul, 'd_achse':self.d_achse,'d_achseFE':self.d_achse_FE, 'd_außen':self.d_außen,
                        'M':self.masse_pro_meter}
        prop_units_default = {'Iz':'m^4', 'EIz':'Nm²', 'd_achse':'m', 'd_achseFE':'m', 'a_außen':'m', 'a_innen':'m', 'd_außen':'m', 'M':'kg'}

        fg, ax = plt.subplots(ncols = len(properties), figsize=(2,5))
        if len(properties) == 1:
            ax = [ax]

        if h_absolut:
            y = self.section_absolute_heights
            y_label = 'h [m]'
        else:
            y = self.hfract
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

        self.section_parameters_mean = {
            'n_sections':len(self.Iz) - 1,
            'n_ebenen':len(self.Iz),
            'Iz': self.compute_sectional_mean(self.Iz),
            'd_achse':self.compute_sectional_mean(self.d_achse),
            'A':self.compute_sectional_mean(self.A),
            'masse_pro_m':self.masse_pro_meter,
            'section_absolute_heights':self.section_absolute_heights,
            'section_heights':self.section_heights
        }
        self.section_parameters = {
            'n_sections':len(self.Iz) - 1,
            'n_ebenen':len(self.Iz),
            'Iz': self.Iz,
            'd_achse':self.d_achse,
            'A':self.A,
            'masse_pro_m':self.masse_pro_meter,
            'section_absolute_heights':self.section_absolute_heights,
            'section_heights':self.section_heights
        }

        print()

    def save_QS_parameters_charakteristisch(self):
        '''
        WICHTIG: uinterschied zu save_section_parameters: die Werte an den Knoten werden gespeichert und nicht die section mittelwerte
        Gedacht um die werte die für einen Bestimmten QS charakterisitisch sind zu speichern um sie dann in ein excel zu packen
        TODO Wandstärke bisher konstant
        '''

        self.querschnitts_werte = {
            'Höhe [m]':self.section_absolute_heights,
            'd_achse [m]':self.d_achse,
            'Iz [m^4]': self.Iz,
            'A [m²]':self.A,
            'U_achse [m]':self.U,
            't [m]':np.asarray([self.wand_stärke]*len(self.d_achse))
            
        }

        return self.querschnitts_werte

    def export_object_to_pkl(self, dest_file):

        with open(dest_file, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print ('\nSaved object in', dest_file)

            

class KreisRing(Querschnitt):

    def __init__(self, cd = 1.1, cp_max =1.8 , lagen_aufbau = None, holz_parameter = {}, nachweis_parameter= {}, 
                hoehen_parameter ={}, einheiten = {}, FE_elements=10) -> None:
        
        super().__init__(lagen_aufbau, holz_parameter, nachweis_parameter, hoehen_parameter, einheiten,FE_elements)


        self.masse_pro_meter = None
        self.name = 'Ring'

        self.A_m= self.compute_area(self.d_achse/2)
        self.A = self.compute_area(self.d_außen/2)-self.compute_area(self.d_innen/2)
        self.U = self.d_achse * np.pi
        self.compute_netto_flächen()

        self.nachweis_parameter= nachweis_parameter

        self.masse_pro_meter = self.compute_sectional_mean(self.A) * self.wichte
        self.V = self.compute_sectional_mean(self.A) * self.section_heights
        self.eigengewicht =  - self.V * self.wichte
        self.gewichtskraft = {'Fx':- self.V * self.wichte*GD.GRAVITY}
        self.initalize_geometrie_parameter_FE()

        self.cd = cd
        self.cp_max = cp_max
        
        self.compute_flächenträgheitsmoment()
        self.save_section_parameters()
        self.save_QS_parameters_charakteristisch()


    def compute_flächenträgheitsmoment(self):

        self.Iz = np.pi / 64 * (self.d_außen**4 - self.d_innen**4)

    def compute_area(self, r):
        area= np.pi * r**2

        return area

    def compute_netto_flächen(self):

        self.A_längslagen = 0
        self.A_querlagen = 0

        for i, lage in enumerate(self.lagen_aufbau):
            t_von_innen = sum([li['ti'] for li in self.lagen_aufbau[:i]])
            r_i_innen = self.d_innen/2 + t_von_innen
            r_i_außen = self.d_innen/2 + t_von_innen + lage['ti']
            if lage['ortho'] == 'X':
                self.A_längslagen += self.compute_area(r_i_außen) - self.compute_area(r_i_innen)
            elif lage['ortho'] == 'Y':
                self.A_längslagen += self.compute_area(r_i_außen) - self.compute_area(r_i_innen)

    def compute_static_moment(self, is_effective = False):
        '''
        Sy berechnene mit nur den Längslagen?!
        my Sy bei z=0 (Viertelkreis)
        Berechnung für Kreisring 
        '''
        self.Sy_max= ((self.d_achse/2)**2)*self.wand_stärke 

    def initalize_geometrie_parameter_FE(self):
        '''
        die Querschnittswerte anhand der höhen parameter gleich auch auf die anzahl an FE Elementen hin interpolieren
        '''
        self.x_fe = np.linspace(0,self.section_absolute_heights[-1], self.FE_nodes)
        self.section_heights_FE = np.diff(self.x_fe)
        geometrie_parameter = [self.d_achse, self.A]
        self.d_achse_FE, self.A_FE = [],[]
        geometrie_parameter_FE = [self.d_achse_FE, self.A_FE]

        for param_i, parameter in enumerate(geometrie_parameter):
            for x in self.x_fe:
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
        self.masse_pro_meter = self.compute_sectional_mean(self.A) * self.wichte

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

        self.section_mass = self.masse_pro_meter * self.section_heights
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
                        'M':cross_section.masse_pro_meter}
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
