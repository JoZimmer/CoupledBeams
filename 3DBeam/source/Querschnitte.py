from math import cos, radians, sin
import pickle 
import matplotlib.pyplot as plt
import numpy as np
from os.path import join  as os_join
import source.utilities.utilities as utils

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

    def __init__(self, d_achse, lagen_aufbau, holz_parameter = {}, nachweis_parameter= {}, hoehen_parameter = {}, einheiten = {}):
        '''
        d_achse: Durchmesser des Turms in der Achse des Querschnitts
        holz_parameter: dict siehe global definitions von einer Bestimmten Holzklasse/Art
        hfract, absolute_height
        '''
        self.d_achse = d_achse
        self.wand_stärke = 0
        self.lagen_aufbau = lagen_aufbau
        self.t_laengslagen = sum([lage['ti'] for lage in lagen_aufbau if lage['ortho'] == 'X'])
        self.t_querlagen = sum([lage['ti'] for lage in lagen_aufbau if lage['ortho'] == 'Y'])

        for lage in lagen_aufbau:
            self.wand_stärke += lage['ti']
        
        self.wichte = holz_parameter['rhok']

        self.d_außen = self.d_achse + self.wand_stärke/2
        self.d_innen = self.d_achse - self.wand_stärke/2

        self.hfract = hoehen_parameter['hfract']
        self.section_absolute_heights = hoehen_parameter['absolute_höhen'] # Höhen koordinate der section 
        self.section_heights = np.diff(self.section_absolute_heights) # Höhe der einzelnen elemente (abstand zwischen höhen)
        self.max_height = self.section_absolute_heights[-1]

        self.holz_parameter = holz_parameter
        self.nachweis_parameter= nachweis_parameter
        self.hoehen_parameter = hoehen_parameter
        self.einheiten = einheiten      

        self.initialize_eimheiten_umrechnung()


    def initialize_eimheiten_umrechnung(self):
        self.einheiten_umrechnung = {'Normalspannung':1}

        einheit_normalspannung = self.einheiten['Kraft'] + '/' + self.einheiten['Länge'] + '²'

        if einheit_normalspannung != self.einheiten['Festigkeit']:
            self.einheiten_umrechnung['Normalspannung'] =  utils.unit_conversion(einheit_normalspannung, self.einheiten['Festigkeit'])
            self.einheiten['Normalspannung'] = self.einheiten['Festigkeit']
        
        else:
            self.einheiten['Normalspannung'] = einheit_normalspannung


    def compute_sectional_mean(self, parameter):
        
        if not isinstance(parameter, list):
            parameter = list(parameter)
        sectional_sum = [j+i for i, j in zip(parameter[:-1], parameter[1:])]
        return np.array(sectional_sum)/2

# _________ EFFEKTIVE TRÄGHEITSMOMENTE _________________________ 

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
                print(lage['ti'])
        self.Iy_eff_platte= eigenanteil+steineranteil
        print(self.Iy_eff_platte)

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
                print(lage['ti'])
        self.Iy_eff_platte= eigenanteil+steineranteil
        print('Iz')
        print(self.Iy_eff_platte)

    def compute_effective_moment_of_inertia(self):
        '''
        sollte Überschrieben werden von der child funktion
        '''
        pass
    
# ________ FESTIGKEITEN/ WIDERTÄNDE ____________________________

    def compute_effektive_festigkeiten_charakteristisch(self):
        #effektive festigkeiten laengs 
        t_querlagen = 0 
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'Y':
                t_querlagen+= lage['ti']

        self.fmk_eff_laengs = (self.wand_stärke - t_querlagen)/self.wand_stärke * self.holz_parameter['fmk']
        self.ft0k_eff_laengs = (self.wand_stärke - t_querlagen)/self.wand_stärke * self.holz_parameter['ft0k']
        self.fc0k_eff_laengs = (self.wand_stärke - t_querlagen)/self.wand_stärke * self.holz_parameter['fc0k']
        self.fvk_eff_laengs = (self.wand_stärke - t_querlagen)/self.wand_stärke * self.holz_parameter['fvk']
    
        #effektive Festigkeiten quer
        t_laengslagen = 0 
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'X':
                t_laengslagen+= lage['ti']

        self.fmk_eff_quer = (self.wand_stärke - t_laengslagen)/self.wand_stärke * self.holz_parameter['fmk']
        self.ft0k_eff_quer = (self.wand_stärke - t_laengslagen)/self.wand_stärke * self.holz_parameter['ft0k']
        self.fc0k_eff_quer = (self.wand_stärke - t_laengslagen)/self.wand_stärke * self.holz_parameter['fc0k']
        self.fvk_eff_quer = (self.wand_stärke - t_laengslagen)/self.wand_stärke * self.holz_parameter['fvk']

    def compute_effektive_festigkeiten_design(self, einwirkungsdauer):


        self.compute_effektive_festigkeiten_charakteristisch()

        #effektive Festigkeite laengs
        self.fmd_eff_laengs = self.fmk_eff_laengs * self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']
        self.ft0d_eff_laengs= self.ft0k_eff_laengs * self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']
        self.fc0d_eff_laengs= self.fc0k_eff_laengs * self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']
        self.fvd_eff_laengs= self.fvk_eff_laengs * self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']
    
        #effektive Festigkeiten Quer
        self.fmd_eff_quer = self.fmk_eff_quer * self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']
        self.ft0d_eff_quer= self.ft0k_eff_quer * self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']
        self.fc0d_eff_quer= self.fc0k_eff_quer * self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']
        self.fvd_eff_quer= self.fvk_eff_quer * self.nachweis_parameter['k_mod'][einwirkungsdauer]*(1/self.nachweis_parameter['gamma_m'])* self.nachweis_parameter['k_sys']
        
        
# _________ NORMALSPANNUNGEN __________________________________
    
    def calculate_normalspannung(self, lasten_design, add_vorspannkraft = False):
        '''
        ergibt dictionaries mit der einwirkungsdauer als key:
            sigma_druck: Spannung aus Normalkraft und moment
            sigma_zug: Spannung aus moment abzüglich Normalkraft(Annahme immer Druckkraft)
            sigma_N: Spannung nur mit Normalkraft berechnet
            sigma_M: Spannung nur mit Moment berechnet
        TODO e Richtig oder sollte das der Außenradius sein?!
        '''
        
        e = self.d_außen/2

        self.sigma_druck, self.sigma_zug, self.sigma_N, self.sigma_M = {},{},{},{}

        for dauer in lasten_design:
            # Wirkungsrichtung von M neutralisieren

            lasten_design[dauer]['Mz'] = abs(lasten_design[dauer]['Mz'])

            self.sigma_druck[dauer] = -(lasten_design[dauer]['Mz'])/self.Iy * e + lasten_design[dauer]['Nx'] / self.A
            self.sigma_zug[dauer] = (lasten_design[dauer]['Mz'])/self.Iy * e + lasten_design[dauer]['Nx'] / self.A

            if add_vorspannkraft:
                self.sigma_druck[dauer] += np.negative(self.sigma_zug[dauer])

            self.sigma_N[dauer] = lasten_design[dauer]['Nx'] / self.A
            self.sigma_M[dauer] = (lasten_design[dauer]['Mz'])/self.Iy * e 

        

    def calculate_ausnutzung_normalspannung(self, lasten_design, add_vorspannkraft = False):
        '''
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

        self.calculate_normalspannung(lasten_design, add_vorspannkraft)

        # erst mal Einheiten klar machen
        einheiten_faktor = self.einheiten_umrechnung['Normalspannung']
        self.sigma_druck_design = self.sigma_druck['egal'] * einheiten_faktor
        self.sigma_zug_design = self.sigma_zug['egal'] * einheiten_faktor
        self.sigma_N_design = self.sigma_N['egal'] * einheiten_faktor
        self.sigma_M_design = self.sigma_M['egal'] * einheiten_faktor

        self.ausnutzung_druck, self.ausnutzung_zug, self.ausnutzung_N, self.ausnutzung_M = 0,0,0,0
        for dauer in lasten_design:
            if dauer == 'egal':
                continue

            self.compute_effektive_festigkeiten_design(dauer)

            self.ausnutzung_druck += abs(self.sigma_druck[dauer]* einheiten_faktor)/self.fc0d_eff_laengs 
            self.ausnutzung_zug += abs(self.sigma_zug[dauer]* einheiten_faktor)/self.ft0d_eff_laengs 
            self.ausnutzung_N += abs(self.sigma_N[dauer]* einheiten_faktor)/self.fc0d_eff_laengs # Annahme nur Normalkraft
            self.ausnutzung_M += abs(self.sigma_M[dauer]* einheiten_faktor)/self.fc0d_eff_laengs 

    
# _________ SCHUBSPANNUNGEN __________________________________

    def calculate_schubspannung_querkraft(self, lasten_design):

        self.compute_effective_moment_of_inertia()
        self.compute_static_moment()

        self.tau_Q = (lasten_design['Q'] * self.Sy_max) / (self.Iy * self.wand_stärke)
        self.einheiten['Schubspannung'] = self.einheiten['Kraft'] + '/' + self.einheiten['Länge'] + '²'
        
    def calculate_schubspannung_torsion(self, lasten_design):
        self.Wt = 2* self.A_m * self.wand_stärke
        self.tau_Mt = lasten_design['Mt'] /self.Wt

    def calculate_ausnutzung_schub(self, lasten_design):
        '''
        TODO Variable geben die definiert was ausgenutzt wird als Druck, Zug, ... bisher ist es druck fc0k
        '''

        self.calculate_schubspannung_querkraft(lasten_design)
        self.calculate_schubspannung_torsion(lasten_design)
        self.compute_effektive_festigkeiten_design('kurz')

        einwirkung = self.tau_Q + self.tau_Mt 


        fd = self.fvd_eff_laengs * utils.unit_conversion(self.einheiten['Festigkeit'], self.einheiten['Normalspannung'])
        self.ausnutzung = 1.35 * einwirkung/fd

        print(self.ausnutzung)

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
        print('ausnutzung_kombiniert')
        print(self.ausnutzunG)


# _________ PLATTENBEANSPRUCHUNG __________________________________


    def calculate_normalspannung_Plattenbeanspruchung_Nebentragrichtung(self, wind_max, width): 
        '''Längsspannungen bei Plattenbeanspruchung aus Windsog und -Druck'''
        
        self.compute_effektive_moment_of_inertia_platte_z()
        M0_k= (wind_max*width**2)/8
        self.sigma_md_platte= M0_k/self.Iy_eff_platte * self.wand_stärke/2

        return self.sigma_md_platte

    def ausnutzung_Plattenbeanspruchung_Nebentragrichtung(self, einwirkungsdauer, wind_max, width):
        
        self.calculate_normalspannung_Plattenbeanspruchung_Nebentragrichtung(wind_max, width)
        self.compute_effektive_festigkeiten_design(einwirkungsdauer)
        fd= self.fmd_eff_quer * utils.unit_conversion(self.einheiten['Festigkeit'], self.einheiten['Normalspannung'])
        ausnutzung= self.sigma_md_platte/fd
        print('ausnutzung Platte quer')
        print(ausnutzung)
        return  ausnutzung

    def calculate_normalspannung_Plattenbeanspruchung_Haupttragrichtung():
        return 0

    def calculate_schub_plattenbeanspruchung(self, wind_max):
        return 0

    def calculate_schub_plattenbeanspruchung_rollschub(self):
        '''Rollschubfestigkeit der schwerpunktnächsten Querlage maßgebend'''
        A_t_net=0 
        I0_net=0
        S0R_net=0

    def compute_drillsteifigkeit(self):
        self.K_xy= np.sqrt(self.Iy_eff_platte*self.holz_parameter['E0_mean']*self.Iz_eff_platte*self.holz_parameter['E90_mean'])

# _________ OBJEKT FUNKTIONEN __________________________________

    def plot_properties_along_height(self, properties:list, units:list=None, E_modul=12000, prop_to_compare=None):
        '''
        werte entlang der Höhe plotten. 
        gewünschte parameter als strings in einer liste geben
        E_modul = optional wenn EI berechnet werden soll
        'Iy','EIy', 'd_achse', 'a_außen', 'a_innen', 'd_außen','M'
        '''
        prop_values = {'Iy':self.Iy,'EIy':self.Iy * E_modul, 'd_achse':self.d_achse, 'a_außen':self.a_außen, 'a_innen':self.a_innen, 'd_außen':self.d_außen,
                        'M':self.masse_pro_meter}
        prop_units_default = {'Iy':'m^4', 'EIy':'Nm²', 'd_achse':'m', 'a_außen':'m', 'a_innen':'m', 'd_außen':'m', 'M':'kg'}

        fg, ax = plt.subplots(ncols = len(properties))
        if len(properties) == 1:
            ax = [ax]

    
        for p_i, prop in enumerate(properties):
            ax[p_i].plot(prop_values[prop], self.hfract)
            if prop_to_compare != None:
                ax[p_i].plot(prop_to_compare, self.hfract)
            ax[p_i].grid()
            ax[p_i].set_xlabel(prop + ' [' + prop_units_default[prop] +']')
        ax[0].set_ylabel('Hfract [-]')

        plt.show()

    def save_section_parameters(self):

        self.section_parameters = {
            'n_sections':len(self.Iy),
            'Iy': self.compute_sectional_mean(self.Iy),
            #'Iy_eff': self.compute_sectional_mean(self.Iy_eff),
            'd_achse':self.compute_sectional_mean(self.d_achse),
            'A':self.compute_sectional_mean(self.A),
            'masse_pro_m':self.masse_pro_meter,
            'section_absolute_heights':self.section_absolute_heights,
            'section_heights':self.section_heights
        }

        # with open(dest_file, 'wb') as handle:
        #     pickle.dump(self.section_parameters, handle, protocol=pickle.HIGHEST_PROTOCOL)

        # print ('\nSaved nEck Daten in', dest_file)

    def export_object_to_pkl(self, dest_file):

        with open(dest_file, 'wb') as handle:
            pickle.dump(self, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print ('\nSaved nEck object in', dest_file)

            

class KreisRing(Querschnitt):

    def __init__(self, d_achse, cd = 1.1, lagen_aufbau = None, holz_parameter = {}, nachweis_parameter= {}, hoehen_parameter ={}, einheiten = {}) -> None:
        
        super().__init__(d_achse, lagen_aufbau, holz_parameter, nachweis_parameter, hoehen_parameter, einheiten)


        self.masse_pro_meter = None
        self.hfract = np.linspace(0,1,10)
        self.name = 'Ring'

        self.A_außen = self.compute_area(self.d_außen/2)
        self.A_innen = self.compute_area(self.d_innen/2)
        self.A_m= self.compute_area(self.d_achse/2)

        self.A = self.A_außen-self.A_innen
        self.nachweis_parameter= nachweis_parameter
        self.masse_pro_meter = self.compute_sectional_mean(self.A) * self.wichte

        self.cd = cd
        
        self.compute_flächenträgheitsmoment()
        self.save_section_parameters()
        # if self.lagen_aufbau:
        #     self.compute_effective_moment_of_inertia()


    def compute_flächenträgheitsmoment(self):

        self.Iy = np.pi / 64 * (self.d_außen**4 - self.d_innen**4)

    def compute_effective_moment_of_inertia(self):
        '''
        Flächenträgheitsmoment der Längslagen
        '''
        for i, lage in enumerate(self.lagen_aufbau):
            querschnitt = KreisRing(lage['di'], lage['ti'], 
                    holz_parameter = self.holz_parameter, hoehen_parameter= self.hoehen_parameter, einheiten =self.einheiten)
            lage['Iy'] = querschnitt.Iy

        self.Iy_eff = 0
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'X':
                self.Iy_eff += lage['Iy']

    def compute_area(self, r):
        area= np.pi * r**2

        return area

    def compute_static_moment(self, is_effective = False):
        '''
        Sy berechnene mit nur den Längslagen?!
        my Sy bei z=0 (Viertelkreis)
        Berechnung für Kreisring 
        '''


        self.Sy_max= ((self.d_achse/2)**2)*self.wand_stärke 
        return self.Sy_max 


        

class nEck(Querschnitt):

    def __init__(self, n_ecken, d_achse, cd = 1.5, lagen_aufbau = None, holz_parameter = {}, nachweis_parameter= {}, hoehen_parameter ={}, einheiten = {}):
        '''
        werte können für einzelne sections oder als arrays gegeben werden
        Geometrie: https://de.wikipedia.org/wiki/Regelm%C3%A4%C3%9Figes_Polygon  
        '''

        super().__init__(d_achse, lagen_aufbau, holz_parameter, nachweis_parameter, hoehen_parameter, einheiten)

        self.n_ecken = n_ecken
        self.name = str(self.n_ecken) + '-Eck'
        self.alpha = 360/self.n_ecken
        self.nachweis_parameter= nachweis_parameter
        self.a_außen = sin(radians(self.alpha/2)) * self.d_außen
        self.a_innen = sin(radians(self.alpha/2)) * self.d_innen

        self.compute_winkel_term()

        Iy_außen = self.compute_flächenträgheitsmoment_neck(self.a_außen)
        Iy_innen = self.compute_flächenträgheitsmoment_neck(self.a_innen)
        
        self.Iy = Iy_außen - Iy_innen

        # if self.lagen_aufbau:
        #     self.compute_effective_moment_of_inertia()

        self.A_außen = self.compute_area_neck(self.d_außen/2)
        self.A_innen = self.compute_area_neck(self.d_innen/2)
        self.A_m= self.compute_area_neck(self.d_achse/2)

        self.A = self.A_außen-self.A_innen
        self.masse_pro_meter = self.compute_sectional_mean(self.A) * self.wichte

        self.cd = cd
        self.save_section_parameters()


    def compute_winkel_term(self):
        
        self.winkel_term = ((2 + cos(radians(self.alpha))) / (1-cos(radians(self.alpha)))**2 ) * sin(radians(self.alpha))

    def compute_area_neck(self, r_u):
        ''' 
        r_u: Umkreis radius des n-ecks
        '''

        A = self.n_ecken/2 *  r_u**2 * sin(2*np.pi/self.n_ecken)

        return A

    def compute_flächenträgheitsmoment_neck(self, a):

        self.Iy = self.n_ecken/96 * a**4 * self.winkel_term
        return self.Iy

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
            lage['Iy'] = nEck_section.Iy

        self.Iy_eff = 0
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'X':
                self.Iy_eff += lage['Iy']

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

    factors = [geometric_objects[-1].Iy_eff/geometric_objects[0].Iy_eff, geometric_objects[-1].Iy_eff/geometric_objects[1].Iy_eff, geometric_objects[1].Iy_eff/geometric_objects[0].Iy_eff]
    factors_label = ['EI12/EI8','EI12/EI10','EI10/EI8']

    for f_i , factor in enumerate(factors):
        print (factors_label[f_i] + ':', round(factor[0],2))

    for e_i, cross_section in enumerate(geometric_objects):
        prop_values = {'Iy':cross_section.Iy,'EIy':cross_section.Iy * E_modul, 'Iy_eff':cross_section.Iy_eff,'EIy_eff':cross_section.Iy_eff * E_modul,
                       'a_außen':cross_section.a_außen, 'a_innen':cross_section.a_innen, 'd_außen':cross_section.d_außen,
                        'M':cross_section.masse_pro_meter}
        prop_units_default = {'Iy':'m^4', 'EIy':'Nm²', 'Iy_eff':'m^4', 'EIy_eff':'Nm²', 'a_außen':'m', 'a_innen':'m', 'd_außen':'m', 'M':'kg'}

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
    
    #f = round(geometric_objects[1].Iy[0]/geometric_objects[0].Iy[0],2)
    #ax[0].plot(0,0, label = 'factor '+ str(f))

    plt.grid()
    plt.legend()
    plt.show()
