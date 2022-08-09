from math import cos, radians, sin
import pickle 
import matplotlib.pyplot as plt
import numpy as np
from os.path import join  as os_join
#import utilities.auxiliary as utils


class Querschnitt(object):

    '''
    Klasse geschlossener Querschnitte aus BSP definiert durch den durchmesser der Achse
    '''

    def __init__(self, d_achse, wand_dicke, lagen_aufbau, holz_parameter = {}, hoehen_parameter = {}, einheiten = {}):
        '''
        d_achse: Durchmesser des Turms in der Achse des Querschnitts
        holz_parameter: dict siehe global definitions von einer Bestimmten Holzklasse/Art
        hfract, absolute_height
        '''
        self.d_achse = d_achse
        self.wand_dicke = wand_dicke
        self.lagen_aufbau = lagen_aufbau
        self.wichte = holz_parameter['rhok']

        self.d_außen = self.d_achse + self.wand_dicke/2
        self.d_innen = self.d_achse - self.wand_dicke/2

        self.hfract = hoehen_parameter['hfract']
        self.section_absolute_heights = hoehen_parameter['absolute_höhen'] # Höhen koordinate der section 
        self.section_heights = np.diff(self.section_absolute_heights) # Höhe der einzelnen elemente (abstand zwischen höhen)

        self.holz_parameter = holz_parameter
        self.hoehen_parameter = hoehen_parameter
        self.einheiten = einheiten
        

    def compute_sectional_mean(self, parameter):
        
        if not isinstance(parameter, list):
            parameter = list(parameter)
        sectional_sum = [j+i for i, j in zip(parameter[:-1], parameter[1:])]
        return np.array(sectional_sum)/2

    def compute_effective_moment_of_inertia(self):
        '''
        sollte Überschrieben werden von der child funktion
        '''
        pass
    




    def compute_effektive_festigkeiten_charakteristisch(self):
        
        t_querlagen = 0 
        for lage in self.lagen_aufbau:
            if lage['ortho'] == 'Y':
                t_querlagen+= lage['ti']

        self.ft0k_eff = (self.wand_dicke - t_querlagen)/self.wand_dicke * self.holz_parameter['ft0k']
        self.fc0k_eff = (self.wand_dicke - t_querlagen)/self.wand_dicke * self.holz_parameter['fc0k']
        self.fvk_eff = (self.wand_dicke - t_querlagen)/self.wand_dicke * self.holz_parameter['fvk']

    def compute_effektive_festigkeiten_design(self):

        self.compute_effektive_festigkeiten_charakteristisch()

        #TODO design parameter kmod, gamma_m etc. als parameter dict übergeben?

    
    def calculate_normalspannung(self, lasten_design, e = 0):
        '''
        TODO mit effectiven Iy oder mit dem vom gesamten QS?!
        oder von einzelnen BSP Element?
        e Richtig? 
        '''
        
        #self.compute_effective_moment_of_inertia()

        e = self.d_achse/2

        self.sigma = (lasten_design['Md']+lasten_design['Mimp'])/self.Iy * e + lasten_design['Nd'] / self.A
        self.einheiten['Normalspannung'] = self.einheiten['Kraft'] + '/' + self.einheiten['Länge'] + '²'

    def calculate_ausnutzung(self, lasten_design):
        '''
        TODO Variable geben die definiert was ausgenutzt wird als Druck, Zug, ... bisher ist es druck fc0k
        '''

        self.calculate_normalspannung(lasten_design)
        # self.calculate_schubspannung_querkraft()
        # self.calculate_schubspannung_torsion()

        einwirkung = self.sigma


        fck = self.holz_parameter['fc0k'] * utils.unit_conversion(self.einheiten['Festigkeit'], self.einheiten['Normalspannung'])
        self.ausnutzung = 1.35 * einwirkung/fck


    def calculate_schubspannung_querkraft(self, Qz):

        self.compute_effective_moment_of_inertia()
        self.compute_static_moment()

        tau_Q = Qz * self.sy_max / (self.Iy * self.t)
        self.einheiten['Schubspannung'] = self.einheiten['Kraft'] + '/' + self.einheiten['Länge'] + '²'
        

    def calculate_schubspannung_torsion(self, Mt):
        self.Wt = 2* self.A_m * self.wand_dicke
        tau_Mt = Mt /self.Wt


    def calculate_ausnutzung_kombiniert(self, LF):

        self.get_forces(LF)
        sigma_d_normal = self.calculate_normalspannung(self.Md, self.Mimp, self.Nd)
        tau_d_quer=0
        tau_d_torsion= 0
        

        fck = self.holz_parameter['fc0k'] * utils.unit_conversion('N/mm²', 'N/m²')


        nu = sigma_d_normal/fck
        self.ausnutzung = nu


    def plot_properties_along_height(self, properties:list, units:list=None, E_modul=12000, prop_to_compare=None, title=''):
        '''
        werte entlang der Höhe plotten. 
        gewünschte parameter als strings in einer liste geben
        E_modul = optional wenn EI berechnet werden soll
        'Iy','EIy', 'd_achse', 'a_außen', 'a_innen', 'd_außen','M'
        '''
        prop_values = {'Iy':self.Iy,'EIy':self.Iy * E_modul, 'd_achse':self.d_achse, 'a_außen':self.a_außen, 'a_innen':self.a_innen, 'd_außen':self.d_außen,
                        'M':self.masse_pro_meter}
        prop_units_default = {'Iy':'m^4', 'EIy':'Nm²', 'd_achse':'m', 'a_außen':'m', 'a_innen':'m', 'd_außen':'m', 'M':'kg'}

        fig, ax = plt.subplots(ncols = len(properties))
        if len(properties) == 1:
            ax = [ax]

        for p_i, prop in enumerate(properties):
            ax[p_i].plot(prop_values[prop], self.hfract)
            if prop_to_compare != None:
                ax[p_i].plot(prop_to_compare, self.hfract)
            ax[p_i].grid()
            ax[p_i].set_xlabel(prop + ' [' + prop_units_default[prop] +']')
        ax[0].set_ylabel('Hfract [-]')
        fig.suptitle(title)

        plt.show()

    def export_to_dict_pkl(self, dest_file):

        dict_to_export = {
            'n_sections':len(self.Iy),
            'Iy': self.compute_sectional_mean(self.Iy),
            'Iy_eff': self.compute_sectional_mean(self.Iy_eff),
            'd_achse':self.compute_sectional_mean(self.d_achse),
            'A':self.compute_sectional_mean(self.A),
            'masse_pro_m':self.masse_pro_meter,
            'section_absolute_heights':self.section_absolute_heights,
            'section_heights':self.section_heights
        }

        with open(dest_file, 'wb') as handle:
            pickle.dump(dict_to_export, handle, protocol=pickle.HIGHEST_PROTOCOL)

        print ('\nSaved nEck Daten in', dest_file)


########## Plattenbeanspruchung ##########



#     def calculate_normalspannung_Plattenbeanspruchung(self, wind_max): 
#         '''Längsspannungen bei Plattenbeanspruchung aus Windsog und -Druck'''
#         width= 3
#         M0_k= (width*wind_max**2)/8
#         W0_net= 
#         f_md= 1,35* 
#         sigma_md= M
#         ausnutzung= sigma_md/f_md

#         return ausnutzung

    def calculate_schub_plattenbeanspruchung(self, wind_max):
        return 0

    def compute_drillsteifigkeit(self):
        self.K_xy= np.sqrt(self.Iy_eff*self.holz_parameter['E0_mean']*self.Iz_eff*self.holz_parameter['E90_mean'])
         

class KreisRing(Querschnitt):

    def __init__(self, d_achse, wand_dicke, lagen_aufbau = None, holz_parameter = {}, hoehen_parameter ={}, einheiten = {}) -> None:
        
        super().__init__(d_achse, wand_dicke, lagen_aufbau, holz_parameter, hoehen_parameter, einheiten)


        self.masse_pro_meter = None
        self.hfract = np.linspace(0,1,10)

        self.A_außen = self.compute_area(self.d_außen/2)
        self.A_innen = self.compute_area(self.d_innen/2)
        self.A_m= self.compute_area(self.d_achse/2)

        self.A = self.A_außen-self.A_innen

        self.masse_pro_meter = self.compute_sectional_mean(self.A) * self.wichte
        
        self.compute_flächenträgheitsmoment()
        if self.lagen_aufbau:
            self.compute_effective_moment_of_inertia()


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


        self.Sy_max= ((self.d_achse/2)**2)*self.wand_dicke 
        return self.Sy_max 
      

class nEck(Querschnitt):

    def __init__(self, n_ecken, d_achse, wand_dicke, lagen_aufbau = None, holz_parameter = {}, hoehen_parameter ={}, einheiten = {}):
        '''
        werte können für einzelne sections oder als arrays gegeben werden
        '''

        super().__init__(d_achse, wand_dicke, lagen_aufbau, holz_parameter, hoehen_parameter, einheiten)

        self.n_ecken = n_ecken
        self.alpha = 360/self.n_ecken

        self.a_außen = sin(radians(self.alpha/2)) * self.d_außen
        self.a_innen = sin(radians(self.alpha/2)) * self.d_innen

        self.compute_winkel_term()

        Iy_außen = self.compute_flächenträgheitsmoment_neck(self.a_außen)
        Iy_innen = self.compute_flächenträgheitsmoment_neck(self.a_innen)
        
        self.Iy = Iy_außen - Iy_innen

        if self.lagen_aufbau:
            self.compute_effective_moment_of_inertia()

        self.A_außen = self.compute_area_neck(self.d_außen/2)
        self.A_innen = self.compute_area_neck(self.d_innen/2)
        self.A_m= self.compute_area_neck(self.d_achse/2)

        self.A = self.A_außen-self.A_innen
        self.masse_pro_meter = self.compute_sectional_mean(self.A) * self.wichte


    def compute_winkel_term(self):
        
        self.winkel_term = ((2 + cos(radians(self.alpha))) / (1-cos(radians(self.alpha)))**2 ) * sin(radians(self.alpha))

    def compute_area_neck(self, r_u):
        ''' 
        r_u: Umkreis radius des n-ecks
        '''

        A = self.n_ecken/2 *  r_u**2 * sin(2*np.pi/self.n_ecken)

        return A

    def compute_flächenträgheitsmoment_neck(self, a):

        Iy = self.n_ecken/96 * a**4 * self.winkel_term
        return Iy

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
            Sy_max+= (z_e*length+((z_before-z_e)/2)*length)*self.wand_dicke
            z_before = z_e 
        self.Sy=Sy_max   
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

if __name__ == '__main__':
    import xlwings as xl
    import utilities.holz as holz
    import utilities.utilities as utils
    excel_file = os_join (*['..','..','RFEM','0704_parametrisches_RFEM_model_12eck_save2.xlsx'])
    wb = xl.Book(excel_file)
    ws = wb.sheets['Formel-Parameter']

    ebenen_radien = utils.read_xl_column(ws, start_cell = 'D90', end_row= 104)
    # boden knoten dazu 
    heights_list = utils.read_xl_column(ws, start_cell='D75', end_row=88)
    heights_list.append(0)
    heights_parameter = {}
    # flip -> von unten nach oben
    heights_parameter['absolute_höhen'] = np.flip(np.array(heights_list))
    heights_parameter['hfract'] = np.flip(np.array(utils.read_xl_column(ws, start_cell='J90', end_row=104)))

    # hier soll definiert werden in welchen Einheiten inputs an Querschnitt ÜBERGEBEN wird
    # 
    einheiten = {'Kraft':'N', 'Moment':'Nm', 'Festigkeit':'N/mm²', 'Länge':'m'}

    n_ecken = 12
    ecken_possible = [8,10,12]
    d_achse = np.flip(np.array(ebenen_radien) *2)
    t_wand = 0.4
    
    lagen_aufbau = [{'ortho':'X','ti':0.08, 'di':d_achse + 0.03 + 0.04 + 0.04},
                    {'ortho':'Y','ti':0.04, 'di':d_achse + 0.03 +0.02},
                    {'ortho':'X','ti':0.06, 'di':d_achse},
                    {'ortho':'Y','ti':0.04, 'di':d_achse -  0.03 - 0.02},
                    {'ortho':'X','ti':0.08, 'di':d_achse - 0.03 -  0.04 -  0.04}]

    kreis_ring = KreisRing(d_achse, 0.4, lagen_aufbau=lagen_aufbau,
                           holz_parameter = holz.charakteristische_werte['BSP_RFEM'], hoehen_parameter= heights_parameter, einheiten=einheiten)
    nEck12_fuß = nEck(12, d_achse, t_wand, lagen_aufbau=lagen_aufbau, 
                          holz_parameter = holz.charakteristische_werte['BSP_RFEM'], hoehen_parameter= heights_parameter, einheiten=einheiten)

    print(kreis_ring.compute_static_moment())
    print(nEck12_fuß.compute_static_moment())
