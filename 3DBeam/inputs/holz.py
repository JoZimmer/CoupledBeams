
'''
TODO alles in eine Klasse Holz überführen
    - Attribute:
        Lagenaufbau
        chrakteristische Werte
        Holzbau sachen 
        Einheiten
    - Lagenstärken anhand der gesamt dicke automatische generieren mit ein paar vorzugsstärken der Bretter
'''

HOLZBAU = {
    'k_mod':{'kurz': 0.9, 'mittel':0.8, 'lang': 0.7, 'ständig':0.6},
    # kdef von Martin Tabelle 3.7
    'k_def':{'NKL1':{
                      30: 0.2, # Belastungsgrad > 30%
                      15:0.1,   # Belastungsgrad < 15%
                      'sonst':0.15
            },
            'NKL2':{
                      30: 0.4, # Belastungsgrad > 30%
                      15:0.3,   # Belastungsgrad < 15%
                      'sonst':0.2
            }
    },
    'gamma_m': 1.25, 
    'k_sys': 1.0,
    'mu':0.5,  # Reibungsbeiwert
}
#
# def calculate_k_sys(Breite_BSP):
#     '''k_sys: Systembeiwert zur erhöhung des Bauteilwiderstands wegen statistischer Effekte gegenüber Brettfestigkeit'''
#     k_sys= 1.08+0.2*Breite_BSP
#     if k_sys >1.2:
#         k_sys= 1.2
#     return k_sys



charakteristische_werte = {
    'units':{'Festigkeit':['N/mm²'], 'Steifigkeit':['N/m²'], 'Rohdichte':['kg/m³'], 'Temp. Dehnung':['K^-1']},
    'C24':{
        'fmk':24,
        'ft0k':14,
        'ft90k':0.4,
        'fc0k':21,
        'fc90k':2.5,
        'fvk':4.0,
        'fvxyk':5.5, # Standartwert 
        'ftornodek':2.5, # Standartwerd
        'fvFek':31, # Scherfestigkeit einer 45x45° Furnierebene Buche
        'rhok':460,
        'E0mean':12000E+06,
        'E90mean':370E+06,
        'Furnierebene':True,
        'rhok_Fe':735, # aus Lechner Baubuche Furnier
        'G0mean':500, # zwischen Lechner BSH 586 und ETA Radiusholz 690
        'G4545':4227,# Aus Dissertation ermittelte Schubsteifigkeit der Furnierebene
        'mu':0.5, #Reibuungskoeffizient
        'alpha':0.02, # schwindmaß je % Holzfuechteänderung
        'alphaT':4E-06, # K^-1 Temperaturausdehnungskoeffizeint
    },
    'C30':{
        'fmk':30,
        'ft0k':18,
        'ft90k':0.4,
        'fc0k':23,
        'fc90k':2.7,
        'fvk':4.0,
        'fvxyk':5.5,  # Standartwert 
        'ftornodek':2.5, 
        'fvFek':31, # Scherfestigkeit einer 45x45° Furnierebene Buche
        'rhok':460,
        'E0mean':12000E+06,
        'E90mean':370E+06,
        'Furnierebene':True,
        'rhok_Fe':735, # aus Lechner Baubuche Furnier
        'G0mean':500, # zwischen Lechner BSH 586 und ETA Radiusholz 690
        'G4545':4227,# Aus Dissertation ermittelte Schubsteifigkeit der Furnierebene
        'mu':0.5, #Reibuungskoeffizient
        'alpha':0.02, # schwindmaß je % Holzfuechteänderung
        'alphaT':4E-06, # K^-1 Temperaturausdehnungskoeffizeint       
    }
}

# Key gibt die gesamt dicke an
# 'a': ist hier nicht aktuell 
class Lagenaufbauten():

    def __init__(self, params:dict=None, typ:str='', gesamtstärke:int=20, t_querlagen:float = 0.04, t_furnier:float=0.009):
        '''
        typ: 'BSP_furnier' oder 'BSP_normal'
        nur 5 schichtige bisher
        gesamtstärke: hier sind ein paar hinterlegt bisher. einzelne Lagenstärken sind eher zufällig gesetzt
        t_querlagen: wenn mit querlagen dann 2 einheitliche Querlagen
        t_furnier: stärker der furnier ebene auch wenn diese aus mehreren furnieren zusammengesetzt ist
        '''
        if params:
            self.tges = params['dicken'][0]
            self.t_querlagen = params['t_querlagen']
            self.t_furnier = params['t_furnier']
        else:
            self.typ = typ # ob 
            self.tges = gesamtstärke
            self.t_querlagen = t_querlagen
            self.t_furnier = t_furnier # 2 Furnierlagen a 4.5 mm -> eine Furnierebene also 9 mm

        if params['mit_furnier']:
            self.lagenaufbauten = {
                    13:[{'ortho':'X','ti':0.04,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.03},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.04,}],

                    15:[{'ortho':'X','ti':0.05,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.03},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.05,}],

                    18:[{'ortho':'X','ti':0.06,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.04},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.06,}],

                    20:[{'ortho':'X','ti':0.07,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.04},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.07,}],

                    22:[{'ortho':'X','ti':0.08,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.04},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.08,}],

                    26:[{'ortho':'X','ti':0.08,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.08},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.08,}],

                    30:[{'ortho':'X','ti':0.10,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.08},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.10,}],

                    38:[{'ortho':'X','ti':0.12,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.12},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.12,}],

                    40:[{'ortho':'X','ti':0.14,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.10},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.14,}],

                    44:[{'ortho':'X','ti':0.16,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.10},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.16,}],

                    48:[{'ortho':'X','ti':0.18,},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.10},
                    {'ortho':'Y','ti':self.t_furnier,},
                    {'ortho':'X','ti':0.18,}],

                }
            
        else:
            self.lagenaufbauten = {
            36:[{'ortho':'X','ti':0.10,},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.08},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.10,}],

            40:[{'ortho':'X','ti':0.12,},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.08},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.12,}],

            44:[{'ortho':'X','ti':0.14,},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.08},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.14,}],

            48:[{'ortho':'X','ti':0.16,},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.08},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.16,}],

            56:[{'ortho':'X','ti':0.20,},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.08},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.20,}],

            64:[{'ortho':'X','ti':0.20,},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.16},
                {'ortho':'Y','ti':self.t_querlagen,},
                {'ortho':'X','ti':0.20,}],
            }
        
    def get(self):
        return self.lagenaufbauten[self.tges]
